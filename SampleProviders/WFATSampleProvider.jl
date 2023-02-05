using ..Targets
using HDF5
using Rotations
using WignerD
using Base.Threads

"Sample provider which yields electron samples through WFAT formula, matching `IonRateMethod=:WFAT`"
struct WFATSampleProvider <: ElectronSampleProvider
    laser           ::Laser;
    target          ::Molecule;     # WFAT only supports [Molecule]
    tSamples        ::AbstractVector;
    ss_pdSamples    ::AbstractVector;
    ss_pzSamples    ::AbstractVector;
    ionRatePrefix   ::Symbol;       # currently supports :ExpRate.
    tunExit         ::Symbol;       # :Para for tunneling, :IpF for over-barrier, automatically specified.
    mol_Ip          ::Real;
    mol_μ           ::Vector;
    wfat_intdata    ::Array;
    wfat_nξMax      ::Int;
    wfat_mMax       ::Int;
    wfat_lMax       ::Int;

    function WFATSampleProvider(;
                                laser               ::Laser,
                                target              ::Molecule,
                                sample_tSpan        ::Tuple{<:Real,<:Real},
                                sample_tSampleNum   ::Int,
                                rate_ionRatePrefix  ::Symbol,
                                ss_pdMax            ::Real,
                                ss_pdNum            ::Int,
                                ss_pzMax            ::Real,
                                ss_pzNum            ::Int,
                                mol_ionOrbitRelHOMO ::Int,
                                kwargs...   # kwargs are surplus params.
                                )
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[WFATSampleProvider] Invalid time sample number $sample_tSampleNum."
        @assert (ss_pdNum>0 && ss_pzNum>0) "[WFATSampleProvider] Invalid pd/pz sample number $ss_pdNum/$ss_pzNum."
        # load WFAT IntData.
        if ! (mol_ionOrbitRelHOMO in MolWFATAvailableIndices(target))
            MolCalcWFATData!(target, mol_ionOrbitRelHOMO)
        end
        μ, intdata = MolWFATData(target, mol_ionOrbitRelHOMO)
        Ip = IonPotential(target, mol_ionOrbitRelHOMO)
        nξMax = size(intdata,1) - 1
        mMax = round(Int,(size(intdata,2)-1)/2)
        lMax = size(intdata,3) - 1
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate])
            error("[WFATSampleProvider] Undefined tunneling rate prefix [$rate_ionRatePrefix].")
        end
        # check Keldysh parameter & over-barrier condition.
        F0 = LaserF0(laser)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        if γ0 ≥ 0.5
            @warn "[WFATSampleProvider] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[WFATSampleProvider] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        F_crit = Ip^2/4/(1-sqrt(Ip/2))
        tunExit = :Para
        if F0 ≥ F_crit*2/3
            @warn "[WFATSampleProvider] Peak electric field strength F0=$F0, reaching 2/3 of over-barrier critical value, weak-field condition not sufficiently satisfied."
        elseif F0 ≥ F_crit
            @warn "[WFATSampleProvider] Peak electric field strength F0=$F0, reaching the over-barrier critical value, weak-field condition unsatisfied."
            tunExit = :IpF
        end
        # finish initialization.
        return new( laser, target,
                    range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                    range(-abs(ss_pdMax),abs(ss_pdMax);length=ss_pdNum), range(-abs(ss_pzMax),abs(ss_pzMax);length=ss_pzNum),
                    rate_ionRatePrefix, tunExit,
                    Ip, μ, # Ionizing Orbit Energy (-Ip), orbital dipole moment μ
                    intdata, nξMax, mMax, lMax
                    )
    end
end

"Gets the total number of batches."
function batchNum(sp::WFATSampleProvider)
    return length(sp.tSamples)
end

"Generates a batch of electrons of `batchId` from `sp` using WFAT method."
function generateElectronBatch(sp::WFATSampleProvider, batchId::Int)
    t = sp.tSamples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft = hypot(Fxt,Fyt)
    φ_field = atan( Fyt, Fxt)   # direction of field vector F.
    φ_exit  = atan(-Fyt,-Fxt)   # direction of tunneling exit, which is opposite to F.
    Ip = sp.mol_Ip
    κ  = sqrt(2Ip)
    μ  = sp.mol_μ
    Z  = MolCharge(sp.target) + 1
    nξMax = sp.wfat_nξMax
    mMax  = sp.wfat_mMax
    if Ft == 0
        return nothing
    end

    # determining tunneling exit position (using ADK's parabolic tunneling exit method)
    r_exit = (Ip + sqrt(Ip^2 - 4*(1-sqrt(Ip/2))*Ft)) / 2Ft

    # determining Euler angles (β,γ) (α contributes nothing to Γ, thus is neglected)
    mα,mβ,mγ = MolRotation(sp.target)
    RotMol = RotZYZ(mγ, mβ, mα)  # the molecule's rotation
    RotLaser = RotMatrix3([[ 0  -sin(φ_field)  cos(φ_field)];
                           [ 0   cos(φ_field)  sin(φ_field)];
                           [-1              0             0]])          # the laser's rotation (directions of x&y can be arbitrary)
    α,β,γ = Rotations.params(RotZYZ(inv(RotLaser)*RotMol))      # the ZYZ Euler angles of the rotations from laser frame (F points to Z axis) to molecular frame.

    # determining tunneling rate Γ.
    # The total rate Γ consists of partial rates of different channels ν=(nξ,m): Γ = ∑ Γ_ν
    # The partial rate consists of structure factor part |G_ν(β,γ)|² and field factor W_ν(F): Γ_ν = |G_ν(β,γ)|²*W_ν(F)
    μ = RotMol * μ  # μ is initially obtained in the MF.
    μ_field = μ[1]*cos(φ_field)+μ[2]*sin(φ_field)   # μ's component along the field F, aka μ_z (here F is not along the Z axis!).
    g_data = zeros(nξMax+1, 2*mMax+1)
    for nξ in 0:nξMax, m in -mMax:mMax
        g_data[nξ+1, m+mMax+1] = _getMolStructFactor_g(sp, nξ, m, β, γ)
    end
    g(nξ,m) = g_data[nξ+1, m+mMax+1]
    ionRate::Function =
        if sp.ionRatePrefix == :ExpRate
            function (F,pd,pz)
                Γsum = 0.
                for nξ in 0:nξMax, m in -mMax:mMax
                    G2 = ( exp(-κ*μ_field)*g(nξ,m) )^2  # G², structural part
                    WF = (κ/2) * (4κ^2/F)^(2Z/κ-2nξ-abs(m)-1) * exp(-2(κ^2+pd^2+pz^2)^1.5/3F)   # W_F, field part
                    Γsum += G2*WF
                end
                return Γsum
            end
        else
            #TODO: Add support for full prefixes.
        end
    dim = 8
    pdNum, pzNum = length(sp.ss_pdSamples), length(sp.ss_pzSamples)
    init = zeros(Float64, dim, pdNum, pzNum) # initial condition
    x0 = r_exit*cos(φ_exit)
    y0 = r_exit*sin(φ_exit)
    z0 = 0.
    @threads for ipd in 1:pdNum
        pd0 = sp.ss_pdSamples[ipd]
        px0 = pd0*-sin(φ_exit)
        py0 = pd0* cos(φ_exit)
        for ipz in 1:pzNum
            pz0 = sp.ss_pzSamples[ipz]
            init[1:8,ipd,ipz] = [x0,y0,z0,px0,py0,pz0,t,ionRate(Ft,pd0,pz0)]
        end
    end
    return reshape(init,dim,:)
end

"Gets the structure factor g of the molecule target of a given channel under a specific Euler angle (β,γ)."
function _getMolStructFactor_g(sp::WFATSampleProvider, nξ::Int, m::Int, β, γ)
    @assert nξ≥0 "[WFATSampleProvider] The nξ must be positive."
    sum = zero(ComplexF64)
    for l in abs(m):sp.wfat_lMax, m_ in -l:l
        sum += sp.wfat_intdata[nξ+1,m+sp.wfat_mMax+1,l+1,m_+l+1] * WignerD.wignerdjmn(l,m,m_,β) * exp(-1im*m_*γ)
    end
    return real(sum)
end
