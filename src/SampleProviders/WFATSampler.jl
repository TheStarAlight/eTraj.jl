using ..Targets
using HDF5
using Rotations
using WignerD
using Base.Threads

"Sample provider which yields electron samples through WFAT formula, matching `IonRateMethod=:WFAT`"
struct WFATSampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::Molecule;     # WFAT only supports [Molecule]
    tSamples        ::AbstractVector;
    ss_kdSamples    ::AbstractVector;
    ss_kzSamples    ::AbstractVector;
    ionRatePrefix   ::Symbol;       # currently supports :ExpRate.
    tunExit         ::Symbol;       # :Para for tunneling, :IpF for over-barrier, automatically specified.
    mol_Ip          ::Real;
    mol_μ           ::Vector;
    wfat_intdata    ::Array;
    wfat_nξMax      ::Int;
    wfat_mMax       ::Int;
    wfat_lMax       ::Int;

    function WFATSampler(;  laser               ::Laser,
                            target              ::Molecule,
                            sample_tSpan        ::Tuple{<:Real,<:Real},
                            sample_tSampleNum   ::Int,
                            rate_ionRatePrefix  ::Symbol,
                            ss_kdMax            ::Real,
                            ss_kdNum            ::Int,
                            ss_kzMax            ::Real,
                            ss_kzNum            ::Int,
                            mol_ionOrbitRelHOMO ::Int,
                            kwargs...   # kwargs are surplus params.
                            )
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[WFATSampler] Invalid time sample number $sample_tSampleNum."
        @assert (ss_kdNum>0 && ss_kzNum>0) "[WFATSampler] Invalid kd/kz sample number $ss_kdNum/$ss_kzNum."
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
            error("[WFATSampler] Undefined tunneling rate prefix [$rate_ionRatePrefix].")
        end
        # check Keldysh parameter & over-barrier condition.
        F0 = LaserF0(laser)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        if γ0 ≥ 0.5
            @warn "[WFATSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[WFATSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        F_crit = Ip^2/4/(1-sqrt(Ip/2))
        tunExit = :Para
        if F0 ≥ F_crit
            @warn "[WFATSampler] Peak electric field strength F0=$F0, reaching the over-barrier critical value, weak-field condition unsatisfied. Tunneling exit method switched from [Para] to [IpF]."
            tunExit = :IpF
        elseif F0 ≥ F_crit*2/3
            @warn "[WFATSampler] Peak electric field strength F0=$F0, reaching 2/3 of over-barrier critical value, weak-field condition not sufficiently satisfied."
        end
        # finish initialization.
        return new( laser, target,
                    range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                    range(-abs(ss_kdMax),abs(ss_kdMax);length=ss_kdNum), range(-abs(ss_kzMax),abs(ss_kzMax);length=ss_kzNum),
                    rate_ionRatePrefix, tunExit,
                    Ip, μ, # Ionizing Orbit Energy (-Ip), orbital dipole moment μ
                    intdata, nξMax, mMax, lMax
                    )
    end
end

"Gets the total number of batches."
function batchNum(sp::WFATSampler)
    return length(sp.tSamples)
end

"Generates a batch of electrons of `batchId` from `sp` using WFAT method."
function generateElectronBatch(sp::WFATSampler, batchId::Int)
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

    # determining tunneling exit position (using ADK's parabolic tunneling exit method if tunExit=:Para)
    r_exit = if sp.tunExit == :Para
        (Ip + sqrt(Ip^2 - 4*(1-sqrt(Ip/2))*Ft)) / 2Ft
    else
        Ip / Ft
    end

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
            function (F,kd,kz)
                Γsum = 0.
                for nξ in 0:nξMax, m in -mMax:mMax
                    G2 = ( exp(-κ*μ_field)*g(nξ,m) )^2  # G², structural part
                    WF = (κ/2) * (4κ^2/F)^(2Z/κ-2nξ-abs(m)-1) * exp(-2(κ^2+kd^2+kz^2)^1.5/3F)   # W_F, field part
                    Γsum += G2*WF
                end
                return Γsum
            end
        else
            #TODO: Add support for full prefixes.
        end
    dim = 8
    kdNum, kzNum = length(sp.ss_kdSamples), length(sp.ss_kzSamples)
    init = zeros(Float64, dim, kdNum, kzNum) # initial condition
    x0 = r_exit*cos(φ_exit)
    y0 = r_exit*sin(φ_exit)
    z0 = 0.
    @threads for ikd in 1:kdNum
        kd0 = sp.ss_kdSamples[ikd]
        kx0 = kd0*-sin(φ_exit)
        ky0 = kd0* cos(φ_exit)
        for ikz in 1:kzNum
            kz0 = sp.ss_kzSamples[ikz]
            init[1:8,ikd,ikz] = [x0,y0,z0,kx0,ky0,kz0,t,ionRate(Ft,kd0,kz0)]
        end
    end
    return reshape(init,dim,:)
end

"Gets the structure factor g of the molecule target of a given channel under a specific Euler angle (β,γ)."
function _getMolStructFactor_g(sp::WFATSampler, nξ::Int, m::Int, β, γ)
    @assert nξ≥0 "[WFATSampler] The nξ must be positive."
    sum = zero(ComplexF64)
    for l in abs(m):sp.wfat_lMax, m_ in -l:l
        sum += sp.wfat_intdata[nξ+1,m+sp.wfat_mMax+1,l+1,m_+l+1] * WignerD.wignerdjmn(l,m,m_,β) * exp(-1im*m_*γ)
    end
    return real(sum)
end
