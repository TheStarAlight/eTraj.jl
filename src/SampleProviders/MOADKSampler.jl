using ..Targets
using HDF5
using Rotations
using WignerD
using Base.Threads

"Sample provider which yields electron samples through MOADK formula, matching `IonRateMethod=:MOADK`."
struct MOADKSampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::Molecule;     # WFAT only supports [Molecule]
    tSamples        ::AbstractVector;
    ss_pdSamples    ::AbstractVector;
    ss_pzSamples    ::AbstractVector;
    ionRatePrefix   ::Symbol;       # currently supports :ExpRate.
    tunExit         ::Symbol;       # :Para for tunneling, :IpF for over-barrier, automatically specified.
    ionOrbitRelHOMO ::Integer;
    ionOrbit_m      ::Integer;

    function MOADKSampler(; laser               ::Laser,
                            target              ::Molecule,
                            sample_tSpan        ::Tuple{<:Real,<:Real},
                            sample_tSampleNum   ::Integer,
                            rate_ionRatePrefix  ::Symbol,
                            ss_pdMax            ::Real,
                            ss_pdNum            ::Integer,
                            ss_pzMax            ::Real,
                            ss_pzNum            ::Integer,
                            mol_ionOrbitRelHOMO ::Integer,
                            moadk_ionOrbit_m    ::Integer,
                            kwargs...   # kwargs are surplus params.
                            )
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[MOADKSampler] Invalid time sample number $sample_tSampleNum."
        @assert (ss_pdNum>0 && ss_pzNum>0) "[MOADKSampler] Invalid pd/pz sample number $ss_pdNum/$ss_pzNum."
        # if coefficients are not available, calculate it.
        if ! (mol_ionOrbitRelHOMO in MolMOADKAvailableIndices(target))
            MolCalcMOADKCoeff!(target, mol_ionOrbitRelHOMO)
        end
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate, :ExpPre])
            error("[MOADKSampler] Unsupported tunneling rate prefix [$rate_ionRatePrefix].")
        end
        # check Keldysh parameter & over-barrier condition.
        Ip = IonPotential(target, mol_ionOrbitRelHOMO)
        F0 = LaserF0(laser)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        if γ0 ≥ 0.5
            @warn "[MOADKSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[MOADKSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        F_crit = Ip^2/4/(1-sqrt(Ip/2))
        tunExit = :Para
        if F0 ≥ F_crit
            @warn "[MOADKSampler] Peak electric field strength F0=$F0, reaching the over-barrier critical value, weak-field condition unsatisfied. Tunneling exit method switched from [Para] to [IpF]."
            tunExit = :IpF
        elseif F0 ≥ F_crit*2/3
            @warn "[MOADKSampler] Peak electric field strength F0=$F0, reaching 2/3 of over-barrier critical value, weak-field condition not sufficiently satisfied."
        end
        # check moadk m.
        @assert moadk_ionOrbit_m≥0 "[MOADKSampler] `moadk_ionOrbit_m` should be non-negative."
        # finish initialization
        return new( laser, target,
                    range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                    range(-abs(ss_pdMax),abs(ss_pdMax);length=ss_pdNum), range(-abs(ss_pzMax),abs(ss_pzMax);length=ss_pzNum),
                    rate_ionRatePrefix, tunExit,
                    mol_ionOrbitRelHOMO, moadk_ionOrbit_m)
    end
end

"Gets the total number of batches."
function batchNum(sp::MOADKSampler)
    return length(sp.tSamples)
end

"Generates a batch of electrons of `batchId` from `sp` using MOADK method."
function generateElectronBatch(sp::MOADKSampler, batchId::Int)
    t = sp.tSamples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft = hypot(Fxt,Fyt)
    if Ft == 0
        return nothing
    end
    φ_field = atan( Fyt, Fxt)   # direction of field vector F.
    φ_exit  = atan(-Fyt,-Fxt)   # direction of tunneling exit, which is opposite to F.
    Ip = IonPotential(sp.target, sp.ionOrbitRelHOMO)
    κ  = sqrt(2Ip)
    Z  = AsympNuclCharge(sp.target)

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
    # the total tunneling rate consists of partial rates of different m' : Γ = ∑ Γ_m'
    # the partial rate consists of a structural part |B_m'|²/(2^|m'|*|m'|!) and a field part W_m'(F) = κ^(-|m'|) * (2κ²/F)^(2Z/κ-|m'|-1) * exp(-2κ³/3F)
    B(m_) = MolMOADKStructureFactor_B(sp.target, sp.ionOrbitRelHOMO, sp.ionOrbit_m, m_, β, γ)
    lMax = size(MolMOADKCoeffs(sp.target, sp.ionOrbitRelHOMO), 1) - 1
    ionRate::Function =
        if sp.ionRatePrefix == :ExpRate || sp.ionRatePrefix == :ExpPre
            function (F,pd,pz)
                Γsum = 0.
                for m_ in -lMax:lMax
                    Γsum += abs2(B(m_))/(2^abs(m_)*factorial(abs(m_))) * κ^(-abs(m_)) * (2κ^2/F)^(2Z/κ-abs(m_)-1) * exp(-2(κ^2+pd^2+pz^2)^1.5/3F)
                end
                return Γsum
            end
        else
            #TODO: Add Full prefix including jac.
        end
    # generating samples
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
