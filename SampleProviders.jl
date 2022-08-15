
module SampleProviders

# include("Lasers.jl")
# include("Targets.jl")
using ..Lasers
using ..Targets
using Base.Threads

export initSampleProvider, ElectronSampleProvider, batchNum, generateElectronBatch


function initSampleProvider(;kwargs...)
    return if kwargs[:ionRateMethod] == :ADK
        ADKSampleProvider(;kwargs...)
    else
        error("Undefined tunneling rate method [$(kwargs[:ionRateMethod])].")
        return
    end
end


abstract type ElectronSampleProvider end

begin :ADK
    "Sample provider which yields electron samples through ADK rate formula, matching `IonRateMethod=:ADK`."
    struct ADKSampleProvider <: ElectronSampleProvider
        laser           ::Laser;
        target          ::SAEAtom;          # ADK only supports [SAEAtom].
        tSamples        ::AbstractRange;
        pdSamples       ::AbstractRange;
        pzSamples       ::AbstractRange;
        phaseMethod     ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
        ionRatePrefix   ::Symbol;           # currently supports :ExpRate.
        ADKTunExit      ::Symbol;           # currently supports :IpF, :FDM, :ADK.
        function ADKSampleProvider(;laser ::Laser,
                                    target::SAEAtom,
                                    sample_tSpan::Tuple{<:Real,<:Real},
                                    sample_tNum ::Int,
                                    sample_pdMax::Real,
                                    sample_pdNum::Int,
                                    sample_pzMax::Real,
                                    sample_pzNum::Int,
                                    simu_phaseMethod  ::Symbol,
                                    rate_ionRatePrefix::Symbol,
                                    adk_ADKTunExit::Symbol,
                                    kwargs...   # kwargs are surplus params.
                                    )
            F0 = LaserF0(laser)
            Ip = target.IonPotential
            γ0 = LaserAngFreq(laser) * sqrt(2Ip) / F0
            # check phase method support.
            if ! (simu_phaseMethod in [:CTMC])
                error("Undefined phase method [$simu_phaseMethod].")
                return
            end
            # check tunneling exit support.
            if ! (adk_ADKTunExit in [:IpF, :FDM, :Para])
                error("Undefined tunneling exit method [$adk_ADKTunExit].")
                return
            end
            # switch tunneling exit method if inappropriate.
            if adk_ADKTunExit == :FDM && Ip^2 < 4F0
                adk_ADKTunExit == :IpF
                @warn "Tunneling exit method has changed from [FDM] to [IpF] due to failure in tunneling exit determination of [FDM] model."
            elseif adk_ADKTunExit == :Para && Ip^2 < 4*(1-sqrt(Ip/2))*F0
                adk_ADKTunExit == :IpF
                @warn "Tunneling exit method has changed from [Para] to [IpF] due to failure in tunneling exit determination of [Para] model."
            end
            # check IonRate prefix support.
            if ! (rate_ionRatePrefix in [:ExpRate])
                error("Undefined tunneling rate prefix [$rate_ionRatePrefix].")
                return
            end
            # check Keldysh parameter.
            if γ0 ≥ 0.5
                @warn "Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
            elseif γ0 ≥ 1.0
                @warn "Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
            end
            # finish initialization.
            return new(laser,target,
                       range(sample_tSpan[1],sample_tSpan[2];length=sample_tNum),
                       range(-abs(sample_pdMax),abs(sample_pdMax);length=sample_pdNum),
                       range(-abs(sample_pzMax),abs(sample_pzMax);length=sample_pzNum),
                       simu_phaseMethod,rate_ionRatePrefix,adk_ADKTunExit)
        end
    end
    "Gets the total number of batches."
    function batchNum(sp::ADKSampleProvider)
        return length(sp.tSamples)
    end
    "Generates a batch of electrons of `batchId` from `sp` using ADK method."
    function generateElectronBatch(sp::ADKSampleProvider, batchId::Int)
        t  = sp.tSamples[batchId]
        Fxt::Function = LaserFx(sp.laser)
        Fyt::Function = LaserFy(sp.laser)
        Fx = Fxt(t)
        Fy = Fyt(t)
        F  = hypot(Fx,Fy)
        φ  = atan(-Fy,-Fx)
        Ip = sp.target.IonPotential
        # determining tunneling exit position.
        r_exit =
            if      sp.ADKTunExit == :IpF
                Ip / F
            elseif  sp.ADKTunExit == :FDM
                (Ip + sqrt(Ip^2 - 4F)) / 2F
            elseif  sp.ADKTunExit == :Para
                (Ip + sqrt(Ip^2 - 4*(1-sqrt(Ip/2))*F)) / 2F
            end
        # determining ADK rate. ADKRate(F,φ,pd,pz)
        ADKRate::Function =
            if      sp.ionRatePrefix == :ExpRate
                ADKRateExp(sp.target)
            else
                #TODO: Add support for other prefixes.
            end
        dim = (sp.phaseMethod == :CTMC) ? 8 : 9 # x,y,z,px,py,pz,t0,rate[,phase]
        pdNum, pzNum = length(sp.pdSamples), length(sp.pzSamples)
        init = zeros(Float64, dim, pdNum, pzNum) # initial condition
        x0 = r_exit*cos(φ)
        y0 = r_exit*sin(φ)
        z0 = 0.
        @threads for ipd in 1:pdNum
            pd0 = sp.pdSamples[ipd]
            px0 = pd0*-sin(φ)
            py0 = pd0* cos(φ)
            for ipz in 1:pzNum
                pz0 = sp.pzSamples[ipz]
                init[1:8,ipd,ipz] = [x0,y0,z0,px0,py0,pz0,t,ADKRate(F,φ,pd0,pz0)]
                if sp.phaseMethod != :CTMC
                    init[9,ipd,ipz] = 0
                end
            end
        end
        return reshape(init,dim,:)
    end
end

end