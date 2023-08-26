
module SampleProviders

using ..Lasers
using ..Targets
using Base.Threads

export init_sampler, ElectronSampleProvider, batch_num, gen_electron_batch

function init_sampler(;kwargs...)
    return if kwargs[:init_cond_method] == :SFA
        SFASampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :SFAAE
        SFAAESampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :ADK
        ADKSampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :MOSFAAE
        MOSFAAESampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :MOADK
        MOADKSampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :WFAT
        WFATSampler(;kwargs...)
    else
        error("[SampleProviders] Undefined initial condition method [$(kwargs[:init_cond_method])].")
        return
    end
end

abstract type ElectronSampleProvider end

include("shared_methods.jl")

include("SFASampler.jl")
include("SFAAESampler.jl")
include("ADKSampler.jl")
include("MOSFAAESampler.jl")
include("MOADKSampler.jl")
include("WFATSampler.jl")

end