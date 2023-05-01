
module SampleProviders

using ..Lasers
using ..Targets
using Base.Threads

export initSampleProvider, ElectronSampleProvider, batchNum, generateElectronBatch

function initSampleProvider(;kwargs...)
    return if kwargs[:ionRateMethod] == :ADK
        ADKSampler(;kwargs...)
    elseif kwargs[:ionRateMethod] == :SFA
        SFASampler(;kwargs...)
    elseif kwargs[:ionRateMethod] == :SFA_AE
        SFAAESampler(;kwargs...)
    elseif kwargs[:ionRateMethod] == :WFAT
        WFATSampler(;kwargs...)
    else
        error("[SampleProviders] Undefined tunneling rate method [$(kwargs[:ionRateMethod])].")
        return
    end
end

abstract type ElectronSampleProvider end

include("ADKSampler.jl")
include("SFASampler.jl")
include("SFAAESampler.jl")
include("WFATSampler.jl")

end