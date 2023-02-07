
module SampleProviders

using ..Lasers
using ..Targets
using Base.Threads

export initSampleProvider, ElectronSampleProvider, batchNum, generateElectronBatch

function initSampleProvider(;kwargs...)
    return if kwargs[:ionRateMethod] == :ADK
        ADKSampleProvider(;kwargs...)
    elseif kwargs[:ionRateMethod] == :SFA
        SFASampleProvider(;kwargs...)
    elseif kwargs[:ionRateMethod] == :SFA_AE
        SFAAESampleProvider(;kwargs...)
    elseif kwargs[:ionRateMethod] == :WFAT
        WFATSampleProvider(;kwargs...)
    else
        error("[SampleProviders] Undefined tunneling rate method [$(kwargs[:ionRateMethod])].")
        return
    end
end

abstract type ElectronSampleProvider end

include("ADKSampleProvider.jl")
include("SFASampleProvider.jl")
include("SFAAESampleProvider.jl")
include("WFATSampleProvider.jl")

end