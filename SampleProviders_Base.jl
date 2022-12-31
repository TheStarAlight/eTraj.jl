
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
    else
        error("[SampleProviders] Undefined tunneling rate method [$(kwargs[:ionRateMethod])].")
        return
    end
end

abstract type ElectronSampleProvider end

include("SampleProviders_ADK.jl")
include("SampleProviders_SFA.jl")

end