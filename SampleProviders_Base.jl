
module SampleProviders

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

include("SampleProviders_ADK.jl")

end