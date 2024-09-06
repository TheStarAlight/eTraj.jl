
module ElectronSamplers

using ..Lasers
using ..Targets

include("imports.jl")

export init_sampler, ElectronSampler, batch_num, gen_electron_batch

function init_sampler(;kwargs...)
    return if kwargs[:init_cond_method] == :SFA || kwargs[:init_cond_method] == :MOSFA
        SFASampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :SFAAE || kwargs[:init_cond_method] == :MOSFAAE
        SFAAESampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :ADK || kwargs[:init_cond_method] == :MOADK
        ADKSampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :WFAT
        WFATSampler(;kwargs...)
    else
        error("[ElectronSamplers] Undefined initial condition method [$(kwargs[:init_cond_method])].")
        return
    end
end

abstract type ElectronSampler end

include("shared_methods.jl")

include("ADKSampler.jl")
include("SFAAESampler.jl")
include("SFASampler.jl")
include("WFATSampler.jl")

end