
module ElectronSamplers

using ..Lasers
using ..Targets

include("imports.jl")

export init_sampler, ElectronSampler, batch_num, gen_electron_batch

function init_sampler(;kwargs...)
    return if kwargs[:init_cond_method] == :SPA || kwargs[:init_cond_method] == :MOSPA || kwargs[:init_cond_method] == :SFA || kwargs[:init_cond_method] == :MOSFA
        SPASampler(;kwargs...)
    elseif kwargs[:init_cond_method] == :SPANE || kwargs[:init_cond_method] == :MOSPANE || kwargs[:init_cond_method] == :SFAAE || kwargs[:init_cond_method] == :MOSFAAE
        SPANESampler(;kwargs...)
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
include("SPANESampler.jl")
include("SPASampler.jl")
include("WFATSampler.jl")

end