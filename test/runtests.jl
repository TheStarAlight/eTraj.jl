using SemiclassicalSFI
using Test
using Base.Threads

if Threads.nthreads() == 1
    @warn "The process is running with single available thread, to enable multi-threading, start julia with parameter `-t N` to use N threads."
end
@info "Running with $(Threads.nthreads()) threads ..."

@testset verbose=true "SemiclassicalSFI" begin
    include("Targets_test.jl")
    include("PySCFMolecularCalculator_test.jl")
    include("Lasers_test.jl")
    include("ElectronSamplers_test.jl")
    include("TrajectorySimulation_test.jl")
end
