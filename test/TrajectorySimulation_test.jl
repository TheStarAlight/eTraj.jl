using SemiclassicalSFI
using Test
using StaticArrays
using OrdinaryDiffEq
using DiffEqGPU

@info "# Testing trajectory simulation ..."

@testset verbose=true "Trajectory simulation" begin
    # l = Lasers.Cos2Laser(2e14,800,2,1.0)
    traj = function (u,p,t)
        F0 = 0.05338027007325633 # Lasers.LaserF0(l)
        ω = 0.05695419065625 # Lasers.AngFreq(l)
        N = 2   # l.cycNum
        φ = 0.0 # l.cep
        ε = 1.0 # l.ellip
        tFx, tFy, tFz = -1*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
        du1 = u[4]
        du2 = u[5]
        du3 = u[6]
        du4 = tFx- F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*sin(ω*t+φ) + 1/N*sin(ω*t/(2N))*cos(ω*t+φ))
        du5 = tFy- F0 * cos(ω*t/(2N)) * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*cos(ω*t+φ) - 1/N*sin(ω*t/(2N))*sin(ω*t+φ)) * ε
        du6 = tFz
        @SVector [du1,du2,du3,du4,du5,du6]
    end
    simu_tFinal = 120.
    simu_relTol = 1e-6
    init = [10.0  9.0  8.0  7.0;    # x
             5.0  1.0  0.0  0.0;    # y
             0.0  0.0  6.0  1.0;    # z
             0.0  1.0  0.0  1.0;    # px
             1.0  1.0 -1.0  2.0;    # py
             0.0  0.0  0.0  3.0;    # pz
             0.0 -5.0  0.0  0.0;]
    final = [-107.586580634  14.769123285  -111.052988864   10.541324293;
               88.693028463  66.732062023  -135.010858308  216.489065097;
                0.0           0.0             0.580939406  357.985889013;
               -0.972128611  -0.017974487    -0.992418369    0.028651930;
                0.852993520   0.666820202    -0.932743294    1.985856245;
                0.0           0.0            -0.047583462    2.972879184;]
    trajODEProb::ODEProblem = ODEProblem(traj, (@SVector zeros(Float64,6)), (0,simu_tFinal))
    initTraj = (prob,i,repeat) -> remake(prob; u0=SVector{6}([init[k,i] for k in 1:6]), tspan=(init[7,i],simu_tFinal))
    ensembleProb::EnsembleProblem = EnsembleProblem(trajODEProb, prob_func=initTraj, safetycopy=false)
    solc = nothing
    solg = nothing
    @info "Testing CPU..."
    @testset verbose=true "CPU" begin
        @test begin
            solc = solve(ensembleProb, OrdinaryDiffEq.Tsit5(), EnsembleThreads(), trajectories=size(init,2), reltol=simu_relTol, save_everystep=false);
            mapreduce(function((k,i),) ≈(final[k,i],solc[i][end][k],rtol=1e-4) end, *, [(k,i) for k in 1:6, i in 1:size(init,2)])
        end
    end
    @info "Testing GPU..."
    @testset verbose=true "GPU" begin
        @test begin
            solg = solve(ensembleProb, DiffEqGPU.GPUTsit5(), DiffEqGPU.EnsembleGPUKernel(0.), trajectories=size(init,2), dt=0.1, adaptive=true, save_everystep=false);
            mapreduce(function((k,i),) ≈(final[k,i],solg[i][end][k],rtol=1e-2) end, *, [(k,i) for k in 1:6, i in 1:size(init,2)])
        end
    end
end