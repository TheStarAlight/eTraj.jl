
"""
Implementation of semiclassical methods in strong field ionization.
"""
module SemiclassicalSFI

using OrdinaryDiffEq
using DiffEqGPU, CUDA
using LinearAlgebra
using StaticArrays
using Parameters
using Base.Threads
using HDF5
using Dates
using ProgressMeter
using Pkg

include("Lasers.jl")
include("Targets.jl")
include("SampleProviders_Base.jl")
using .Lasers
using .Targets
using .SampleProviders

export performSFI, Lasers, Targets

"""
Performs a semiclassical simulation with given parameters.

# Parameters

## Required params. for all methods:
- `ionRateMethod = <:ADK|:SFA|:SFA_AE>`         : Method of determining ionization rate. Currently only supports ADK.
- `laser::Laser`                                : Parameters of the laser field.
- `target::Target`                              : Parameters of the target.
- `sample_tSpan = (start,stop)`                 : Time span in which electrons are sampled.
- `sample_tSampleNum`                           : Number of time samples.
- `simu_tFinal`                                 : Time when every trajectory simulation ends.
- `finalMomentum_pMax = (pxMax,pyMax,pzMax)`    : Boundaries of final momentum spectrum collecting in three dimensions.
- `finalMomentum_pNum = (pxNum,pyNum,pzNum)`    : Numbers of final momentum spectrum collecting in three dimensions.

## Required params. for step-sampling methods:
- `ss_pdMax`    : Boundary of pd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_pdNum`    : Number of pd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_pzMax`    : Boundary of pz (momentum's component along propagation direction (z ax.)) samples.
- `ss_pzNum`    : Number of pz (momentum's component along propagation direction (z ax.)) samples.

## Required params. for Monte-Carlo-sampling methods:
- `mc_tBatchSize`   : Number of electron samples in a single time sample.
- `mc_ptMax`        : Maximum value of momentum's transversal component (perpendicular to field direction).

## Optional params. for all methods:
- `save_fileName`                           : Output HDF5 file name.
- `save_3D_momentumSpec = false`            : Determines whether 3D momentum spectrum is saved.
- `simu_phaseMethod = <:CTMC|:QTMC|:SCTS>`  : Method of classical trajectories' phase.
- `simu_relTol = 1e-6`                      : Relative error tolerance when solving classical trajectories.
- `simu_nondipole = false`                  : Determines whether non-dipole effect is taken account in the simulation (currently not supported).
- `simu_GPU = false`                        : Determines whether GPU acceleration in trajectory simulation is used, requires `DiffEqGPU` up to v1.18.
- `rate_monteCarlo = false`                 : Determines whether Monte-Carlo sampling is used when generating electron samples.
- `rate_ionRatePrefix = <:ExpRate>`         : Prefix of the exponential term in the ionization rate.
- `rydberg_collect = false`                 : Determines whether rydberg final states are collected.
- `rydberg_prinQNMax`                       : Maximum principle quantum number n to be collected.

## Optional params. for ADK methods:
- `adk_ADKTunExit = <:IpF|:FDM|:Para>`      : Tunneling exit method for ADK methods (when `ionRateMethod==:ADK`).

"""
function performSFI(; # some abbrs.:  req. = required, opt. = optional, params. = parameters.
                        #* req. params. for all methods
                    ionRateMethod       ::Symbol,
                    laser               ::Laser,
                    target              ::Target,
                    sample_tSpan        ::Tuple{<:Real,<:Real},
                    sample_tSampleNum   ::Int,
                    simu_tFinal         ::Real,
                    finalMomentum_pMax  ::Tuple{<:Real,<:Real,<:Real},
                    finalMomentum_pNum  ::Tuple{<:Int,<:Int,<:Int},
                        #* req. params. for step-sampling (ss) methods
                    ss_pdMax            ::Real = 0.,
                    ss_pdNum            ::Int  = 0 ,
                    ss_pzMax            ::Real = 0.,
                    ss_pzNum            ::Int  = 0 ,
                        #* req. params. for Monte-Carlo (mc) methods
                    mc_tBatchSize       ::Int  = 0 ,
                    mc_ptMax            ::Real = 0.,
                        #* opt. params. for all methods
                    save_fileName       ::String = defaultFileName(),
                    save_3D_momentumSpec::Bool   = false,
                    simu_phaseMethod    ::Symbol = :CTMC,
                    simu_relTol         ::Real   = 1e-6,
                    simu_nondipole      ::Bool   = false,
                    simu_GPU            ::Bool   = false,
                    rate_monteCarlo     ::Bool   = false,
                    rate_ionRatePrefix  ::Symbol = :ExpRate,
                    rydberg_collect     ::Bool   = false,
                    rydberg_prinQNMax   ::Int    = 0,
                        #* opt. params. for ADK methods
                    adk_ADKTunExit      ::Symbol = :IpF
                    )
    #* pack up all parameters.
    kwargs = Dict{Symbol,Any}()
    @pack! kwargs= (ionRateMethod, laser, target, sample_tSpan, sample_tSampleNum, simu_tFinal, finalMomentum_pMax, finalMomentum_pNum,
                    ss_pdMax, ss_pdNum, ss_pzMax, ss_pzNum,
                    mc_tBatchSize, mc_ptMax,
                    simu_phaseMethod, simu_relTol, simu_nondipole, simu_GPU, rate_monteCarlo, rate_ionRatePrefix, rydberg_collect, rydberg_prinQNMax,
                    adk_ADKTunExit)
    #* compatibility check
    # GPU acceleration requires [DiffEqGPU] up to v1.18
    if simu_GPU
        dep = Pkg.dependencies()
        for (k,v::Pkg.API.PackageInfo) in dep
            if v.name == "DiffEqGPU"
                if v.version < VersionNumber("1.18")
                    error("Package [DiffEqGPU]'s version is lower than required version [v1.18].")
                end
            end
        end
    end
    #* initialize sample provider.
    sp::ElectronSampleProvider = initSampleProvider(;kwargs...)
    #* launch electrons and summarize.
    #   * prepare storage
    # ionization amplitude (ionAmpFinal is for final data, ionAmpConnect is for temporary cache)
    ionProbFinal, ionProbSumTemp, ionProbCollect =
        if simu_phaseMethod == :CTMC
            zeros(Float64, finalMomentum_pNum),     zeros(Float64, finalMomentum_pNum),     zeros(Float64, tuple(finalMomentum_pNum...,nthreads()))
        else
            zeros(ComplexF64, finalMomentum_pNum),  zeros(ComplexF64, finalMomentum_pNum),  zeros(ComplexF64, tuple(finalMomentum_pNum...,nthreads()))
        end
    # rydberg amplitude
    rydProbFinal, rydProbSumTemp, rydProbCollect =
        if rydberg_collect
            if simu_phaseMethod == :CTMC
                zeros(Float64, rydberg_prinQNMax, rydberg_prinQNMax, 2*rydberg_prinQNMax+1),
                zeros(Float64, rydberg_prinQNMax, rydberg_prinQNMax, 2*rydberg_prinQNMax+1),
                zeros(Float64, rydberg_prinQNMax, rydberg_prinQNMax, 2*rydberg_prinQNMax+1, nthreads())
            else
                zeros(ComplexF64, rydberg_prinQNMax, rydberg_prinQNMax, 2*rydberg_prinQNMax+1),
                zeros(ComplexF64, rydberg_prinQNMax, rydberg_prinQNMax, 2*rydberg_prinQNMax+1),
                zeros(ComplexF64, rydberg_prinQNMax, rydberg_prinQNMax, 2*rydberg_prinQNMax+1, nthreads())
            end
        else
            nothing, nothing, nothing
        end
    # classical rates
    classicalRates = Dict{Symbol,Float64}()
        classicalRates[:ion]                = 0.
        classicalRates[:ion_uncollected]    = 0.
        classicalRates[:ryd]                = 0.
        classicalRates[:ryd_uncollected]    = 0.
    #   * launch electrons and collect
    prog1 = ProgressUnknown(dt=0.2, desc="Launching electrons and collecting...", color = :cyan, spinner = true)
    prog2 = Progress(batchNum(sp); dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1)
    for batchId in 1:batchNum(sp)
        init = generateElectronBatch(sp, batchId)
        launchAndCollect!(init, ionProbFinal, ionProbSumTemp, ionProbCollect,
                                rydProbFinal, rydProbSumTemp, rydProbCollect,
                                classicalRates; kwargs...)
        next!(prog1,spinner=raw"-\|/"); next!(prog2);
    end
    finish!(prog1); finish!(prog2);
    if simu_phaseMethod != :CTMC
        ionProbFinal = abs2.(ionProbFinal)
        if rydberg_collect
            rydProbFinal = abs2.(rydProbFinal)
        end
    end
    #* save as HDF5.
    if isfile(save_fileName)
        @warn "File \"$save_fileName\" already exists. Saving at \"$(defaultFileName())\"."
        save_fileName = "$(defaultFileName()).h5"
    end
    #TODO: add support to save simulation abstract.
    h5write(save_fileName, "px", collect(range(-finalMomentum_pMax[1],finalMomentum_pMax[1], length=finalMomentum_pNum[1])))
    h5write(save_fileName, "py", collect(range(-finalMomentum_pMax[2],finalMomentum_pMax[2], length=finalMomentum_pNum[2])))
    if save_3D_momentumSpec
        h5write(save_fileName, "pz", collect(range(-finalMomentum_pMax[3],finalMomentum_pMax[3], length=finalMomentum_pNum[3])))
        h5write(save_fileName, "ionProb", ionProbFinal)
    end
    h5write(save_fileName, "ionProb_xy", reshape(sum(ionProbFinal, dims=3),size(ionProbFinal)[1:2]))
    if rate_ionRatePrefix == :ExpFull   # will get wrong total ionization rate unless full prefix is included.
        h5write(save_fileName, "classicalIonRate",              classicalRates[:ion])
        h5write(save_fileName, "classicalIonRateUncollected",   classicalRates[:ion_uncollected])
    end
    if rydberg_collect
        h5write(save_fileName, "rydProb", rydProbFinal)
        if rate_ionRatePrefix == :ExpFull
            h5write(save_fileName, "classicalRydRate",            classicalRates[:ryd])
            h5write(save_fileName, "classicalRydRateUncollected", classicalRates[:ryd_uncollected])
        end
    end
end


function launchAndCollect!( init,
                            ionProbFinal,
                            ionProbSumTemp,
                            ionProbCollect,
                            rydProbFinal,
                            rydProbSumTemp,
                            rydProbCollect,
                            classicalRates;
                            laser               ::Laser,
                            target              ::Target,
                            simu_tFinal         ::Real,
                            simu_phaseMethod    ::Symbol,
                            simu_relTol         ::Real,
                            simu_nondipole      ::Bool,
                            simu_GPU            ::Bool,
                            rydberg_collect     ::Bool,
                            rydberg_prinQNMax   ::Int,
                            finalMomentum_pMax  ::Tuple{<:Real,<:Real,<:Real},
                            finalMomentum_pNum  ::Tuple{<:Int,<:Int,<:Int},
                            kwargs...   # kwargs are surplus params.
                            )
    Ip                  = IonPotential(target)
    Fx::Function        = LaserFx(laser)
    Fy::Function        = LaserFy(laser)
    targetF::Function   = TargetForce(target)
    targetP::Function   = TargetPotential(target)
    NuclCharge          = AsympNuclCharge(target)
    # including external function call is infeasible in GPU, thus the external targetF & targetP are replaced by pure Coulomb ones.
    traj::Function =
        if ! simu_nondipole
            if simu_phaseMethod == :CTMC
                function traj_dipole_ctmc(u,p,t)
                    # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                    tFx, tFy, tFz = -NuclCharge*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
                    du1 = u[4]
                    du2 = u[5]
                    du3 = u[6]
                    du4 = tFx-Fx(t)
                    du5 = tFy-Fy(t)
                    du6 = tFz
                    @SVector [du1,du2,du3,du4,du5,du6]
                end
            elseif simu_phaseMethod == :QTMC
                function traj_dipole_qtmc(u,p,t)
                    # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                    tFx, tFy, tFz = -NuclCharge*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
                    du1 = u[4]
                    du2 = u[5]
                    du3 = u[6]
                    du4 = tFx-Fx(t)
                    du5 = tFy-Fy(t)
                    du6 = tFz
                    # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                    du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - NuclCharge*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-0.5))
                    @SVector [du1,du2,du3,du4,du5,du6,du7]
                end
            elseif simu_phaseMethod == :SCTS
                function traj_dipole_scts(u,p,t)
                    # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                    tFx, tFy, tFz = -NuclCharge*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
                    du1 = u[4]
                    du2 = u[5]
                    du3 = u[6]
                    du4 = tFx-Fx(t)
                    du5 = tFy-Fy(t)
                    du6 = tFz
                    # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                    du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - NuclCharge*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-0.5) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                    @SVector [du1,du2,du3,du4,du5,du6,du7]
                end
            end
        else
            #TODO: add support for nondipole simulation.
        end
    batchSize = size(init,2)
    classRates_ion              = zeros(nthreads())
    classRates_ion_uncollected  = zeros(nthreads())
    classRates_ryd              = zeros(nthreads())
    classRates_ryd_uncollected  = zeros(nthreads())

    # create ODE problem and solve the ensemble.
    probDim = (simu_phaseMethod == :CTMC) ? 6 : 7
    trajODEProb::ODEProblem = ODEProblem(traj, (@SVector zeros(Float64,probDim)), (0,simu_tFinal))
    initTraj::Function =
        if probDim == 6
            (prob,i,repeat) -> remake(prob; u0=SVector{6}([init[k,i] for k in 1:6]),     tspan = (init[7,i],simu_tFinal))
        else
            (prob,i,repeat) -> remake(prob; u0=SVector{7}([init[k,i] for k in [1:6;9]]), tspan = (init[7,i],simu_tFinal))
        end
    ensembleProb::EnsembleProblem = EnsembleProblem(trajODEProb, prob_func=initTraj, safetycopy=false)
    sol =
        if ! simu_GPU
            solve(ensembleProb, OrdinaryDiffEq.Tsit5(), EnsembleThreads(), trajectories=batchSize, reltol=simu_relTol, save_everystep=false)
        else
            solve(ensembleProb, DiffEqGPU.GPUTsit5(), DiffEqGPU.EnsembleGPUKernel(0.), trajectories=batchSize, dt=0.1, adaptive=true, save_everystep=false)
        end
    # collect and summarize.
    @threads for i in 1:batchSize
        x0,y0,z0,px0,py0,pz0 = sol.u[i][ 1 ][1:6]
        x, y, z, px, py, pz  = sol.u[i][end][1:6]
        phase = (simu_phaseMethod == :CTMC) ? (0.) : (sol.u[i][end][7])
        if simu_phaseMethod == :SCTS # asymptotic Coulomb phase correction term in SCTS
            sqrtb = (2Ip)^(-0.25)
            g = sqrt(1+2Ip*((y*pz-z*py)^2+(z*px-x*pz)^2+(x*py-y*px)^2))
            phase -= px0*x0+py0*y0+pz0*z0 + AsympNuclCharge(target)*sqrtb*(log(g)+asinh((x*py+y*py+z*pz)/(g*sqrtb)))
        end
        E_inf = (px^2+py^2+pz^2)/2 + targetP(x,y,z)
        r_vec = [x, y, z ]
        p_vec = [px,py,pz]
        L_vec = r_vec × p_vec
        L2    = sum(abs2.(L_vec))
        if E_inf ≥ 0    # finally ionized.
            classRates_ion[threadid()] += init[8,i]
            p_inf = sqrt(2E_inf)
            a_vec = p_vec × L_vec - AsympNuclCharge(target) * r_vec ./ norm(r_vec)
            p_inf_vec = (p_inf/(1+p_inf^2*L2)) .* (p_inf .* (L_vec×a_vec) - a_vec)
            pxIdx = round(Int, (p_inf_vec[1]+finalMomentum_pMax[1])/(finalMomentum_pMax[1]/finalMomentum_pNum[1]*2))
            pyIdx = round(Int, (p_inf_vec[2]+finalMomentum_pMax[2])/(finalMomentum_pMax[2]/finalMomentum_pNum[2]*2))
            pzIdx = round(Int, (p_inf_vec[3]+finalMomentum_pMax[3])/(finalMomentum_pMax[3]/finalMomentum_pNum[3]*2))
            if checkbounds(Bool, ionProbCollect, pxIdx,pyIdx,pzIdx, threadid())
                if simu_phaseMethod == :CTMC
                    ionProbCollect[pxIdx,pyIdx,pzIdx, threadid()] += init[8,i] # ionRate
                else
                    ionProbCollect[pxIdx,pyIdx,pzIdx, threadid()] += sqrt(init[8,i])*exp(1im*init[9,i]) # sqrt(ionRate)*phaseFactor
                end
            else
                classRates_ion_uncollected[threadid()] += init[8,i]
            end
        else            # finally become rydberg.
            classRates_ryd[threadid()] += init[8,i]
            if rydberg_collect
                n = round(Int, AsympNuclCharge(target) / sqrt(-2E_inf))
                l = round(Int, (sqrt(1.0+4L2)-1.0)/2)
                m = round(Int, L_vec[3])
                nIdx = n
                lIdx = l+1
                mIdx = m+rydberg_prinQNMax
                if checkbounds(Bool, rydProbCollect, nIdx,lIdx,mIdx, threadid())
                    if simu_phaseMethod == :CTMC
                        rydProbCollect[nIdx,lIdx,mIdx,threadid()] += init[8,i]
                    else
                        rydProbCollect[nIdx,lIdx,mIdx,threadid()] += sqrt(init[8,i])*exp(1im*init[9,i])
                    end
                else
                    classRates_ryd_uncollected[threadid()] += init[8,i]
                end
            else
                classRates_ryd_uncollected[threadid()] += init[8,i]
            end
        end
    end
    sum!(ionProbSumTemp, ionProbCollect)
    ionProbFinal .+= ionProbSumTemp
    if rydberg_collect
        sum!(rydProbSumTemp, rydProbCollect)
        rydProbFinal .+= rydProbSumTemp
    end
    classicalRates[:ion]                += sum(classRates_ion)
    classicalRates[:ion_uncollected]    += sum(classRates_ion_uncollected)
    classicalRates[:ryd]                += sum(classRates_ryd)
    classicalRates[:ryd_uncollected]    += sum(classRates_ryd_uncollected)
end

function defaultFileName()
    Y,M,D = yearmonthday(now())
    h,m,s = hour(now()), minute(now()), second(now())
    return "SCSFI-$(string(Y,pad=4))$(string(M,pad=2))$(string(D,pad=2))-$(string(h,pad=2))$(string(m,pad=2))$(string(s,pad=2)).h5"
end

end
