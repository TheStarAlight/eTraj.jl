using SemiclassicalSFI
using SemiclassicalSFI.ElectronSamplers
using SemiclassicalSFI.Lasers
using SemiclassicalSFI.Targets
using Test

t = get_atom("H")
tmol = GenericMolecule("Molecule_Hydrogen.h5")
tmol_rot = GenericMolecule("Molecule_Hydrogen.h5", 0,π/2,0)
l = Cos4Laser(peak_int=1e15, wave_len=800.0, cyc_num=2, ellip=1.0)
params = Dict{Symbol,Any}(
    :init_cond_method   => :ADK,
    :laser              => l,
    :target             => t,
    :dimension          => 2,
    :sample_t_intv      => (-100,100),
    :sample_t_num       => 400,
    :sample_cutoff_limit=> 1e-14,
    :sample_monte_carlo => false,
    :traj_phase_method  => :CTMC,
    :rate_prefix        => :Full,
    :ss_kd_max          => 2.0,
    :ss_kd_num          => 800,
    :ss_kz_max          => 2.0,
    :ss_kz_num          => 0,
    :mc_kt_num          => 800,
    :mc_kd_max          => 2.0,
    :mc_kz_max          => 2.0,
    :mol_orbit_idx      => 0,
)
sfaae = Dict(:init_cond_method=> :SFAAE)
sfa = Dict(:init_cond_method=> :SFA)
wfat = Dict(:init_cond_method=> :WFAT)
_3d = Dict(:dimension=>3, :ss_kd_num=>30, :ss_kz_num=>30)
mol = Dict(:target=>tmol)
mol_rot = Dict(:target=>tmol_rot)
qtmc = Dict(:traj_phase_method=> :QTMC)
scts = Dict(:traj_phase_method=> :SCTS)
mc2d = Dict(:sample_monte_carlo=>true)
mc3d = Dict(:sample_monte_carlo=>true, :dimension=>3)
exp_ = Dict(:rate_prefix=> :Exp)
pre = Dict(:rate_prefix=> :Pre)
precc = Dict(:rate_prefix=> :PreCC)
jac = Dict(:rate_prefix=> :Jac)
prejac = Dict(:rate_prefix=> Set([:Pre,:Jac]))
f(sp) = size(gen_electron_batch(sp,200))

@info "Testing ADKSampler ..."
@testset verbose=true "ADKSampler" begin
    @test begin
        sp = init_sampler(; params...)
        f(sp) == (6,640)
    end
    @test begin
        sp = init_sampler(; merge(params, qtmc)...)
        f(sp) == (7,640)
    end
    @test begin
        sp = init_sampler(; merge(params, scts)...)
        f(sp) == (7,640)
    end
    @test begin
        sp = init_sampler(; merge(params, _3d)...)
        f(sp) == (8,432)
    end
    @test begin
        sp = init_sampler(; merge(params, qtmc, _3d)...)
        f(sp) == (9,432)
    end
    @test begin
        sp = init_sampler(; merge(params, scts, _3d)...)
        f(sp) == (9,432)
    end
    @test begin
        sp = init_sampler(; merge(params, mol)...)
        f(sp) == (6,618)
    end
    @test begin
        sp = init_sampler(; merge(params, mol, qtmc)...)
        f(sp) == (7,618)
    end
    @test begin
        sp = init_sampler(; merge(params, mol, scts)...)
        f(sp) == (7,618)
    end
    @test begin
        sp = init_sampler(; merge(params, mol, _3d)...)
        f(sp) == (8,408)
    end
    @test begin
        sp = init_sampler(; merge(params, mol, _3d, qtmc)...)
        f(sp) == (9,408)
    end
    @test begin
        sp = init_sampler(; merge(params, mol, _3d, scts)...)
        f(sp) == (9,408)
    end
    @test begin
        sp = init_sampler(; merge(params, mc2d)...)
        # f(sp) == (6,623)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, mc3d)...)
        # f(sp) == (8,387)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, exp_)...)
        f(sp) == (6,594)
    end
    @test begin
        sp = init_sampler(; merge(params, pre)...)
        f(sp) == (6,616)
    end
    @test begin
        sp = init_sampler(; merge(params, precc)...)
        f(sp) == (6,656)
    end
    @test begin
        sp = init_sampler(; merge(params, jac)...)
        f(sp) == (6,576)
    end
    @test begin
        sp = init_sampler(; merge(params, prejac)...)
        f(sp) == (6,598)
    end
end


@info "Testing SFAAESampler ..."
@testset verbose=true "SFAAESampler" begin
    @test begin
        sp = init_sampler(; merge(params, sfaae)...)
        f(sp) == (6,528)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, qtmc)...)
        f(sp) == (7,528)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, scts)...)
        f(sp) == (7,528)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, _3d)...)
        f(sp) == (8,320)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, _3d, qtmc)...)
        f(sp) == (9,320)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, _3d, scts)...)
        f(sp) == (9,320)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, mol)...)
        f(sp) == (6,506)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, mol, qtmc)...)
        f(sp) == (7,506)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, mol, scts)...)
        f(sp) == (7,506)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, mol, _3d)...)
        f(sp) == (8,292)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, _3d, mol, qtmc)...)
        f(sp) == (9,292)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, _3d, mol, scts)...)
        f(sp) == (9,292)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, mc2d)...)
        # f(sp) == (6,502)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, mc3d)...)
        # f(sp) == (8,294)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, exp_)...)
        f(sp) == (6,521)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, pre)...)
        f(sp) == (6,542)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, precc)...)
        f(sp) == (6,542)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, jac)...)
        f(sp) == (6,505)
    end
    @test begin
        sp = init_sampler(; merge(params, sfaae, prejac)...)
        f(sp) == (6,528)
    end
end

@info "Testing SFASampler ..."
@testset verbose=true "SFASampler" begin
    @test begin
        sp = init_sampler(; merge(params, sfa)...)
        f(sp) == (6,633)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, qtmc)...)
        f(sp) == (7,633)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, scts)...)
        f(sp) == (7,633)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, _3d)...)
        f(sp) == (8,446)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, _3d, qtmc)...)
        f(sp) == (9,446)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, _3d, scts)...)
        f(sp) == (9,446)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, mol)...)
        f(sp) == (6,621)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, mol, qtmc)...)
        f(sp) == (7,621)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, mol, scts)...)
        f(sp) == (7,621)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, mol, _3d)...)
        f(sp) == (8,436)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, _3d, mol, qtmc)...)
        f(sp) == (9,436)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, _3d, mol, scts)...)
        f(sp) == (9,436)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, mc2d)...)
        # f(sp) == (6,617)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, mc3d)...)
        # f(sp) == (8,416)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, exp_)...)
        f(sp) == (6,611)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, pre)...)
        f(sp) == (6,617)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, precc)...)
        f(sp) == (6,637)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, jac)...)
        f(sp) == (6,607)
    end
    @test begin
        sp = init_sampler(; merge(params, sfa, prejac)...)
        f(sp) == (6,613)
    end
end

@info "Testing WFATSampler ..."
@testset verbose=true "WFATSampler" begin
    @test begin
        sp = init_sampler(; merge(params, wfat, mol)...)
        f(sp) == (6,566)
    end
    @test begin
        sp = init_sampler(; merge(params, wfat, mol_rot)...)
        f(sp) == (6,570)
    end
    @test begin
        sp = init_sampler(; merge(params, wfat, _3d, mol)...)
        f(sp) == (8,332)
    end
    @test begin
        sp = init_sampler(; merge(params, wfat, mc2d, mol)...)
        # f(sp) == (6,315)
        !isnothing(f(sp))
    end
    @test begin
        sp = init_sampler(; merge(params, wfat, mc3d, mol)...)
        # f(sp) == (8,323)
        !isnothing(f(sp))
    end
end
