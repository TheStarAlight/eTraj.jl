using SemiclassicalSFI
using SemiclassicalSFI.Targets
using SemiclassicalSFI.Lasers
using Test

@info "# Testing TrajectorySimulation ..."

@testset verbose=true "TrajectorySimulation" begin

    tmpdir = mktempdir()

    @info "Testing ADK-CTMC ..."
    @testset "ADK-CTMC" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-CTMC_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-CTMC (3D) ..."
    @testset "ADK-CTMC (3D)" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 3,
                sample_t_intv       = (-80,80),
                sample_t_num        = 2000,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5,1.0),
                final_p_num         = (250,250,100),
                ss_kd_max           = 2.0,
                ss_kd_num           = 200,
                ss_kz_max           = 1.0,
                ss_kz_num           = 100,
                output_path         = "$tmpdir/test_ADK-CTMC_3D_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-CTMC (Monte-Carlo) ..."
    @testset "ADK-CTMC (Monte-Carlo)" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                sample_monte_carlo  = true,
                mc_kt_num           = 8000,
                mc_kd_max           = 2.0,
                output_path         = "$tmpdir/test_ADK-CTMC_MonteCarlo_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-CTMC (SAEAtom) ..."
    @testset "ADK-CTMC (SAEAtom)" begin
        t = get_atom("He")
        l = Cos4Laser(peak_int=8e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-CTMC_He_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-CTMC (Bichromatic) ..."
    @testset "ADK-CTMC (Bichromatic)" begin
        t = get_atom("H")
        l1 = Cos4Laser(peak_int=1e15, wave_len=800., cyc_num=10.0, ellip=1.)
        l2 = Cos4Laser(peak_int=1e15, wave_len=400., cyc_num=20.0, ellip=-1.)
        l = BichromaticLaser(l1=l1,l2=l2)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-500,500),
                sample_t_num        = 2000,
                traj_t_final        = 550,
                final_p_max         = (4.0,4.0),
                final_p_num         = (400,400),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-CTMC_Bichromatic_400+800_CP_CR.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = [:Pre, :Jac]
            )
            true
        end
    end

    @info "Testing ADK-CTMC (HDF5 output) ..."
    @testset "ADK-CTMC (HDF5 output)" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_fmt          = :h5,
                output_path         = "$tmpdir/test_ADK-CTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-QTMC ..."
    @testset "ADK-QTMC" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-QTMC_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-SCTS ..."
    @testset "ADK-SCTS" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-SCTS_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SPANE-CTMC ..."
    @testset "SFA-SPANE-CTMC" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPANE,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPANE-CTMC_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SPANE-QTMC ..."
    @testset "SFA-SPANE-QTMC" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPANE,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPANE-QTMC_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SPANE-SCTS ..."
    @testset "SFA-SPANE-SCTS" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPANE,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPANE-SCTS_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SPA-CTMC ..."
    @testset "SFA-SPA-CTMC" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPA,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPA-CTMC_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SPA-QTMC ..."
    @testset "SFA-SPA-QTMC" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPA,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPA-QTMC_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SPA-SCTS ..."
    @testset "SFA-SPA-SCTS" begin
        t = get_atom("H")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPA,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPA-SCTS_4e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOADK-CTMC ..."
    @testset "MOADK-CTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-CTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOADK-QTMC ..."
    @testset "MOADK-QTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-QTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOADK-SCTS ..."
    @testset "MOADK-SCTS" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_ADK-SCTS_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SPANE-CTMC ..."
    @testset "MOSFA-SPANE-CTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPANE,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPANE-CTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SPANE-QTMC ..."
    @testset "MOSFA-SPANE-QTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPANE,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPANE-QTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SPANE-SCTS ..."
    @testset "MOSFA-SPANE-SCTS" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPANE,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPANE-SCTS_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SPA-CTMC ..."
    @testset "MOSFA-SPA-CTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SPA,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPA-CTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SPA-QTMC ..."
    @testset "MOSFA-SPA-QTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPA-QTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SPA-SCTS ..."
    @testset "MOSFA-SPA-SCTS" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_SPA-SCTS_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end


    @info "Testing WFAT-CTMC ..."
    @testset "WFAT-CTMC" begin
        t = LoadMolecule("Molecule_Hydrogen.jld2")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            perform_traj_simulation(
                init_cond_method    = :WFAT,
                laser               = l,
                target              = t,
                dimension           = 2,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.5,2.5),
                final_p_num         = (250,250),
                ss_kd_max           = 2.0,
                ss_kd_num           = 800,
                output_path         = "$tmpdir/test_WFAT-CTMC_Hydrogen_3e14_800nm_2cyc_CP.jld2",
                traj_rtol           = 1e-6
            )
            true
        end
    end

end