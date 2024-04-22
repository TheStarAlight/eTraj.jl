using SemiclassicalSFI
using SemiclassicalSFI.Targets
using SemiclassicalSFI.Lasers
using Test
using Base.Threads

@info "# Testing main method performSFI ..."

@testset verbose=true "performSFI" begin

    if Threads.nthreads() == 1
        @warn "The process is running with single available thread, to enable multi-threading, start julia with parameter `-t N` to use N threads."
    end

    rm("./performSFI_test_output/", recursive=true, force=true)
    mkdir("./performSFI_test_output/")

    @info "Testing ADK-CTMC ..."
    @testset "ADK-CTMC" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_ADK-CTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-CTMC (Large) ..."
    @testset "ADK-CTMC (Large)" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 1000,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.6),
                final_p_num         = (400,400,60),
                ss_kd_max           = 2.0,
                ss_kd_num           = 2000,
                ss_kz_max           = 0.1,
                ss_kz_num           = 100,
                save_path           = "./performSFI_test_output/test_ADK-CTMC_Large_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-CTMC (Monte-Carlo) ..."
    @testset "ADK-CTMC (Monte-Carlo)" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                mc_kt_num           = 8000,
                mc_kd_max           = 2.0,
                mc_kz_max           = 0.1,
                save_path           = "./performSFI_test_output/test_ADK-CTMC_MonteCarlo_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full,
                sample_monte_carlo  = true
            )
            true
        end
    end

    @info "Testing ADK-QTMC ..."
    @testset "ADK-QTMC" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_ADK-QTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing ADK-SCTS ..."
    @testset "ADK-SCTS" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_ADK-SCTS_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFAAE-CTMC ..."
    @testset "SFAAE-CTMC" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_SFAAE-CTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFAAE-QTMC ..."
    @testset "SFAAE-QTMC" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_SFAAE-QTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFAAE-SCTS ..."
    @testset "SFAAE-SCTS" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_SFAAE-SCTS_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-CTMC ..."
    @testset "SFA-CTMC" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_SFA-CTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-QTMC ..."
    @testset "SFA-QTMC" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_SFA-QTMC_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-SCTS ..."
    @testset "SFA-SCTS" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_SFA-SCTS_4e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOADK-CTMC ..."
    @testset "MOADK-CTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOADK-CTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOADK-QTMC ..."
    @testset "MOADK-QTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOADK-QTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOADK-SCTS ..."
    @testset "MOADK-SCTS" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOADK-SCTS_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFAAE-CTMC ..."
    @testset "MOSFAAE-CTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOSFAAE-CTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFAAE-QTMC ..."
    @testset "MOSFAAE-QTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOSFAAE-QTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFAAE-SCTS ..."
    @testset "MOSFAAE-SCTS" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOSFAAE-SCTS_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-CTMC ..."
    @testset "MOSFA-CTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOSFA-CTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :CTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-QTMC ..."
    @testset "MOSFA-QTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOSFA-QTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :QTMC,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MOSFA-SCTS ..."
    @testset "MOSFA-SCTS" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_MOSFA-SCTS_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_phase_method   = :SCTS,
                traj_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end


    @info "Testing WFAT-CTMC ..."
    @testset "WFAT-CTMC" begin
        t = GenericMolecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=3e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :WFAT,
                laser               = l,
                target              = t,
                sample_t_intv       = (-80,80),
                sample_t_num        = 400,
                traj_t_final        = 120,
                final_p_max         = (2.0,2.0,0.2),
                final_p_num         = (200,200,20),
                ss_kd_max           = 2.0,
                ss_kd_num           = 400,
                ss_kz_max           = 0.1,
                ss_kz_num           = 20,
                save_path           = "./performSFI_test_output/test_WFAT-CTMC_Hydrogen_3e14_800nm_2cyc_CP.h5",
                traj_rtol           = 1e-6
            )
            true
        end
    end

end