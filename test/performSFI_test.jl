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

    rm("test_ADK.h5", force=true)
    rm("test_SFAAE.h5", force=true)
    rm("test_SFA.h5", force=true)
    rm("test_MOADK.h5", force=true)
    rm("test_WFAT.h5", force=true)

    @info "Testing ADK ..."
    @testset "ADK" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :ADK,
                laser               = l,
                target              = t,
                sample_t_span       = (-80,80),
                sample_t_num        = 400,
                simu_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "test_ADK.h5",
                simu_phase_method   = :CTMC,
                simu_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA-AE ..."
    @testset "SFA-AE" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFAAE,
                laser               = l,
                target              = t,
                sample_t_span       = (-80,80),
                sample_t_num        = 400,
                simu_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "test_SFAAE.h5",
                simu_phase_method   = :CTMC,
                simu_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing SFA ..."
    @testset "SFA" begin
        t = HAtom()
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :SFA,
                laser               = l,
                target              = t,
                sample_t_span       = (-80,80),
                sample_t_num        = 400,
                simu_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "test_SFA.h5",
                simu_phase_method   = :CTMC,
                simu_rtol           = 1e-6,
                rate_prefix         = :Full
            )
            true
        end
    end

    @info "Testing MO-ADK ..."
    @testset "MO-ADK" begin
        t = Molecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :MOADK,
                laser               = l,
                target              = t,
                sample_t_span       = (-80,80),
                sample_t_num        = 400,
                simu_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "test_MOADK.h5",
                simu_rtol           = 1e-6
            )
            true
        end
    end

    @info "Testing WFAT ..."
    @testset "WFAT" begin
        t = Molecule("Molecule_Hydrogen.h5")
        l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2, ellip=1.0)
        @test begin
            performSFI(
                init_cond_method    = :WFAT,
                laser               = l,
                target              = t,
                sample_t_span       = (-80,80),
                sample_t_num        = 400,
                simu_t_final        = 120,
                final_p_max         = (2.0,2.0,2.0),
                final_p_num         = (200,200,1),
                ss_kd_max           = 2.0,
                ss_kd_num           = 100,
                ss_kz_max           = 2.0,
                ss_kz_num           = 100,
                save_path           = "test_WFAT.h5",
                simu_rtol           = 1e-6
            )
            true
        end
    end

end