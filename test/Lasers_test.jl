using SemiclassicalSFI
using SemiclassicalSFI.Lasers
using Test

@info "# Testing Lasers ..."

@testset verbose=true "Lasers" begin

    @info "Testing Cos2Laser ..."
    @testset verbose=true "Cos2Laser" begin
        l1 = Cos2Laser(peakInt=4e14, waveLen=800., cycNum=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = Cos2Laser(4e14,800.,2.,1.,π/2,π,10.)
        @test l1 == l2
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNum(l1)        == 2.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing Cos4Laser ..."
    @testset verbose=true "Cos4Laser" begin
        l1 = Cos4Laser(peakInt=4e14, waveLen=800., cycNum=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = Cos4Laser(4e14,800.,2.,1.,π/2,π,10.)
        @test l1 == l2
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNum(l1)        == 2.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing GaussianLaser ..."
    @testset verbose=true "GaussianLaser" begin
        l1 = GaussianLaser(peakInt=4e14, waveLen=800., spreadCycNum=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = GaussianLaser(4e14, 800., 2., 1., π/2, π, 10.)
        @test l1 == l2
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test SpreadCycNum(l1)  == 2.
        @test SpreadDuration(l1)== SpreadCycNum(l1) * Period(l1)
        @test FWHM_Duration(l1) == SpreadDuration(l1) * (2*sqrt(2*log(2)))
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing TrapezoidalLaser ..."
    @testset verbose=true "TrapezoidalLaser" begin
        l1 = TrapezoidalLaser(peakInt=4e14, waveLen=800., cycNumTurnOn=2., cycNumTurnOff=2., cycNumConst=6., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = TrapezoidalLaser(4e14, 800., 2., 2., 6., 1., π/2, π, 10.)
        @test l1 == l2
        @test begin
            show(l1)
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNumTurnOn(l1)  == 2.
        @test CycNumTurnOff(l1) == 2.
        @test CycNumConst(l1)   == 6.
        @test CycNumTotal(l1)   == 10.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

end