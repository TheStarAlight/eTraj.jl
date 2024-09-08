using SemiclassicalSFI
using SemiclassicalSFI.Lasers
using Unitful
using Test

@info "# Testing Lasers ..."

@testset verbose=true "Lasers" begin

    @info "Testing Cos2Laser ..."
    @testset verbose=true "Cos2Laser" begin
        l1 = Cos2Laser(peak_int=4e14, wave_len=800., cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = Cos2Laser(peak_int=4e14, ang_freq=0.05695419065625, cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l3 = Cos2Laser(peak_int=4e14, wave_len=800., duration=220.63996467273427, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l4 = Cos2Laser(peak_int=4e14u"W/cm^2", wave_len=800.0u"nm", duration=5.337025523633234u"fs", ellip=1., azi=π/2, cep=π, t_shift=241.88843265767443u"as")
        l5 = Cos2Laser(peak_int=4e14u"W/cm^2", ang_freq=1.5498024802806117u"eV", duration=5.337025523633234u"fs", ellip=1., azi=π/2, cep=π, t_shift=241.88843265767443u"as")
        l = Cos2Laser(4e14,800.,2.,1.,π/2,π,10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test l == l4
        @test l == l5
        @test begin
            show(l1)
            println()
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNum(l1)        == 2.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing Cos4Laser ..."
    @testset verbose=true "Cos4Laser" begin
        l1 = Cos4Laser(peak_int=4e14, wave_len=800., cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = Cos4Laser(peak_int=4e14, ang_freq=0.05695419065625, cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l3 = Cos4Laser(peak_int=4e14, wave_len=800., duration=220.63996467273427, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l4 = Cos4Laser(peak_int=4e14u"W/cm^2", wave_len=800.0u"nm", duration=5.337025523633234u"fs", ellip=1., azi=π/2, cep=π, t_shift=241.88843265767443u"as")
        l5 = Cos4Laser(peak_int=4e14u"W/cm^2", ang_freq=1.5498024802806117u"eV", duration=5.337025523633234u"fs", ellip=1., azi=π/2, cep=π, t_shift=241.88843265767443u"as")
        l = Cos4Laser(4e14,800.,2.,1.,π/2,π,10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test l == l4
        @test l == l5
        @test begin
            show(l1)
            println()
            true
        end
        @test PeakInt(l1)       == 4e14
        @test WaveLen(l1)       == 800.
        @test CycNum(l1)        == 2.
        @test Ellipticity(l1)   == 1.
        @test Azimuth(l1)       == π/2
        @test AngFreq(l1)       == 45.563352525 / WaveLen(l1)
        @test Period(l1)        == 2π / AngFreq(l1)
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing GaussianLaser ..."
    @testset verbose=true "GaussianLaser" begin
        l1 = GaussianLaser(peak_int=4e14, wave_len=800., spread_cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l2 = GaussianLaser(peak_int=4e14, ang_freq=0.05695419065625, spread_cyc_num=2., ellip=1., azi=π/2, cep=π, t_shift=10.)
        l3 = GaussianLaser(peak_int=4e14, wave_len=800., spread_duration=220.63996467273427, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l4 = GaussianLaser(peak_int=4e14, wave_len=800., FWHM_duration=519.5674115462751, ellip=1., azi=π/2, cep=π, t_shift=10.)
        l5 = GaussianLaser(peak_int=4e14u"W/cm^2", wave_len=800.0u"nm", spread_duration=5.337025523633234u"fs", ellip=1., azi=π/2, cep=π, t_shift=10.)
        l6 = GaussianLaser(peak_int=4e14u"W/cm^2", wave_len=800.0u"nm", FWHM_duration=12.567734683893338u"fs", ellip=1., azi=π/2, cep=π, t_shift=10.)
        l = GaussianLaser(4e14, 800., 2., 1., π/2, π, 10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test l == l4
        @test l == l5
        @test l == l6
        @test begin
            show(l1)
            println()
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
        @test CEP(l1)           == π
        @test TimeShift(l1)     == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing TrapezoidalLaser ..."
    @testset verbose=true "TrapezoidalLaser" begin
        l1 = TrapezoidalLaser(peak_int=4e14, wave_len=800., cyc_num_turn_on=2., cyc_num_turn_off=2., cyc_num_const=6., ellip=1., azi=π/2, cep=π, t_turn_on=10.)
        l2 = TrapezoidalLaser(peak_int=4e14, ang_freq=0.05695419065625, cyc_num_turn_on=2., cyc_num_turn_off=2., cyc_num_const=6., ellip=1., azi=π/2, cep=π, t_turn_on=10.)
        l3 = TrapezoidalLaser(peak_int=4e14u"W/cm^2", wave_len=800.0u"nm", cyc_num_turn_on=2., cyc_num_turn_off=2., cyc_num_const=6., ellip=1., azi=π/2, cep=π, t_turn_on=10.)
        l4 = TrapezoidalLaser(peak_int=4e14u"W/cm^2", ang_freq=1.5498024802806117u"eV", cyc_num_turn_on=2., cyc_num_turn_off=2., cyc_num_const=6., ellip=1., azi=π/2, cep=π, t_turn_on=0.24188843265767443u"fs")
        l = TrapezoidalLaser(4e14, 800., 2., 2., 6., 1., π/2, π, 10.)
        @test l == l1
        @test l == l2
        @test l == l3
        @test l == l4
        @test begin
            show(l1)
            println()
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
        @test CEP(l1)           == π
        @test TimeTurnOn(l1)    == 10.0
        @test LaserF0(l1)       == sqrt(PeakInt(l1)/(1.0+Ellipticity(l1)^2)/3.50944521e16)
        @test LaserA0(l1)       == LaserF0(l1) / AngFreq(l1)
    end

    @info "Testing BichromaticLaser ..."
    @testset verbose=true "BichromaticLaser" begin
        l1 = Cos4Laser(peak_int=4e14, wave_len=800., cyc_num=6., ellip=1.)
        l2 = Cos4Laser(peak_int=4e14, wave_len=400., cyc_num=6., ellip=-1.)
        l_ = BichromaticLaser(l1=l1, l2=l2, delay=50.0)
        l__ = BichromaticLaser(l1=l1, l2=l2, delay=1.2094421632883721u"fs")
        l = BichromaticLaser(l1,l2,50.0)
        @test l === l_
        @test l === l__
        @test Laser1(l) === l[1] === l1
        @test Laser2(l) === l[2] === l2
        @test Delay21(l) == 50.0
        @test begin
            show(l)
            println()
            true
        end
    end

end