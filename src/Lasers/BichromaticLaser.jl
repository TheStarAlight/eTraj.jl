
"""
    struct BichromaticLaser <: Laser

Represents a bichromatic laser field which consists of two `MonochromaticLaser`s.
"""
struct BichromaticLaser <: Laser
    laser1::MonochromaticLaser
    laser2::MonochromaticLaser
    delay
end

"""
    BichromaticLaser(l1::MonochromaticLaser, l2::MonochromaticLaser [,delay=0.0]) <: Laser

Initializes a new `BichromaticLaser` with two `MonochromaticLaser`s.

## Parameters
- `l1, l2::MonochromaticLaser`  : Two `MonochromaticLaser`s.
- `delay`                       : Time delay of `l2` respective to `l1` (numerically in **a.u.** or a `Unitful.Quantity`).

## Examples
```julia-repl
julia> using eTraj.Units

julia> l = BichromaticLaser(l1=Cos4Laser(peak_int=1.0PW/cm^2, wave_len=800nm, cyc_num=10, ellip=1), l2=Cos4Laser(peak_int=1.0PW/cm^2, wave_len=400nm, cyc_num=20, ellip=-1), delay=0.5fs)
[BichromaticLaser] delay Δt = 20.67 a.u. (0.50 fs)
├ [MonochromaticLaser] Envelope cos⁴, peak intensity 1.0e+15 W/cm², wavelen=800 nm, 10 cycle(s), ε=1 [circularly polarized]
└ [MonochromaticLaser] Envelope cos⁴, peak intensity 1.0e+15 W/cm², wavelen=400 nm, 20 cycle(s), ε=-1 [circularly polarized]
```
"""
function BichromaticLaser(;l1::MonochromaticLaser, l2::MonochromaticLaser, delay=0.0)
    (delay isa Quantity) && (delay = auconvert(delay).val)
    BichromaticLaser(l1,l2,delay)
end

function Laser1(l::BichromaticLaser)
    return l.laser1
end

function Laser2(l::BichromaticLaser)
    return l.laser2
end

function Base.getindex(l::BichromaticLaser, i::Integer)
    return if i == 1
        l.laser1
    elseif i == 2
        l.laser2
    else
        error("Index of the laser should be either 1 or 2.")
    end
end

function Delay21(l::BichromaticLaser)
    return l.delay
end

function LaserAx(l::BichromaticLaser)
    Δt = l.delay
    Ax1 = LaserAx(l.laser1)
    Ax2 = LaserAx(l.laser2)
    return t -> Ax1(t) + Ax2(t-Δt)
end

function LaserAy(l::BichromaticLaser)
    Δt = l.delay
    Ay1 = LaserAy(l.laser1)
    Ay2 = LaserAy(l.laser2)
    return t -> Ay1(t) + Ay2(t-Δt)
end

function LaserFx(l::BichromaticLaser)
    Δt = l.delay
    Fx1 = LaserFx(l.laser1)
    Fx2 = LaserFx(l.laser2)
    return t -> Fx1(t) + Fx2(t-Δt)
end

function LaserFy(l::BichromaticLaser)
    Δt = l.delay
    Fy1 = LaserFy(l.laser1)
    Fy2 = LaserFy(l.laser2)
    return t -> Fy1(t) + Fy2(t-Δt)
end

function Base.show(io::IO, l::BichromaticLaser)
    print(io, "[BichromaticLaser]")
    if l.delay != 0
        if isinteger(l.delay)
            @printf(io, " delay Δt = %i a.u. (%.2f fs)", l.delay, l.delay*24.19e-3)
        else
            @printf(io, " delay Δt = %.2f a.u. (%.2f fs)", l.delay, l.delay*24.19e-3)
        end
    end
    println(io)
    print(io, " ├ "); show(io, l[1]); println(io)
    print(io, " └ "); show(io, l[2])
end

function Serialize(l::BichromaticLaser)
    dict = OrderedDict{Symbol,Any}()
    type = typeof(l)
    laser1 = Serialize(l.laser1)
    laser2 = Serialize(l.laser2)
    delay21 = l.delay
    @pack! dict = (type, laser1, laser2, delay21)
    return dict
end
