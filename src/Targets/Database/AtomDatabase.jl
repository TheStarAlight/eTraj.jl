# This file stores commonly used atoms.

"""
    get_atom(name::String, charge::Integer=0; kwargs...)

Constructs a `HydrogenLikeAtom` or `SAEAtom` in the database based on the `name` and `charge` given, other initialization parameters can be passed via `kwargs`.
For `HydrogenLikeAtom`, the init params include `Ip`,`Z` and `name`; for `SAEAtom`, the init params include `Ip`,`Z`,`a1`,`b1`,`a2`,`b2`,`a3`,`b3`,`l` and `name`.
The available keys can be found by invoking [`get_available_atoms`](@ref).

## Examples

```jldoctest

julia> get_atom("H")
[HydrogenLikeAtom] Atom H, Ip=0.5000, Z=1, soft_core=0.2000

julia> get_atom("He")
[SAEAtom] Atom He, Ip=0.9036, Z=1

julia> get_atom("Ar"; m=1, asymp_coeff=0.950) # asymp_coeff from [Phys.-Uspekhi 47, 855 (2004)] (l=0 case)
[SAEAtom] Atom Ar, Ip=0.5792, Z=1

julia> get_atom("Xe", 1)
[SAEAtom] Atom Xe⁺, Ip=0.7708, Z=2

julia> get_atom("Xe", 2)
ERROR: [get_atom] Atom name key `Xe2p` not found in database.
[...]

```
See also [`get_available_atoms`](@ref), [`HydrogenLikeAtom`](@ref) and [`SAEAtom`](@ref).

"""
function get_atom(name::String, charge::Integer=0; kwargs...)
    key = name
    if charge > 0
        key *= string(charge)*"p"
    elseif charge < 0
        key *= string(-charge)*"n"
    end
    if key in keys(HydLikeAtomData)
        return HydrogenLikeAtom(;HydLikeAtomData[key]..., kwargs...)
    end
    if key in keys(SAEAtomData)
        return SAEAtom(;SAEAtomData[key]..., kwargs...)
    end
    error("[get_atom] Atom name key `$key` not found in database.")
end

"Gets the available atom keys in the atom database, which are used to access atom objects using [`get_atom`](@ref)."
function get_available_atoms()
    return vcat(collect(keys(HydLikeAtomData)),collect(keys(SAEAtomData)))
end

const HydLikeAtomData = OrderedDict{String, Any}(
    "H" => OrderedDict{Symbol, Any}(
        :Ip => 0.5,
        :Z  => 1,
        :name => "H"
    ),
    "He1p" => OrderedDict{Symbol, Any}(
        :Ip => 1.0,
        :Z  => 2,
        :name => "He⁺"
    ),
    "Li2p" => OrderedDict{Symbol, Any}(
        :Ip => 1.5,
        :Z  => 3,
        :name => "Li²⁺"
    ),
    )
const SAEAtomData = OrderedDict{String, Any}(
    "He" => OrderedDict{Symbol, Any}(
        :Ip => 0.9035698802,
        :Z  => 1,
        :a1 =>  1.230723,
        :b1 =>  0.6620055,
        :a2 => -1.325040,
        :b2 =>  1.236224,
        :a3 => -0.2307230,
        :b3 =>  0.4804286,
        :l  => 0,   # 1s²
        :name => "He"
    ),
    "Ne" => OrderedDict{Symbol, Any}(
        :Ip => 0.79248225,
        :Z  => 1,
        :a1 =>  8.069321,
        :b1 =>  2.148149,
        :a2 => -3.569610,
        :b2 =>  1.985509,
        :a3 =>  0.9306791,
        :b3 =>  0.6018819,
        :l  => 1,   # [He]2s²2p⁶
        :name => "Ne"
    ),
    "Ne1p" => OrderedDict{Symbol, Any}(
        :Ip => 1.505361,
        :Z  => 2,
        :a1 =>  8.042896,
        :b1 =>  2.714635,
        :a2 =>  0.5064300,
        :b2 =>  0.9818268,
        :a3 => -0.04289633,
        :b3 =>  0.4012269,
        :l  => 1,   # [He]2s²2p⁵
        :name => "Ne⁺"
    ),
    "Ne2p" => OrderedDict{Symbol, Any}(
        :Ip => 2.3307635,
        :Z  => 3,
        :a1 => 7.000000,
        :b1 => 3.007574,
        :a2 => 0.6521344,
        :b2 => 1.048617,
        :a3 => 0.0,
        :b3 => 0.0,
        :l  => 1,   # [He]2s²2p⁴
        :name => "Ne²⁺"
    ),
    "Ar" => OrderedDict{Symbol, Any}(
        :Ip => 0.579155055,
        :Z  => 1,
        :a1 =>  16.03862,
        :b1 =>   2.007423,
        :a2 => -25.54278,
        :b2 =>   4.524663,
        :a3 =>   0.9613848,
        :b3 =>   0.4426447,
        :l  => 1,   # [Ne]3s²3p⁶
        :name => "Ar"
    ),
    "Ar1p" => OrderedDict{Symbol, Any}(
        :Ip => 1.0153715,
        :Z  => 2,
        :a1 =>  14.98869,
        :b1 =>   2.217293,
        :a2 => -23.60558,
        :b2 =>   4.584761,
        :a3 =>   1.011314,
        :b3 =>   0.5510695,
        :l  => 1,   # [Ne]3s²3p⁵
        :name => "Ar⁺"
    ),
    "Ar2p" => OrderedDict{Symbol, Any}(
        :Ip => 1.497,
        :Z  => 3,
        :a1 =>  15.00000,
        :b1 =>  3.588835,
        :a2 =>  5.951071,
        :b2 =>  1.692914,
        :a3 =>  0.0,
        :b3 =>  0.0,
        :l  => 1,   # [Ne]3s²3p⁴
        :name => "Ar²⁺"
    ),
    "V" => OrderedDict{Symbol, Any}(
        :Ip => 0.2479178,
        :Z  => 1,
        :a1 =>  21.27861,
        :b1 =>   1.810071,
        :a2 => -30.02658,
        :b2 =>   3.201500,
        :a3 =>   0.7213853,
        :b3 =>   0.4534235,
        :l  => 2,   # [Ar]4s²3d³
        :name => "V"
    ),
    "Ni" => OrderedDict{Symbol, Any}(
        :Ip => 0.28076035,
        :Z  => 1,
        :a1 =>  26.81231,
        :b1 =>   1.937326,
        :a2 => -39.55201,
        :b2 =>   3.161884,
        :a3 =>   0.1876890,
        :b3 =>   0.2626299,
        :l  => 2,   # [Ar]4s²3d⁸
        :name => "Ni"
    ),
    "Kr" => OrderedDict{Symbol, Any}(
        :Ip => 1.02895202,
        :Z  => 1,
        :a1 =>  23.6188,
        :b1 =>  10.8391,
        :a2 => 109.956,
        :b2 =>   6.55432,
        :a3 =>  11.8312,
        :b3 =>   1.27723,
        :l  => 1,   # [Ar]4s²3d¹⁰4p⁶
        :name => "Kr"
    ),
    "Kr1p" => OrderedDict{Symbol, Any}(
        :Ip => 1.790416,
        :Z  => 2,
        :a1 =>  23.1460,
        :b1 =>  11.0416,
        :a2 => 110.447,
        :b2 =>   6.67063,
        :a3 =>  10.8540,
        :b3 =>   1.41524,
        :l  => 1,   # [Ar]4s²3d¹⁰4p⁵
        :name => "Kr⁺"
    ),
    "Rb" => OrderedDict{Symbol, Any}(
        :Ip => 0.153506625,
        :Z  => 1,
        :a1 =>  24.02321,
        :b1 =>  11.10663,
        :a2 => 115.1998,
        :b2 =>   6.629123,
        :a3 =>  11.97679,
        :b3 =>   1.244882,
        :l  => 0,   # [Kr]5s¹
        :name => "Rb"
    ),
    "Nb" => OrderedDict{Symbol, Any}(
        :Ip => 0.248383,
        :Z  => 1,
        :a1 =>  4.000000,
        :b1 =>  1.379674,
        :a2 => -6.356689,
        :b2 =>  2.245839,
        :a3 =>  0.0,
        :b3 =>  0.0,
        :l  => 2,   # [Kr]5s¹4d⁴
        :name => "Nb"
    ),
    "Pd" => OrderedDict{Symbol, Any}(
        :Ip => 0.3063732,
        :Z  => 1,
        :a1 =>  43.83606,
        :b1 =>   2.395037,
        :a2 => -93.35653,
        :b2 =>   4.907321,
        :a3 =>   1.163943,
        :b3 =>   0.3767380,
        :l  => 2,   # [Kr]4d¹⁰
        :name => "Pd"
    ),
    "Xe" => OrderedDict{Symbol, Any}(
        :Ip => 0.445763535,
        :Z  => 1,
        :a1 =>  51.35554,
        :b1 =>   2.111554,
        :a2 => -99.92747,
        :b2 =>   3.737221,
        :a3 =>   1.644457,
        :b3 =>   0.4306465,
        :l  => 1,   # [Kr]5s²4d¹⁰5p⁶
        :name => "Xe"
    ),
    "Xe1p" => OrderedDict{Symbol, Any}(
        :Ip => 0.7708,
        :Z  => 2,
        :a1 =>  37.49194,
        :b1 =>   9.452474,
        :a2 => 103.1967,
        :b2 =>   5.305537,
        :a3 =>  14.50806,
        :b3 =>   1.315683,
        :l  => 1,   # [Kr]5s²4d¹⁰5p⁵
        :name => "Xe⁺"
    ),
    "Ta" => OrderedDict{Symbol, Any}(
        :Ip => 0.27744165,
        :Z  => 1,
        :a1 =>   67.70192,
        :b1 =>    2.565180,
        :a2 => -135.5191,
        :b2 =>    4.088800,
        :a3 =>    4.298078,
        :b3 =>    0.7303277,
        :l  => 2,   # [Xe]4f¹⁴5d³
        :name => "Ta"
    ),
)