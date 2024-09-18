
using Printf: @printf, @sprintf
using StaticArrays: @SVector
using OrderedCollections: OrderedDict
using Parameters: @pack!, @unpack
using SpecialFunctions: gamma
using Unitful: Quantity, uconvert, eV, @u_str
using UnitfulAtomic: auconvert

# GenericMolecule
using JLD2
using WignerD: wignerdjmn
using Dates: yearmonthday, now, hour, minute, second
using Rotations: RotZYZ

# PySCFMolecularCalculator
using SphericalHarmonics: sphericalharmonic
using Folds
using Base.Threads
using Einsum: @einsum, @einsimd
using LsqFit: curve_fit, coef, confidence_interval
using libcint_jll
using PyCall: pyimport
using ProgressMeter: Progress, ProgressUnknown, BarGlyphs, next!, finish!
