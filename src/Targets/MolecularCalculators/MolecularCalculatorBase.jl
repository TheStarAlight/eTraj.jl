
@doc """
    abstract type MolecularCalculatorBase

Supertype of all molecular calculators.
"""
abstract type MolecularCalculatorBase end

include("PySCFMolecularCalculator.jl")
