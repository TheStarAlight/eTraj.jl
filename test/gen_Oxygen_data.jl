# gen_Oxygen_data.jl
# This script generates the Molecule_Oxygen.jld2

# Tested on AMD Ryzen 9 7950X with 24 threads.
# WFAT calculation of each orbital takes ~3 mins, asymp coeff takes ~12 secs.

using SemiclassicalSFI
using SemiclassicalSFI.Targets

molO2 = GenericMolecule(atoms=["O","O"], atom_coords=[0 0 -0.605; 0 0 0.605]*u"Å", charge=0, spin=1, name="Oxygen")
MolInitCalculator!(molO2, basis="cc-pVQZ")
@time MolCalcWFATData!(molO2, orbit_ridx=(1,0), grid_rNum=200, grid_θNum=90, grid_ϕNum=90, sf_nξMax=6, sf_mMax=6, sf_lMax=6)
@time MolCalcWFATData!(molO2, orbit_ridx=(1,-1), grid_rNum=200, grid_θNum=90, grid_ϕNum=90, sf_nξMax=6, sf_mMax=6, sf_lMax=6)
@time MolCalcAsympCoeff!(molO2, orbit_ridx=(1,0), grid_rNum=200, grid_rReg=(3,8), grid_θNum=90, grid_ϕNum=90, l_max=6)
@time MolCalcAsympCoeff!(molO2, orbit_ridx=(1,-1), grid_rNum=200, grid_rReg=(3,8), grid_θNum=90, grid_ϕNum=90, l_max=6)
MolSaveDataAs!(molO2, "Molecule_Oxygen.jld2")