# gen_Hydrogen_data.jl
# This script generates the Molecule_Hydrogen.jld2

# Tested on AMD Ryzen 9 7950X with 24 threads.
# WFAT calculation takes ~2 mins, asymp_coeff takes ~4 secs.

using SemiclassicalSFI
using SemiclassicalSFI.Targets

molH2 = GenericMolecule(atoms=["H","H"], atom_coords=[0 0 -0.375; 0 0 0.375], charge=0, name="Hydrogen")
MolInitCalculator!(molH2, basis="cc-pVQZ")
@time MolCalcWFATData!(molH2, orbit_ridx=0, grid_rNum=400, grid_θNum=90, grid_ϕNum=90, sf_nξMax=6, sf_mMax=6, sf_lMax=6)
@time MolCalcAsympCoeff!(molH2, orbit_ridx=0, grid_rNum=90, grid_rReg=(3,8), grid_θNum=90, grid_ϕNum=90, l_max=6)
MolSaveDataAs!(molH2, "Molecule_Hydrogen.jld2")