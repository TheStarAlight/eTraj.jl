# gen_Hydrogen_data.jl
# This script generates the Molecule_Hydrogen.h5

using SemiclassicalSFI
using SemiclassicalSFI.Targets

mol = GenericMolecule(atoms=["H","H"], atom_coords=[0 0 -0.375; 0 0 0.375], charge=0, name="Hydrogen")
MolCalcEnergyData!(mol, MCType=PySCFMolecularCalculator, basis="cc-pVTZ")
MolCalcWFATData!(mol, MCType=PySCFMolecularCalculator, grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=3, sf_mMax=3, sf_lMax=6)
MolCalcAsympCoeff!(mol, MCType=PySCFMolecularCalculator, grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
MolSaveDataAs!(mol, "Molecule_Hydrogen.h5")