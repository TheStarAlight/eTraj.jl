# gen_mol_data.jl
# This script is not included in the program.
# This script generates the molecule data.

using eTraj
using eTraj.Targets
using Unitful

mols = Dict(
"Hydrogen"          => GenericMolecule(atoms=["H","H"], atom_coords=[0 0 -0.37072; 0 0 0.37072]*u"Å", charge=0, name="Hydrogen (H₂)"),
"Nitrogen"          => GenericMolecule(atoms=["N","N"], atom_coords=[0 0 -0.54885; 0 0 0.54885]*u"Å", charge=0, name="Nitrogen (N₂)"),
"Oxygen"            => GenericMolecule(atoms=["O","O"], atom_coords=[0 0 -0.6037; 0 0 0.6037]*u"Å", charge=0, spin=1, name="Oxygen (O₂)"),
"Carbon Monoxide"   => GenericMolecule(atoms=["C","O"], atom_coords=[0 0 -0.56415; 0 0 0.56415]*u"Å", charge=0, name="Carbon Monoxide (CO)"),
"Nitric Oxide"      => GenericMolecule(atoms=["N","O"], atom_coords=[0 0 -0.5753; 0 0 0.5753]*u"Å", charge=0, spin=1/2, name="Nitric Oxide (NO)"),
"Hydrochloric Acid" => GenericMolecule(atoms=["H","Cl"], atom_coords=[0 0 -0.6373; 0 0 0.6373]*u"Å", charge=0, name="Hydrochloric Acid (HCl)"),
"Carbon Dioxide"    => GenericMolecule(atoms=["O","C","O"], atom_coords=[0 0 -1.1600; 0 0 0; 0 0 1.1600]*u"Å", charge=0, name="Carbon Dioxide (CO₂)"),
"Sulfur Dioxide"    => GenericMolecule(atoms=["O","S","O"], atom_coords=[-1.2349 0 0.72264; 0 0 0; 1.2349 0 0.72264]*u"Å", charge=0, name="Sulfur Dioxide (SO₂)"),     # on xz plane, bisector of ∠OSO towards +z axis. r(SO) = 1.4308 Å, ∠OSO = 119.329°
"Water"             => GenericMolecule(atoms=["H","O","H"], atom_coords=[-0.7571 0 0.5861; 0 0 0; 0.7571 0 0.5861]*u"Å", charge=0, name="Water (H₂O)"),   # on xz plane, bisector of ∠HOH towards +z axis. r(HO) = 0.9575 Å, ∠HOH = 104.51°
"Ammonia"           => GenericMolecule(atoms=["N","H","H","H"], atom_coords=[0 0 0; 0.9375 0 0.3810; -0.4688 0.8119 0.3810; -0.4688 -0.8119 0.3810]*u"Å", charge=0, name="Ammonia (NH₃)"),    # N is in the origin, with the symmetric axis towards +z, and an H atom on xz plane. r(NH) = 1.012 Å, ∠HNH = 106.7°.
"Acetylene"         => GenericMolecule(atoms=["H","C","C","H"], atom_coords=[0 0 -1.6615; 0 0 -0.6015; 0 0 0.6015; 0 0 1.6615]*u"Å", charge=0, name="Acetylene (C₂H₂)"), # r(CC) = 1.203 Å, r(CH) = 1.060 Å
"Methane"           => GenericMolecule(atoms=["C","H","H","H","H"], atom_coords=[0 0 0; 0 0 1.0870; 1.0250 0 -0.36182; -0.51251 -0.88769 -0.36182; -0.51251 0.88769 -0.36182]*u"Å", charge=0, name="Methane (CH₄)"), # C is in the origin, with one C-H bond towards +z, and another H atom on xz plane. r(CH) = 1.0870 Å, ∠HNH = 109.5°.
"Benzene"           => GenericMolecule(atoms=["C","C","C","C","C","C","H","H","H","H","H","H"], atom_coords=[-0.700 0 1.2124; 0.700 0 1.2124; -0.700 0 -1.2124; 0.700 0 -1.2124; 1.399 0 0; -1.399 0 0; -1.24 0 1.2124; 1.24 0 1.2124; -1.24 0 -1.2124; 1.24 0 -1.2124; -2.48 0 0; 2.48 0 0]*u"Å", charge=0, name="Benzene (C₆H₆)") # on xz plane, with two opposite C-H bonds on x axis. r(CC) = 1.399 Å, r(CH) = 1.101 Å
)

param_prec3 = Dict(:grid_rNum=>400, :grid_θNum=>90, :grid_ϕNum=>90, :sf_nξMax=>6, :sf_mMax=>6, :sf_lMax=>6)
param_prec2 = Dict(:grid_rNum=>300, :grid_θNum=>75, :grid_ϕNum=>75, :sf_nξMax=>6, :sf_mMax=>6, :sf_lMax=>6)
param_prec1 = Dict(:grid_rNum=>200, :grid_θNum=>60, :grid_ϕNum=>60, :sf_nξMax=>6, :sf_mMax=>6, :sf_lMax=>6)

# Tested on AMD Ryzen 9 7950X with 24 threads on WSL Ubuntu 22.04.1 LTS with PySCF 2.3.0

molH2 = mols["Hydrogen"]
MolInitCalculator!(molH2, basis="cc-pVQZ")
@time MolCalcAsympCoeff!(molH2, 0; param_prec3...)
@time MolCalcWFATData!(molH2, 0; param_prec3...)
MolSaveDataAs!(molH2, "Molecule_Hydrogen.jld2")

molN2 = mols["Nitrogen"]
MolInitCalculator!(molN2, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molN2, 0; param_prec3...)
@time MolCalcAsympCoeff!(molN2,-1; param_prec3...)
@time MolCalcAsympCoeff!(molN2,-2; param_prec3...)
@time MolCalcWFATData!(molN2, 0; param_prec3...)
@time MolCalcWFATData!(molN2,-1; param_prec3...)
@time MolCalcWFATData!(molN2,-2; param_prec3...)
MolSaveDataAs!(molN2, "Molecule_Nitrogen.jld2")

molO2 = mols["Oxygen"]
MolInitCalculator!(molO2, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molO2,(1, 0); param_prec3...)
@time MolCalcAsympCoeff!(molO2,(1,-1); param_prec3...)
@time MolCalcWFATData!(molO2,(1, 0); param_prec3...)
@time MolCalcWFATData!(molO2,(1,-1); param_prec3...)
MolSaveDataAs!(molO2, "Molecule_Oxygen.jld2")

molCO = mols["Carbon Monoxide"]
MolInitCalculator!(molCO, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molCO, 0; param_prec3...)
@time MolCalcAsympCoeff!(molCO,-1; param_prec3...)
@time MolCalcAsympCoeff!(molCO,-2; param_prec3...)
@time MolCalcWFATData!(molCO, 0; param_prec3...)
@time MolCalcWFATData!(molCO,-1; param_prec3...)
@time MolCalcWFATData!(molCO,-2; param_prec3...)
MolSaveDataAs!(molCO, "Molecule_CarbonMonoxide.jld2")

molNO = mols["Nitric Oxide"]
MolInitCalculator!(molNO, basis="cc-pVTZ")  # NO's α-HOMO & α-LUMO are degenerate.
@time MolCalcAsympCoeff!(molNO, (1,0); param_prec3...)
@time MolCalcAsympCoeff!(molNO, (1,1); param_prec3...)
@time MolCalcWFATData!(molNO, (1,0); param_prec3...)
@time MolCalcWFATData!(molNO, (1,1); swap_HOMO_LUMO=true, param_prec3...)   # turn on swap_HOMO_LUMO because WFAT involves interactions between orbitals.
MolSaveDataAs!(molNO, "Molecule_NitricOxide.jld2")

molHCl = mols["Hydrochloric Acid"]
MolInitCalculator!(molHCl, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molHCl, 0; param_prec3...)
@time MolCalcAsympCoeff!(molHCl,-1; param_prec3...)
@time MolCalcWFATData!(molHCl, 0; param_prec3...)
@time MolCalcWFATData!(molHCl,-1; param_prec3...)
MolSaveDataAs!(molHCl, "Molecule_HydrochloricAcid.jld2")

molCO2 = mols["Carbon Dioxide"]
MolInitCalculator!(molCO2, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molCO2, 0; param_prec3...)
@time MolCalcAsympCoeff!(molCO2,-1; param_prec3...)
@time MolCalcWFATData!(molCO2, 0; param_prec2...)
@time MolCalcWFATData!(molCO2,-1; param_prec2...)
MolSaveDataAs!(molCO2, "Molecule_CarbonDioxide.jld2")

molSO2 = mols["Sulfur Dioxide"]
MolInitCalculator!(molSO2, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molSO2, 0; merge(param_prec3, Dict(:grid_rReg=>(5,10)))...)
@time MolCalcAsympCoeff!(molSO2,-1; param_prec3...)
@time MolCalcAsympCoeff!(molSO2,-2; param_prec3...)
@time MolCalcWFATData!(molSO2, 0; param_prec2...)
@time MolCalcWFATData!(molSO2,-1; param_prec2...)
@time MolCalcWFATData!(molSO2,-2; param_prec2...)
MolSaveDataAs!(molSO2, "Molecule_SulfurDioxide.jld2")

molH2O = mols["Water"]
MolInitCalculator!(molH2O, basis="aug-cc-pVDZ")
@time MolCalcAsympCoeff!(molH2O, 0; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molH2O,-1; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molH2O,-2; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcWFATData!(molH2O, 0; param_prec2...)
@time MolCalcWFATData!(molH2O,-1; param_prec2...)
@time MolCalcWFATData!(molH2O,-2; param_prec2...)
MolSaveDataAs!(molH2O, "Molecule_Water.jld2")

molNH3 = mols["Ammonia"]
MolInitCalculator!(molNH3, basis="aug-cc-pVTZ")
@time MolCalcAsympCoeff!(molNH3, 0; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molNH3,-1; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molNH3,-2; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcWFATData!(molNH3, 0; param_prec2...)
@time MolCalcWFATData!(molNH3,-1; param_prec2...)
@time MolCalcWFATData!(molNH3,-2; param_prec2...)
MolSaveDataAs!(molNH3, "Molecule_Ammonia.jld2")

molC2H2 = mols["Acetylene"]
MolInitCalculator!(molC2H2, basis="cc-pVDZ")
@time MolCalcAsympCoeff!(molC2H2, 0; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molC2H2,-1; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcWFATData!(molC2H2, 0; param_prec2...)
@time MolCalcWFATData!(molC2H2,-1; param_prec2...)
MolSaveDataAs!(molC2H2, "Molecule_Acetylene.jld2")

molCH4 = mols["Methane"]
MolInitCalculator!(molCH4, basis="cc-pVDZ")
@time MolCalcAsympCoeff!(molCH4, 0; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molCH4,-1; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molCH4,-2; merge(param_prec3, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcWFATData!(molCH4, 0; param_prec2...)
@time MolCalcWFATData!(molCH4,-1; param_prec2...)
@time MolCalcWFATData!(molCH4,-2; param_prec2...)
MolSaveDataAs!(molCH4, "Molecule_Methane.jld2")

molC6H6 = mols["Benzene"]
MolInitCalculator!(molC6H6, basis="cc-pVTZ")
@time MolCalcAsympCoeff!(molC6H6, 0; merge(param_prec2, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molC6H6,-1; merge(param_prec2, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcAsympCoeff!(molC6H6,-2; merge(param_prec2, Dict(:grid_rReg=>(3,10)))...)
@time MolCalcWFATData!(molC6H6, 0; param_prec1...)
@time MolCalcWFATData!(molC6H6,-1; param_prec1...)
@time MolCalcWFATData!(molC6H6,-2; param_prec1...)
MolSaveDataAs!(molC6H6, "Molecule_Benzene.jld2")

