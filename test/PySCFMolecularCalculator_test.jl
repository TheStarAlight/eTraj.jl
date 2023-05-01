
using SemiclassicalSFI
using Test

@info "# Testing PySCFMolecularCalculator ..."

@testset verbose=true "PySCFMolecularCalculator" begin
    mc = nothing
    mol = Targets.Molecule(["H","H"], [0 0 -0.375; 0 0 0.375], 0, "Hydrogen")
    #* initialization and energy calculation
    @test begin
        Targets.MolCalcEnergyData!(mol, Targets.MolecularCalculators.PySCFMolecularCalculator, basis="sto-3g")
        Targets.MolHOMOEnergy(mol) ≈ -0.57443656
    end
    #* WFAT structure factor calculation
    @test begin
        Targets.MolCalcWFATData!(mol, grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=1, sf_mMax=1, sf_lMax=6)
        Targets.MolWFATStructureFactor_G(mol,0,0,0,0.0,0.0) ≈ 2.274804156 && Targets.MolWFATStructureFactor_G(mol,0,0,0,π/2,0.0) ≈ 2.053653020 && Targets.MolWFATStructureFactor_G(mol,0,0,1,π/4,0.0) ≈ -0.211109323
    end
end
