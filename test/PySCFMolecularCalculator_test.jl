
using SemiclassicalSFI
using SemiclassicalSFI.Targets
using Test

@info "# Testing PySCFMolecularCalculator ..."

if Sys.iswindows()
    @warn "The OS is Windows, the PySCFMolecularCalculator is inoperable."
    @testset verbose=true "PySCFMolecularCalculator" begin
        @test_skip 0==0
    end
else

    @testset verbose=true "PySCFMolecularCalculator" begin
        mol = GenericMolecule(atoms=["H","H"], atom_coords=[0 0 -0.375; 0 0 0.375], charge=0, name="Hydrogen")
        #* initialization and energy calculation
        @testset verbose=true "Energy" begin
            MolCalcEnergyData!(mol, MCType=PySCFMolecularCalculator, basis="cc-pVTZ")
            @test MolHOMOEnergy(mol) ≈ -0.591627497396511
            @test MolEnergyLevels(mol)[MolHOMOIndex(mol)+1] ≈ 0.16584338781619645
        end
        #* WFAT structure factor calculation
        @testset verbose=true "WFAT" begin
            MolCalcWFATData!(mol, MCType=PySCFMolecularCalculator, grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=3, sf_mMax=3, sf_lMax=6)
            @test MolWFATStructureFactor_G(mol,0,0,0,0.0,0.0) ≈ 1.850581827979092 && MolWFATStructureFactor_G(mol,0,0,0,π/2,0.0) ≈ 1.6049709989987115 && MolWFATStructureFactor_G(mol,0,0,1,π/4,0.0) ≈ -0.23551608527272053
        end
        #* Asymptotic coefficients calculation
        @testset verbose=true "Asymptotic Coeff" begin
            MolCalcAsympCoeff!(mol, MCType=PySCFMolecularCalculator, grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
            @test MolAsympCoeff(mol)[1,1] ≈ 1.0902102938182585 && MolAsympCoeff(mol)[3,3] ≈ 0.10307388691997624
        end
    end

end
