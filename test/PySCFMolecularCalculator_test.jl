
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
        mol = Molecule(["H","H"], [0 0 -0.375; 0 0 0.375], 0, "Hydrogen")
        #* initialization and energy calculation
        @testset verbose=true "Energy" begin
            MolCalcEnergyData!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, basis="pc-1")
            @test MolHOMOEnergy(mol) ≈ -0.590975727629967
            @test MolEnergyLevels(mol)[MolHOMOIndex(mol)+1] ≈ 0.17418844361014238
        end
        #* WFAT structure factor calculation
        @testset verbose=true "WFAT" begin
            MolCalcWFATData!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=1, sf_mMax=1, sf_lMax=6)
            @test MolWFATStructureFactor_G(mol,0,0,0,0.0,0.0) ≈ 1.9040322797077103 && MolWFATStructureFactor_G(mol,0,0,0,π/2,0.0) ≈ 2.053653020 && MolWFATStructureFactor_G(mol,0,0,1,π/4,0.0) ≈ -0.24183919604716864
        end
        #* Asymptotic coefficients calculation
        @testset verbose=true "Asymptotic Coeff" begin
            MolCalcAsympCoeff!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
            @test MolAsympCoeff(mol)[1,1] ≈ 1.1146439877360017 && MolAsympCoeff(mol)[3,3] ≈ 0.12640164662355585
        end
    end

end
