
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
            MolCalcEnergyData!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, basis="sto-3g")
            @test MolHOMOEnergy(mol) ≈ -0.57443656
            @test MolEnergyLevels(mol)[MolHOMOIndex(mol)+1] ≈ 0.6609100513
        end
        #* WFAT structure factor calculation
        @testset verbose=true "WFAT" begin
            MolCalcWFATData!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=1, sf_mMax=1, sf_lMax=6)
            @test MolWFATStructureFactor_G(mol,0,0,0,0.0,0.0) ≈ 2.274804156 && MolWFATStructureFactor_G(mol,0,0,0,π/2,0.0) ≈ 2.053653020 && MolWFATStructureFactor_G(mol,0,0,1,π/4,0.0) ≈ -0.211109323
        end
        #* Asymptotic coefficients calculation
        @testset verbose=true "Asymptotic Coeff" begin
            MolCalcAsympCoeff!(mol, MCType=MolecularCalculators.PySCFMolecularCalculator, grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
            @test MolAsympCoeff(mol)[1,1] ≈ 2.0786574842 && MolAsympCoeff(mol)[3,3] ≈ 0.25253676
            @test abs2(MolMOADKStructureFactor_B(mol,0,0,0.0,0.0)) ≈ 3.45000461 && abs2(MolMOADKStructureFactor_B(mol,0,0,π/2,0.0)) ≈ 1.648543830
        end
    end

end
