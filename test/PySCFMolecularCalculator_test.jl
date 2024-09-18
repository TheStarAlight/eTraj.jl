
using eTraj
using eTraj.Targets
using Unitful
using Test

@info "# Testing PySCFMolecularCalculator ..."  # PySCF version 2.3.0

if Sys.iswindows()
    @warn "The OS is Windows, the PySCFMolecularCalculator is inoperable, skipping."
    @testset verbose=true "PySCFMolecularCalculator" begin
        @test_skip 0==0
    end
else
    @testset verbose=true "PySCFMolecularCalculator" begin
        @testset verbose=true "RHF" begin
            mol = GenericMolecule(atoms=["H","H"], atom_coords=[0 0 -0.375; 0 0 0.375]*u"Å", charge=0, name="Hydrogen")
            #* initialization and energy calculation
            MolInitCalculator!(mol, basis="cc-pVTZ")
            @test MolHOMOEnergy(mol) ≈ -0.591627497396511
            @test MolEnergyLevel(mol,1) ≈ 0.1658433878161975
            #* WFAT structure factor calculation
            MolCalcWFATData!(mol,0; grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=3, sf_mMax=3, sf_lMax=6)
            @test MolWFATStructureFactor_G(mol,0,0,0,0.0,0.0) ≈ 1.850581827979092
            @test MolWFATStructureFactor_G(mol,0,0,0,π/2,0.0) ≈ 1.6049709989987115
            @test MolWFATStructureFactor_G(mol,0,0,1,π/4,0.0) ≈ -0.23551608527272053
            #* Asymptotic coefficients calculation
            MolCalcAsympCoeff!(mol,0; grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
            @test MolAsympCoeff(mol,0)[1,1] ≈ 1.0902102938182585
            @test MolAsympCoeff(mol,0)[3,3] ≈ 0.10307388691997624
        end
        @testset verbose=true "UHF" begin
            mol = GenericMolecule(atoms=["O","O"], atom_coords=[0 0 -0.605; 0 0 0.605]*u"Å", charge=0, spin=1, name="Oxygen")
            #* initialization and energy calculation
            MolInitCalculator!(mol, basis="cc-pVDZ")
            @test MolHOMOEnergy(mol,1) ≈ -0.5501877457376334
            @test MolEnergyLevel(mol,(1,-1)) ≈ MolHOMOEnergy(mol,1)
            #* WFAT structure factor calculation
            MolCalcWFATData!(mol,(1, 0); grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=3, sf_mMax=3, sf_lMax=6)
            MolCalcWFATData!(mol,(1,-1); grid_rNum=100, grid_θNum=30, grid_ϕNum=30, sf_nξMax=3, sf_mMax=3, sf_lMax=6)
            @test MolWFATStructureFactor_G(mol,(1, 0),0,1,0.0,0.0) |> abs2 ≈ 8.681792866824848
            @test MolWFATStructureFactor_G(mol,(1, 0),0,1,π/4,0.0) |> abs2 ≈ 3.6641554825114815
            @test MolWFATStructureFactor_G(mol,(1, 0),0,1,π/2,0.0) |> abs2 ≈ 0.0                atol=1e-10
            @test MolWFATStructureFactor_G(mol,(1,-1),0,1,0.0,0.0) |> abs2 ≈ 8.681792866824809
            #* Asymptotic coefficients calculation
            MolCalcAsympCoeff!(mol,(1, 0); grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
            MolCalcAsympCoeff!(mol,(1,-1); grid_rNum=50, grid_rReg=(3,8), grid_θNum=30, grid_ϕNum=30, l_max=6)
            @test MolAsympCoeff(mol,(1, 0))[3,2] ≈ -0.7004560834188851im
            @test MolAsympCoeff(mol,(1, 0))[3,4] ≈ -0.7004560834189529im
            @test MolAsympCoeff(mol,(1,-1))[3,2] ≈ -0.7004560834191522
            @test MolAsympCoeff(mol,(1,-1))[3,4] ≈  0.7004560834191522
        end
    end
end
