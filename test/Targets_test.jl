using SemiclassicalSFI
using SemiclassicalSFI.Targets
using Unitful
using Test

@info "# Testing Targets ..."

@testset verbose=true "Targets" begin

    @info "Testing HydrogenLikeAtom ..."
    @testset verbose=true "HydrogenLikeAtom" begin
        t1 = HydrogenLikeAtom(Ip=0.5, Z=1, l=0, m=0, asymp_coeff=:hartree, quan_ax_θ=π/2, quan_ax_ϕ=π, soft_core=0.2, name="H")
        t1_= HydrogenLikeAtom(Ip=13.605693123044498u"eV", Z=1, l=0, m=0, asymp_coeff=:hartree, quan_ax_θ=π/2, quan_ax_ϕ=π, soft_core=0.2, name="H")
        t2 = HydrogenLikeAtom(0.5,1,0,0,1.0, π/2,π, 0.2,"H")
        t3 = get_atom("H"; quan_ax_θ=π/2, quan_ax_ϕ=π)
        @test t1 == t1_ == t2 == t3
        @test begin
            show(t1)
            println()
            true
        end
        @test IonPotential(t1)          == 0.5
        @test AsympNuclCharge(t1)       == 1
        @test AngularQuantumNumber(t1)  == 0
        @test MagneticQuantumNumber(t1) == 0
        @test QuantizationAxisOrientaion(t1) == (π/2,π)
        @test AsympCoeff(t1)            == 1.0
        @test SoftCore(t1)              == 0.2
        @test TargetName(t1)            == "H"
        @test TargetPotential(t1)(1.0,1.0,1.0) ≈ -1/sqrt(3.2)
        @test reduce(*, TargetForce(t1)(1.0,1.0,1.0) .≈ (-1.0,-1.0,-1.0) .* (3.2)^(-1.5))
    end

    @info "Testing SAEAtom ..."
    @testset verbose=true "SAEAtom" begin
        t1 = SAEAtom(Ip=0.9035698802, Z=1, l=0, m=0, asymp_coeff=:hartree, quan_ax_θ=π/2, quan_ax_ϕ=π, a1=1.230723, b1=0.6620055, a2=-1.325040, b2=1.236224, a3=-0.2307230, b3=0.4804286, name="He")
        t1_= SAEAtom(Ip=24.58738901045456u"eV", Z=1, l=0, m=0, asymp_coeff=:hartree, quan_ax_θ=π/2, quan_ax_ϕ=π, a1=1.230723, b1=0.6620055, a2=-1.325040, b2=1.236224, a3=-0.2307230, b3=0.4804286, name="He")
        t2 = SAEAtom(0.9035698802, 1, 0, 0, 0.9124458182520485, π/2, π, 1.230723, 0.6620055, -1.325040, 1.236224, -0.2307230, 0.4804286, "He")
        t3 = get_atom("He"; quan_ax_θ=π/2, quan_ax_ϕ=π)
        @test t1 == t2 == t3
        @test begin
            show(t1)
            println()
            true
        end
        @test IonPotential(t1)          == 0.9035698802
        @test AsympNuclCharge(t1)       == 1
        @test AngularQuantumNumber(t1)  == 0
        @test MagneticQuantumNumber(t1) == 0
        @test QuantizationAxisOrientaion(t1) == (π/2,π)
        @test AsympCoeff(t1)            == 0.9124458182520485
        @test TargetName(t1)            == "He"
        @test TargetPotential(t1)(1.0,1.0,1.0) ≈ - (1.0 + 1.230723*exp(-0.6620055*1.73205081) + -1.325040*1.73205081*exp(-1.236224*1.73205081) + -0.2307230*exp(-0.4804286*1.73205081)) / 1.73205081
        @test reduce(*, TargetForce(t1)(1.0,1.0,1.0) .≈ (-1.0,-1.0,-1.0) .* (1.73205081^(-3) * (1.0 + 1.230723*(1+0.6620055*1.73205081)*exp(-0.6620055*1.73205081) + -0.2307230*(1+0.4804286*1.73205081)*exp(-0.4804286*1.73205081)) + -1.325040*1.236224/1.73205081 * exp(-1.236224*1.73205081)))
    end

    @info "Testing Atom Database ..."
    @testset verbose=true "Atom Database" begin
        t1 = get_atom("H")
        @test TargetName(t1)        == "H"
        @test AsympNuclCharge(t1)   == 1
        t2 = get_atom("Li", 2)
        @test TargetName(t2)        == "Li²⁺"
        @test AsympNuclCharge(t2)   == 3
        t3 = get_atom("Li2p")
        @test t3 == t2
        t4 = get_atom("Xe")
        @test TargetName(t4)        == "Xe"
        @test AsympNuclCharge(t4)   == 1
        t5 = get_atom("Xe", 1)
        @test TargetName(t5)        == "Xe⁺"
        @test AsympNuclCharge(t5)   == 2
        t6 = get_atom("Xe1p")
        @test t6 == t5
        t7 = get_atom("Xe"; m=1, asymp_coeff=1.0)
        @test MagneticQuantumNumber(t7) == 1
        @test AsympCoeff(t7) == 1.0
    end

    @info "Testing GenericMolecule ..."
    @testset verbose=true "GenericMolecule" begin
        # fresh init
        m1 = GenericMolecule(atoms=["C","C","C","C","C","C","H","H","H","H","H","H"],
                    atom_coords=[0 0 1.4 ; 1.212 0 0.7 ; 1.212 0 -0.7 ; 0 0 -1.4 ; -1.212 0 -0.7 ; -1.212 0 0.7 ; 0 0 2.48; 2.148 0 1.24; 2.148 0 -1.24; 0 0 -2.48; -2.148 0 -1.24; -2.148 0 1.24],
                    charge=0, name="Benzene")
        @test MolAtoms(m1) == ["C","C","C","C","C","C","H","H","H","H","H","H"]
        @test MolAtomCoords(m1) == [0 0 1.4 ; 1.212 0 0.7 ; 1.212 0 -0.7 ; 0 0 -1.4 ; -1.212 0 -0.7 ; -1.212 0 0.7 ; 0 0 2.48; 2.148 0 1.24; 2.148 0 -1.24; 0 0 -2.48; -2.148 0 -1.24; -2.148 0 1.24]
        @test MolCharge(m1) == 0
        @test TargetName(m1) == "Benzene"
        # init from existing file
        m2 = LoadMolecule("Molecule_Hydrogen.jld2")
        @test begin
            show(m2)
            println()
            true
        end
        @test MolAtoms(m2)      == ["H","H"]
        @test MolAtomCoords(m2) == [0.0 0.0 -0.375; 0.0 0.0 0.375]
        @test MolCharge(m2)     == 0
        @test MolEnergyDataAvailable(m2) == true
        @test MolEnergyLevels(m2)[1] ≈ -0.591819099386013
        @test MolHOMOEnergy(m2) ≈ -0.591819099386013
        @test MolWFATAvailableIndices(m2) == Set([0])
        begin
            μ,IntData = MolWFATData(m2,0)
            @test reduce(*, abs.(μ) .<= [1e-10,1e-10,1e-10])
            @test IntData[1,1,1,1] ≈ 1.7276641139621103
        end
        @test abs2(MolWFATStructureFactor_G(m2,0,0,0,0.0,0.0)) ≈ 3.594605830985567
        @test abs2(MolWFATStructureFactor_G(m2,0,0,0,π/2,0.0)) ≈ 2.707676191146582
        @test MolAsympCoeffAvailableIndices(m2) == Set([0])
        @test MolAsympCoeff(m2, 0)[1,1] ≈ 1.0843994462656439
        begin
            @test begin
                SetMolRotation!(m2,π,π/2,π/3)
                reduce(*, MolRotation(m2) .≈ (π,π/2,π/3))
            end
            @test begin
                SetMolRotation!(m2,(π/2,π/3,π/4))
                reduce(*, MolRotation(m2) .≈ (π/2,π/3,π/4))
            end
        end
        @test MolExportAtomInfo(m2) == "H [0.0, 0.0, -0.375]; H [0.0, 0.0, 0.375]"
        @test IonPotential(m2)      == -MolHOMOEnergy(m2)
        @test IonPotential(m2,0)    == -MolHOMOEnergy(m2)
        @test AsympNuclCharge(m2)   == 1
        @test TargetName(m2)        == "Hydrogen"
        @test TargetPotential(m2)(1.0,1.0,1.0) == -0.5
        @test reduce(*, TargetForce(m2)(1.0,1.0,1.0) .== (-1/8,-1/8,-1/8))
        @test begin
            tmp_file, tmp_io = mktemp()
            close(tmp_io)
            MolSaveDataAs!(m2, tmp_file, true)
            true
        end
    end
end