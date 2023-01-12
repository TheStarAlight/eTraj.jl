using ...Targets
using Base.Threads
using Interpolations
using PyCall
using CUDA

"An interface of molecular calculation using PySCF."
mutable struct PySCFMolecularCalculator <: MolecularCalculatorBase
    "The molecule to be calculated."
    mol::Molecule;
    "The basis function used for calculation."
    basis::String;
    "The molecular orbital from which the electron is ionized. Candidates are `:HOMO`, `:HOMO-1`, `:HOMO-2`."
    orbitIdx::Symbol;

    # private variables.
    "PySCF library."
    _pyscf;
    "The PySCF molecule object."
    _pymol;
    "The PySCF computation task object."
    _pytask;


    """
    Initializes an instance of `PySCFMolecularCalculator` with given parameter.
    # Parameters
    - `mol::Molecule`   : The molecule to be calculated.
    - `basis::String`   : Basis set used for calculation (default `pcseg-3`).
    - `orbitIdx::Symbol`: The molecular orbital from which the electron is ionized (default `:HOMO`). Candidates are `:HOMO`, `:HOMO-1`, `:HOMO-2`.
    """
    function PySCFMolecularCalculator(mol::Molecule, basis::String="pcseg-3", orbitIdx::Symbol=:HOMO)
        if ! orbitIdx in [:HOMO, :HOMO-1, :HOMO-2]
            error("[PySCFMolecularCalculator] orbitIdx $orbitIdx not supported.")
        end
        mc::PySCFMolecularCalculator = new(mol,basis,orbitIdx)
        try
            mc._pyscf  = pyimport("pyscf")
            mc._pymol  = mc._pyscf.gto.M(atom=exportMolAtomInfo(mol), charge=MolCharge(mol), basis=basis)
            mc._pytask = mc._pyscf.scf.RHF(mc._pymol)
            @info "[PySCFMolecularCalculator] Running molecular calculation..."
            time = @elapsed mc._pytask.run()
            @info "Finished initialization [taking $time second(s)]."
        catch
            @error "[PySCFMolecularCalculator] Encountered error when calling pyscf."
            rethrow()
        end
        return mc
    end
end

function calcStructFactor(; molCalc::PySCFMolecularCalculator,
                            grid_rNum::Int  = 150,
                            grid_rMax::Real = 10.,
                            grid_θNum::Int  = 90,
                            grid_ϕNum::Int  = 90,
                            EulerGrid_βNum::Int,
                            EulerGrid_γNum::Int,
                            GPU_enabled::Bool = false)
    # == PROCEDURE ==
    # 0. Obtain the wavefunction and coefficients (finished in the initialization).
    # 1. Calculate the Hartree-Fock potential.
    #   1-1. Calc. the nuclear potential.
    #   1-2. Calc. the Coulomb part of inter-electron interaction potential.
    #   1-3. Calc. the exchange part of inter-electron interaction potential.
    # 2. Calculate the integral for each Euler angle parameter (β,γ) (α is omitted) to get structure factor.


    mol = molCalc._pymol    # storing the molecule's info and the basis's info.
    task = molCalc._pytask  # storing the calculation result.

    #* Defining the spherical grids.
    grid_rMin = 0.01
    r_grid = range(start=grid_rMin, stop=grid_rMax, length=grid_rNum)
    θ_grid = range(start=0., stop=π, length=grid_θNum)
    ϕ_grid = range(start=0., stop=2π, length=grid_ϕNum)
    # converting to cartesian grid points.
    # x = r sinθ cosϕ, y = r sinθ sinϕ, z = r cosθ.
    pt_x = zeros(grid_rNum*grid_θNum*grid_ϕNum)
    pt_y = zeros(grid_rNum*grid_θNum*grid_ϕNum)
    pt_z = zeros(grid_rNum*grid_θNum*grid_ϕNum)
    @threads for i in 1:grid_rNum*grid_θNum*grid_ϕNum
        iϕ = i % grid_ϕNum + 1
        iθ = floor(Int, i/grid_ϕNum) % grid_θNum + 1
        ir = floor(Int, i/grid_ϕNum/grid_θNum) + 1
        r = r_grid[ir]
        θ = θ_grid[iθ]
        ϕ = ϕ_grid[iϕ]
        pt_x[i] = r*sin(θ)*cos(ϕ)
        pt_y[i] = r*sin(θ)*sin(ϕ)
        pt_z[i] = r*cos(θ)
    end

    #* 1. Calculate the HF potential.
    

end
