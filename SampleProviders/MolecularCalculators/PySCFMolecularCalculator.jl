using ...Targets
using PyCall

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
        catch
            @error "[PySCFMolecularCalculator] Encountered error when calling pyscf."
            rethrow()
        end
        return mc
    end
end
