using ...Targets
using Base.Threads
using SpecialFunctions
using CUDA
using Interpolations
using PyCall

"An interface of molecular calculation using PySCF."
mutable struct PySCFMolecularCalculator <: MolecularCalculatorBase
    "The molecule to be calculated."
    mol::Molecule;
    "The basis function used for calculation."
    basis::String;
    "The molecular orbital from which the electron is ionized. Candidates are `:HOMO`, `:HOMO-1`, `:HOMO-2`, `:LUMO`, `:LUMO+1`, `:LUMO+2`."
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
    - `orbitIdx::Symbol`: The molecular orbital from which the electron is ionized (default `:HOMO`). Candidates are `:HOMO`, `:HOMO-1`, `:HOMO-2`, `:LUMO`, `:LUMO+1`, `:LUMO+2`.
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
                            grid_rNum::Int  = 100,
                            grid_rMax::Real = 10.,
                            grid_θNum::Int  = 60,
                            grid_ϕNum::Int  = 60,
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

    #* Preprocess molecular information

    mol   = molCalc.mol
    pymol = molCalc._pymol   # storing the molecule's info and the basis's info.
    task  = molCalc._pytask  # storing the calculation result.

    Num_AO      = pymol.nao  # Number of atomic orbits (aka AO) or Gaussian basis.
    Num_Elec    = pymol.nelectron   # Total number of electrons.
    mo_coeff    = task.mo_coeff     # Linear combination coefficients of AO to make up MO.
    mo_occ      = Int.(task.mo_occ) # Ground state electron occupation number of each MO.
    HOMO_orbitIdx = begin
        idx = 1
        for i = 1:Num_AO
            if mo_occ[i]==0
                idx = i
                break
            end
        end
        idx
    end
    orbitIdx = HOMO_orbitIdx + begin
        if molCalc.orbitIdx == :HOMO
            0
        elseif molCalc.orbitIdx == :HOMO-1
            -1
        elseif molCalc.orbitIdx == :HOMO-2
            -2
        elseif molCalc.orbitIdx == :LUMO
            +1
        elseif molCalc.orbitIdx == :LUMO+1
            +2
        elseif molCalc.orbitIdx == :LUMO+2
            +3
        end
    end
    Ip  = -task.mo_energy[orbitIdx]
    if Ip ≤ 0
        error("[PySCFMolecularCalculator] The energy of the selected molecular orbit $(molCalc.orbitIdx) is positive.")
    end

    #* Define the spherical grids.
    grid_rMin = 0.001
    r_grid = range(start=grid_rMin, stop=grid_rMax, length=grid_rNum)
    θ_grid = range(start=0., stop=π, length=grid_θNum)
    ϕ_grid = range(start=0., stop=2π, length=grid_ϕNum)
    N = grid_rNum*grid_θNum*grid_ϕNum

    "Returns the spherical coordinate of the given index of the given point."
    function ptIdx2sphCoord(i::Int)
        iϕ = i % grid_ϕNum + 1
        iθ = floor(Int, i/grid_ϕNum) % grid_θNum + 1
        ir = floor(Int, i/grid_ϕNum/grid_θNum) + 1
        return (r_grid[ir],θ_grid[iθ],ϕ_grid[iϕ])
    end

    # converting to cartesian grid points.
    # x = r sinθ cosϕ, y = r sinθ sinϕ, z = r cosθ.
    pt_x = zeros(N)
    pt_y = zeros(N)
    pt_z = zeros(N)
    pt_xyz = hcat(pt_x,pt_y,pt_z)
    dV   = zeros(N) # dV = r²sinθ drdθdϕ, used in the integration.
    @threads for i in 1:N
        r,θ,ϕ = ptIdx2sphCoord(i)
        pt_x[i] = r*sin(θ)*cos(ϕ)
        pt_y[i] = r*sin(θ)*sin(ϕ)
        pt_z[i] = r*cos(θ)
        dV[i]   = r^2*sin(θ)*(r_grid[2]-r_grid[1])*(θ_grid[2]-θ_grid[1])*(ϕ_grid[2]-ϕ_grid[1])
    end

    #* 1. Calculate the HF potential.
    V_c = zeros(N)



    #* 2. Calculate the integral

    #*  2.0 Define some special functions
    "Kummer's confluent hypergeometric function M(a,b,z) = ₁F₁(a;b;z)."
    function M(a,b,z)
        S₀, S₁, j = 1, 1+a*z/b, 1
        while abs(S₀-S₁) > 1e-10 || j ≤ 1
            rⱼ = (a+j)/((b+j)*(j+1))
            S₀, S₁ = S₁, S₁+(S₁-S₀)*rⱼ*z
            j += 1
        end
        return S₁
    end
    "Factorial."
    function fact(N::Int)
        F = 1
        n = 2
        while n≤N
            F*=n
            n+=1
        end
        return F
    end
    "Normalization coefficients ω_l^ν for radial function R_l^ν of Ω_{lm'}^ν. (n→n_ξ, m_→m')"
    function ω(n,l,m,m_,Z,κ)
        F1 = ((-1)^(l+(abs(m)-m)/2+1)) * (2^(l+3/2)) * (κ^(Z/κ-(abs(m)+1)/2-n))
        F2 = sqrt((2l+1)*fact(l+m)*fact(l-m)*fact(abs(m)+n)*fact(n)) * fact(l)/fact(2l+1)
        F3 = 0  # Factor3 is a sum over k from 0 to min(n,l-|m|).
        for k in 1:min(n,l-abs(m))
            F3 += gamma(l+1-Z/κ+n-k) / (fact(k)*fact(l-k)*fact(abs(m)+k)*fact(l-abs(m)-k)*fact(n-k))
        end
        return F1*F2*F3
    end
    "Radial function R_l^ν of Ω_{lm'}^ν. (n→n_ξ, m_→m')"
    function R()
        
    end

    χi = pymol.eval_gto("GTOval",pt_xyz)    # Size: N×Num_AO. Wavefunction of all AOs by calling eval_gto.
    orbit_coeff = (@view mo_coeff[orbitIdx,:])'    # Select the coefficients related to the interested MO.
    if GPU_enabled
        χi          = CuArray(χi)
        orbit_coeff = CuArray(orbit_coeff)
    end
    ψ0 = χi * orbit_coeff   # Calculate wavefunction of the interested MO (matmul operation).


end
