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
    "The molecular orbital from which the electron is ionized."
    orbitIdx::String;

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
    - `orbitIdx::Symbol`: The molecular orbital from which the electron is ionized (default `"HOMO"`). Candidates are `"HOMO"`, `"HOMO-1"`, `"HOMO-2"`, `"LUMO"`, `"LUMO+1"`, `"LUMO+2"`.
    """
    function PySCFMolecularCalculator(mol::Molecule, basis::String="pcseg-3", orbitIdx::String="HOMO")
        if ! (orbitIdx in ["HOMO", "HOMO-1", "HOMO-2", "LUMO", "LUMO+1", "LUMO+2"])
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

"""
Calculates the structure factors in WFAT of the given molecule.
"""
function calcStructFactor(; molCalc::PySCFMolecularCalculator,
                            grid_rNum::Int  = 100,
                            grid_rMax::Real = 10.,
                            grid_θNum::Int  = 60,
                            grid_ϕNum::Int  = 60,
                            EulerGrid_βNum::Int = 30,
                            EulerGrid_γNum::Int = 30,
                            int_lMax  ::Int = 10,
                            int_npmMax::Int = 6)
    # == PROCEDURE ==
    # 0. Obtain the wavefunction and coefficients (finished in the initialization).
    # 1. Calculate the effective core potential.
    # 2. Calculate the integral for each Euler angle parameter (β,γ) (α is omitted) to get structure factor.

    #* Preprocess molecular information

    mol   = molCalc.mol
    pymol = molCalc._pymol   # storing the molecule's info and the basis's info.
    task  = molCalc._pytask  # storing the calculation result.

    mo_coeff    = task.mo_coeff     # Linear combination coefficients of AO to make up MO.
    mo_occ      = Int.(task.mo_occ) # Ground state electron occupation number of each MO.
    num_elec    = pymol.nelectron   # Total number of electrons.
    num_AO      = size(mo_coeff,1)  # Number of atomic orbits (aka AO) or Gaussian basis. (pymol.nao doesn't return an interger!)
    HOMO_orbitIdx = begin
        idx = 1
        for i = 1:num_AO
            if mo_occ[i]==0
                idx = i
                break
            end
        end
        idx-1
    end
    orbitIdx = HOMO_orbitIdx + begin
        if molCalc.orbitIdx == "HOMO"
            0
        elseif molCalc.orbitIdx == "HOMO-1"
            -1
        elseif molCalc.orbitIdx == "HOMO-2"
            -2
        elseif molCalc.orbitIdx == "LUMO"
            +1
        elseif molCalc.orbitIdx == "LUMO+1"
            +2
        elseif molCalc.orbitIdx == "LUMO+2"
            +3
        end
    end

    Z   = MolCharge(molCalc.mol)+1
    Ip  = -task.mo_energy[orbitIdx]
    if Ip ≤ 0
        error("[PySCFMolecularCalculator] The energy of the selected molecular orbit $(molCalc.orbitIdx) is positive.")
    end
    κ   = sqrt(2Ip)

    #* Define the spherical grids.
    grid_rMin = 0.001
    r_grid = range(start=grid_rMin, stop=grid_rMax, length=grid_rNum)
    θ_grid = range(start=0., stop=π, length=grid_θNum)
    ϕ_grid = range(start=0., stop=2π, length=grid_ϕNum)
    N = grid_rNum*grid_θNum*grid_ϕNum

    "Returns the spherical coordinate of the given index of the given point."
    function ptIdx2sphCoord(i::Int)
        iϕ = (i-1) % grid_ϕNum + 1
        iθ = (ceil(Int, i/grid_ϕNum)-1) % grid_θNum + 1
        ir = ceil(Int, i/grid_ϕNum/grid_θNum)
        return (r_grid[ir],θ_grid[iθ],ϕ_grid[iϕ])
    end

    # converting to cartesian grid points.
    # x = r sinθ cosϕ, y = r sinθ sinϕ, z = r cosθ.
    pt_x = zeros(N)
    pt_y = zeros(N)
    pt_z = zeros(N)
    dV   = zeros(N) # dV = r²sinθ drdθdϕ, used in the integration.
    @threads for i in 1:N
        r,θ,ϕ = ptIdx2sphCoord(i)
        pt_x[i] = r*sin(θ)*cos(ϕ)
        pt_y[i] = r*sin(θ)*sin(ϕ)
        pt_z[i] = r*cos(θ)
        dV[i]   = r^2*sin(θ)*(r_grid[2]-r_grid[1])*(θ_grid[2]-θ_grid[1])*(ϕ_grid[2]-ϕ_grid[1])
    end
    pt_xyz = hcat(pt_x,pt_y,pt_z)

    #* 1. Calculate the effective core potential.

    χi = pymol.eval_gto("GTOval",pt_xyz)    # Size: N×Num_AO. Wavefunction of all AOs by calling eval_gto.
    V_HF_coeff = task.get_veff()[:,orbitIdx]
    V_c = χi * V_HF_coeff           # now it is the HF potential.
    @threads for ir in 1:grid_rNum  # add Z/r and now it is the complete core potential.
        V_c[(ir-1)*grid_ϕNum*grid_θNum+1:ir*grid_ϕNum*grid_θNum] .+= Z / r_grid[ir]
    end

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
        F = 1; n = 2
        while n≤N
            F*=n; n+=1
        end
        return F
    end
    "Normalization coefficient ω_l^ν for radial function R_l^ν of Ω_{lm'}^ν. (n→n_ξ)"
    function ω(n,l,m,Z,κ)
        F1 = ((-1)^(l+(abs(m)-m)/2+1)) * (2^(l+3/2)) * (κ^(Z/κ-(abs(m)+1)/2-n))
        F2 = sqrt(1.0*(2l+1)*fact(l+m)*fact(l-m)*fact(abs(m)+n)*fact(n)) * fact(l)/fact(2l+1)     # 1.0 to avoid overflow
        F3 = 0  # Factor3 is a sum over k from 0 to min(n,l-|m|).
        for k in 1:min(n,l-abs(m))
            F3 += gamma(l+1-Z/κ+n-k) / (fact(k)*fact(l-k)*fact(abs(m)+k)*fact(l-abs(m)-k)*fact(n-k))
        end
        return F1*F2*F3
    end
    "Radial function R_l^ν of Ω_{lm'}^ν EXCLUDING the normalization constant ω. (n→n_ξ)"
    function R_(l,Z,κ,r)
        return (κ*r)^l*exp(-κ*r)*M(l+1-Z/κ,2l+2,2κ*r)
    end

    # create pre-computation data to accelerate.
    R_precomp = zeros(int_npmMax,int_lMax+1,int_npmMax,grid_rNum)  # gets the R_nlm(r) by calling R_precomp[n,l+1,abs(m)+1,r_idx] for m≥0, as for m<0, times (-1)^m.
    @threads for l in 0:int_lMax
        R_precomp[1,l+1,1,:] = map(r->R_(l,Z,κ,r), r_grid)
        for n in 1:int_npmMax
        for m in 0:min(int_npmMax-n, l)
            R_precomp[n,l+1,m+1,:] = R_precomp[1,l+1,1,:] .* ω(n,l,m,Z,κ)
        end
        end
    end

    #*  2.1 Calculate the wavefunction
    orbit_coeff = @view mo_coeff[:,orbitIdx]    # Select the coefficients related to the interested MO.
    ψ0 = χi * orbit_coeff   # Calculate wavefunction of the interested MO (matmul operation).

end
