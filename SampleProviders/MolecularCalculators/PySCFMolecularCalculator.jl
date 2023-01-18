using ...Targets
using Base.Threads
using SpecialFunctions
using SphericalHarmonics
using Rotations
using Interpolations
using Folds
using TensorOperations
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

    pyscf = molCalc._pyscf
    pyscf_df = pyimport("pyscf.df")
    mol   = molCalc.mol
    pymol = molCalc._pymol   # storing the molecule's info and the basis's info.
    task  = molCalc._pytask  # storing the calculation result.

    mo_coeff    = task.mo_coeff     # Linear combination coefficients of AO to make up MO.
    mo_occ      = Int.(task.mo_occ) # Ground state electron occupation number of each MO.
    den_mat     = task.make_rdm1()  # Density matrix.
    num_atom    = pymol.natm        # Total number of atoms.
    num_elec    = pymol.nelectron   # Total number of electrons.
    num_AO      = size(mo_coeff,1)  # Number of atomic orbits (aka AO) or Gaussian basis. (pymol.nao doesn't return an interger!)
    HOMO_orbitIdx = begin
        idx = 1
        for i = 1:num_AO
            if mo_occ[i]==0
                idx = i; break
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
    grid_dr = (grid_rMax-grid_rMin)/(grid_rNum-1)
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

    #* 1. Calculate the effective core potential
    # Vc is composed of:
    # 1. Asymptotic Coulomb potential Z/r
    # 2. Nuclear potential Vnuc
    # 3. Inter-electron interaction: direct part Vd
    # 4. Inter-electron interaction: exchange part Vex
    # The direct expression of Vex couldn't be obtained, but Vex*ψ0 is obtainable.
    # In the following scheme, Z/r, Vnuc, Vd would be first calculated and directly added to Vc_ψ0,
    # and then would be multiplied by ψ0.
    # Vc_ψ0 = (Z/r + Vnuc + Vd) * ψ0 + Vex_ψ0.
    Vc_ψ0 = zeros(N)

    #*  1.1 Calculate the wavefunction
    χi = pymol.eval_gto("GTOval",pt_xyz)        # Size: N×Num_AO. Wavefunction of all AOs by calling eval_gto.
    orbit_coeff = @view mo_coeff[:,orbitIdx]    # Select the coefficients related to the interested MO.
    ψ0 = χi * orbit_coeff                       # Calculate wavefunction of the interested MO (matmul operation).
    ψ0 ./= sqrt(Folds.mapreduce(i->abs2(ψ0[i])*dV[i], +, 1:N))    # Normalization

    #*  1.2 Calculate the asymptotic Coulomb potential Z/r
    @threads for ir in 1:grid_rNum
        Vc_ψ0[(ir-1)*grid_ϕNum*grid_θNum+1:ir*grid_ϕNum*grid_θNum] .+= Z / r_grid[ir]
    end

    #*  1.3 Calculate the nuclear potential Vnuc
    atomCharges = pymol.atom_charges()
    atomCoords  = pymol.atom_coords()   # Size: Num_Atoms×3
    for iatm in 1:num_atom
        Folds.map(
            function (i)
                Vc_ψ0[i] -= atomCharges[iatm]/sqrt((pt_xyz[i,1]-atomCoords[iatm,1])^2+(pt_xyz[i,2]-atomCoords[iatm,2])^2+(pt_xyz[i,3]-atomCoords[iatm,3])^2)
            end
        , 1:N)
    end

    #*  1.4 Calculate the inter-electron interaction: Vd & Vex
    batch_size = 5000  # the integral takes huge memory and thus needs to be performed in batches.
    batch_num = ceil(Int, N/batch_size)
    Vd = nothing
    @threads for i in 1:batch_num
        pt_idx = if i < batch_num
            CartesianIndices(((i-1)*batch_size+1:i*batch_size))
        else
            CartesianIndices(((i-1)*batch_size+1:N))
        end
        fakemol = pyscf.gto.fakemol_for_charges(pt_xyz[pt_idx,:])
        I = pyscf.df.incore.aux_e2(pymol, fakemol)      # Size: num_AO × num_AO × batch_size, I_ab = ∫dr' (χa(r'-Ra)*χb(r'-Rb))/|r-r'|
        @tensoropt Vd[3] := I[1,2,3] * den_mat[1,2]     # Vd[k] := I[i,j,k] * den_mat[i,j]

    end

    #* 2. Calculate the integral
    #*  2.1 Define some special functions
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
        for k in 0:min(n,l-abs(m))
            F3 += gamma(l+1-Z/κ+n-k) / (fact(k)*fact(l-k)*fact(abs(m)+k)*fact(l-abs(m)-k)*fact(n-k))
        end
        return F1*F2*F3
    end
    "Radial function R_l^ν of Ω_{lm'}^ν EXCLUDING the normalization constant ω. (n→n_ξ)"
    function R_(l,Z,κ,r)
        return (κ*r)^l*exp(-κ*r)*M(l+1-Z/κ,2l+2,2κ*r)
    end

    #*  2.2 Create pre-computation data to accelerate.
    R_precomp_data = zeros(int_npmMax+1,int_lMax+1,int_npmMax+1,grid_rNum)  # gets the R_nlm(r) by calling R_precomp[n+1,l+1,abs(m)+1,r_idx] for m≥0, as for m<0, times (-1)^m.
    @threads for l in 0:int_lMax
        R_precomp_data[1,l+1,1,:] = map(r->R_(l,Z,κ,r), r_grid)
        for n in 0:int_npmMax
        for m in 0:min(int_npmMax-n, l)
            R_precomp_data[n+1,l+1,m+1,:] = R_precomp_data[1,l+1,1,:] .* ω(n,l,m,Z,κ)
        end
        end
    end
    angularGrid = [(θ_grid[iθ],ϕ_grid[iϕ]) for iθ in 1:grid_θNum, iϕ in 1:grid_ϕNum] # obtain the interpolation of Y_lm by indexing (l+1,abs(m)+1) for m≥0, as for m<0, times exp(-2im*m*ϕ).
    Y_precomp_data = Array{AbstractInterpolation}(undef, int_lMax+1, int_npmMax+1)
    for l in 0:int_lMax
        for m in 0:min(int_npmMax, l)
            Y_precomp_data[l+1,m+1] = linear_interpolation((θ_grid,ϕ_grid), map(ang->SphericalHarmonics.sphericalharmonic(ang...;l=l,m=m), angularGrid))
        end
    end
    function Ων_precomp(n,m,(r,θ,ϕ))
        sum = 0
        abs(m)>int_lMax && return 0.0
        for l in abs(m):int_lMax
            if m≥0
                sum += R_precomp_data[n+1,l+1,m+1,round(Int,(r-grid_rMin)/grid_dr)+1]                 * Y_precomp_data[l+1,m+1](θ,ϕ)
            else
                sum += R_precomp_data[n+1,l+1,abs(m)+1,round(Int,(r-grid_rMin)/grid_dr)+1] * (-1)^m   * Y_precomp_data[l+1,abs(m)+1](θ,ϕ) * exp(-2im*m*ϕ)
            end
        end
        return sum
    end

    #*  2.3 Calculate the integral
    "Converts the cartesian coordinate to spherical one."
    function cart2sph(vec)
        x = vec[1]
        y = vec[2]
        z = vec[3]
        r = sqrt(x^2+y^2+z^2)
        θ = acos(z/r)
        ϕ = atan(y,x)
        (ϕ<0) && (ϕ+=2π)
        return [r,θ,ϕ]
    end
    βlist = 0:π/18:π
    glist = zeros(size(βlist))
    for iβ in eachindex(βlist)
        rot = inv(Rotations.RotZXZ(0,βlist[iβ],0))
        glist[iβ] = Folds.mapreduce(i->conj(Ων_precomp(0,0,Tuple(cart2sph(rot*pt_xyz[i,1:3]))))*Vc_ψ0[i]*dV[i], +, 1:N)
        @info "β=$(βlist[iβ]), g=$(glist[iβ])"
    end
    return βlist,glist
end
