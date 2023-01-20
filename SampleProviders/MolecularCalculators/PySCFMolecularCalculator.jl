using ...Targets
using Base.Threads
using SpecialFunctions
using SphericalHarmonics
using Folds
using Einsum
using PyCall
using HDF5
using Dates
using ProgressMeter

"An interface of molecular calculation using PySCF."
mutable struct PySCFMolecularCalculator <: MolecularCalculatorBase
    "The molecule to be calculated."
    mol::Molecule;
    "The basis function used for calculation."
    basis::String;
    "The molecular orbital from which the electron is ionized."
    orbitIdx::String;

    "PySCF library."
    _pyscf;
    "The PySCF molecule object."
    _pymol;
    "The PySCF computation task object."
    _pytask;

    "HOMO energy."
    HOMO_energy;
    "Dipole momentum vector in the MF."
    dip_momentum;

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
        @info "[PySCFMolecularCalculator] Running molecular calculation..."
        try
            mc._pyscf  = pyimport("pyscf")
            mc._pymol  = mc._pyscf.gto.M(atom=exportMolAtomInfo(mol), charge=MolCharge(mol), basis=basis)
            mc._pytask = mc._pyscf.scf.RHF(mc._pymol)
            mc._pytask.chkfile = nothing
            time = @elapsed mc._pytask.run()
        catch
            @error "[PySCFMolecularCalculator] Encountered error when calling pyscf."
            rethrow()
        end
        if ! mc._pytask.converged
            error("[PySCFMolecularCalculator] SCF calculation unable to converge.")
        end
        mc.HOMO_energy = mc._pytask.mo_energy[HOMOIndex(mc)]
        mc.dip_momentum = mc._pytask.dip_moment(unit="AU")
        @info "Finished initialization [taking $time second(s)]."
        return mc
    end
end

"Gets the index of the HOMO."
function HOMOIndex(mc::PySCFMolecularCalculator)
    task     = mc._pytask
    mo_occ   = Int.(task.mo_occ)
    num_AO   = size(task.mo_coeff,1)
    idx = 1
    for i = 1:num_AO
        if mo_occ[i]==0
            idx = i; break
        end
    end
    return idx-1    # idx is LUMO, while idx-1 is HOMO
end

"Gets the energy level of the molecule's HOMO (in a.u.)."
function HOMOEnergy(mc::PySCFMolecularCalculator)
    return mc.HOMO_energy
end

"""
Gets the energy level of a specific molecular orbital (MO) of the molecule (in a.u.).
- `orbitIdx::Int`   : Index of the orbital.
"""
function EnergyLevel(mc::PySCFMolecularCalculator, orbitIdx::Int)
    return mc._pytask.mo_energy[orbitIdx]
end

"Gets the energy levels of all the molecular orbitals (MO) of the molecule (in a.u.)."
function EnergyLevels(mc::PySCFMolecularCalculator)
    return mc._pytask.mo_energy
end

"Gets the permanent dipole momentum vector of the molecule in the molecular frame (MF) (in a.u.)."
function DipoleMomentum(mc::PySCFMolecularCalculator)
    return mc.dip_momentum
end

"""
Calculates the DATA used in structure factor calculation in WFAT of the given molecule.

# Parameters
- `molCalc`     : The molecular calculator.
- `grid_rNum`   : The number of radial grid (default 200).
- `grid_rMax`   : The maximum radius of the radial grid (default 10.0).
- `grid_θNum`   : The number of angular grid in the θ direction (default 60).
- `grid_ϕNum`   : The number of angular grid in the ϕ direction (default 60).
- `sf_nξMax`    : The maximum number of nξ used in calculation (default 3).
- `sf_mMax`     : The maximum number of |m| used in calculation (default 3).
- `sf_lMax`     : The maximum angular quantum number l used in calculation (default 6).
- `save_path`   : The path to save the data (default "WFAT_StructFactor_CalcData_yyyymmdd_hhmmss.h5").
                   An empty string indicates saving as default filename; `nothing` indicates no saving.
"""
function calcStructFactorData(; molCalc::PySCFMolecularCalculator,
                                grid_rNum::Int  = 200,
                                grid_rMax::Real = 10.,
                                grid_θNum::Int  = 60,
                                grid_ϕNum::Int  = 60,
                                sf_nξMax ::Int = 3,
                                sf_mMax  ::Int = 3,
                                sf_lMax  ::Int = 6,
                                save_path::String = "")
    # == PROCEDURE ==
    # 0. Obtain the coefficients (finished in the initialization).
    # 1. Calculate the effective core potential.
    # 2. Calculate the integrals and save them as output.

    #* Preprocess molecular information

    pyscf = molCalc._pyscf
    pyscf_df = pyimport("pyscf.df")
    pymol = molCalc._pymol   # storing the molecule's info and the basis's info.
    task  = molCalc._pytask  # storing the calculation result.

    mol         = molCalc.mol
    mo_coeff    = task.mo_coeff     # Linear combination coefficients of AO to make up MO.
    mo_occ      = Int.(task.mo_occ) # Ground state electron occupation number of each MO.
    num_atom    = pymol.natm        # Total number of atoms.
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
        error("[PySCFMolecularCalculator] The energy of the selected molecular orbital $(molCalc.orbitIdx) is positive.")
    end
    κ   = sqrt(2Ip)

    #* Define the spherical grids.
    grid_rMin = 0.001
    grid_dr = (grid_rMax-grid_rMin)/(grid_rNum-1)
    r_grid = range(start=grid_rMin, stop=grid_rMax, length=grid_rNum)
    θ_grid = range(start=0., stop=π, length=grid_θNum)
    ϕ_grid = range(start=0., stop=2π, length=grid_ϕNum)
    N = grid_rNum*grid_θNum*grid_ϕNum

    "Returns the spherical grid indices of the given index of the point."
    function ptIdx2sphCoordIdx(i::Int)
        iϕ = (i-1) % grid_ϕNum + 1
        iθ = (ceil(Int, i/grid_ϕNum)-1) % grid_θNum + 1
        ir = ceil(Int, i/grid_ϕNum/grid_θNum)
        return (ir,iθ,iϕ)
    end

    # converting to cartesian grid points.
    # x = r sinθ cosϕ, y = r sinθ sinϕ, z = r cosθ.
    pt_x = zeros(N)
    pt_y = zeros(N)
    pt_z = zeros(N)
    dV   = zeros(N) # dV = r²sinθ drdθdϕ, used in the integration.
    dr   = r_grid[2]-r_grid[1]
    dθ   = θ_grid[2]-θ_grid[1]
    dϕ   = ϕ_grid[2]-ϕ_grid[1]
    @threads for i in 1:N
        ir,iθ,iϕ = ptIdx2sphCoordIdx(i)
        r,θ,ϕ = r_grid[ir],θ_grid[iθ],ϕ_grid[iϕ]
        pt_x[i] = r*sin(θ)*cos(ϕ)
        pt_y[i] = r*sin(θ)*sin(ϕ)
        pt_z[i] = r*cos(θ)
        dV[i]   = r^2*sin(θ)*dr*dθ*dϕ
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
                Vc_ψ0[i] -= atomCharges[iatm]/sqrt((pt_xyz[i,1]-atomCoords[iatm,1])^2+(pt_xyz[i,2]-atomCoords[iatm,2])^2+(pt_xyz[i,3]-atomCoords[iatm,3])^2+1e-5)
            end
        , 1:N)
    end

    #*  1.4 Calculate the inter-electron interaction: Vd & Vex
    batch_size = 5000  # the integral takes huge memory and thus needs to be performed in batches.
    batch_num = ceil(Int, N/batch_size)
    den_mat = zeros(num_AO,num_AO); i,α,β=0,0,0 # den_mat_αβ = ∑_i c[i,α]*c[i,β]
    occupied_mo_coeff = mo_coeff[:,1:HOMO_orbitIdx]
    @einsum den_mat[α,β] := 2 * occupied_mo_coeff[α,i] * occupied_mo_coeff[β,i]

    prog11 = ProgressUnknown(dt=0.2, desc="Calculating the effective potential...", color = :cyan, spinner = true)
    prog12 = Progress(N; dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1)

    # @threads for i in 1:batch_num
    for i in 1:batch_num    #TODO: Multi-threading is disabled because PyCall.jl doesn't support it. Directly calling libcint might be a solution.
        pt_idx = if i < batch_num
            CartesianIndices((((i-1)*batch_size+1): i*batch_size,))
        else
            CartesianIndices((((i-1)*batch_size+1): N,))
        end
        fakemol = pyscf.gto.fakemol_for_charges(pt_xyz[pt_idx,:])
        I = pyscf_df.incore.aux_e2(pymol, fakemol)      # Size: num_AO × num_AO × batch_size, I_αβ = ∫dr' (χα(r'-Rα)*χβ(r'-Rβ))/|r-r'|
        #* Calculate Vd
        Vd = zeros(size(pt_idx,1)); α,β,i_pt = 0,0,0     # α,β are indices of the basis, i_pt is the index of the points.
        @einsum Vd[i_pt] := I[α,β,i_pt] * den_mat[α,β]     # Einstein's summation notation is used.
        Vc_ψ0[pt_idx] .+= Vd            # now it is (Z/r + Vnuc + Vd).
        Vc_ψ0[pt_idx] .*= ψ0[pt_idx]    # now it is (Z/r + Vnuc + Vd) * ψ0.
        #* Calculate Vex_ψ0
        χα = χi[pt_idx,:]   # basis function values on the selected points.
        Vex_ψ0 = zeros(size(pt_idx,1)); α,β,γ,i_pt = 0,0,0,0
        @einsimd Vex_ψ0[i_pt] := -1/2 * den_mat[α,β] * χα[i_pt,α] * orbit_coeff[γ] * I[β,γ,i_pt]
        Vc_ψ0[pt_idx] .+= Vex_ψ0    # finished building Vc_ψ0.

        next!(prog11,spinner=raw"-\|/"); update!(prog12,pt_idx.indices[1][end]);
    end
    finish!(prog11); finish!(prog12); println()

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
    "Normalization coefficient ω_l^ν for radial function R_l^ν of Ω_{lm'}^ν."
    function ω(nξ,m,l,Z,κ)
        F1 = ((-1)^(l+(abs(m)-m)/2+1)) * (2^(l+3/2)) * (κ^(Z/κ-(abs(m)+1)/2-nξ))
        F2 = sqrt(1.0*(2l+1)*fact(l+m)*fact(l-m)*fact(abs(m)+nξ)*fact(nξ)) * fact(l)/fact(2l+1)     # 1.0 to avoid overflow
        F3 = 0  # Factor3 is a sum over k from 0 to min(n,l-|m|).
        for k in 0:min(nξ,l-abs(m))
            F3 += gamma(l+1-Z/κ+nξ-k) / (fact(k)*fact(l-k)*fact(abs(m)+k)*fact(l-abs(m)-k)*fact(nξ-k))
        end
        return F1*F2*F3
    end
    "Radial function R_l^ν of Ω_{lm'}^ν EXCLUDING the normalization constant ω."
    function R_(l,Z,κ,r)
        return (κ*r)^l*exp(-κ*r)*M(l+1-Z/κ,2l+2,2κ*r)
    end

    #*  2.2 Create pre-computation data to accelerate.
    R_precomp_data = zeros(sf_nξMax+1, sf_mMax+1, sf_lMax+1, grid_rNum)
    # gets the R_nml(r) by indexing [n+1,abs(m)+1,l+1,r_idx] for m≥0, as for m<0, times (-1)^m.
    @threads for l in 0:sf_lMax
        R_precomp_data[1,1,l+1,:] = map(r->R_(l,Z,κ,r), r_grid)
        for nξ in 0:sf_nξMax
        for m  in 0:sf_mMax
            R_precomp_data[nξ+1,m+1,l+1,:] = R_precomp_data[1,1,l+1,:] .* ω(nξ,m,l,Z,κ)
        end
        end
    end
    Y_precomp_data = zeros(ComplexF64, sf_lMax+1, sf_lMax+1, grid_θNum, grid_ϕNum)  # for given l, -l ≤ m' ≤ l.
    # obtain the Y_lm'(θ,ϕ) by indexing [l+1,abs(m)+1,θ_idx,ϕ_idx] for m≥0, as for m<0, times exp(-2im*m*ϕ).
    @threads for l in 0:sf_lMax
    for m_ in 0:l
        for iθ in eachindex(θ_grid)
        for iϕ in eachindex(ϕ_grid)
            Y_precomp_data[l+1,m_+1,iθ,iϕ] = SphericalHarmonics.sphericalharmonic(θ_grid[iθ], ϕ_grid[iϕ]; l=l, m=m_)
        end; end
    end; end
    "Utilizes the pre-computed data to calculate the Ω_{lm'}^{ν}, where ν=(nξ,m)."
    @inline function Ωνl_precomp(nξ,m,l,m_, i_pt)
        ir,iθ,iϕ = ptIdx2sphCoordIdx(i_pt)
        abs(m_)>l && return 0.0
        F1 = if m≥0
            R_precomp_data[nξ+1,m+1,l+1,ir]
        else
            R_precomp_data[nξ+1,abs(m)+1,l+1,ir] * (-1)^m
        end
        F2 = if m_≥0
            Y_precomp_data[l+1,m_+1,iθ,iϕ]
        else
            Y_precomp_data[l+1,abs(m_)+1,iθ,iϕ] * exp(-2im*m_*ϕ_grid[iϕ])
        end
        return F1*F2
    end

    #*  2.3 Calculate the integral

    prog21 = ProgressUnknown(dt=0.2, desc="Calculating the integrals...", color = :cyan, spinner = true)
    prog22 = Progress((sf_nξMax+1)*(2*sf_mMax+1)*(sf_lMax+1)^2; dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1)

    # `IntData` would store the final data: The integral I(nξ,m,l,m')=∫Ω(nξ,m,l,m')*Vc_ψ0(r)*dV.
    # nξ=0,1,⋯,nξMax;  m=0,±1,⋯,±mMax;  l=0,1,⋯,lMax;  m'=-l,-l+1,⋯,0,1,⋯,l.
    # Obtain I(nξ,m,l,m') by indexing [nξ+1, m+mMax+1, l+1, m'+l+1]
    IntData = zeros(sf_nξMax+1,2*sf_mMax+1,sf_lMax+1,2*sf_lMax+1)
    Vc_ψ0_dV = Vc_ψ0 .* dV
    for nξ in 0:sf_nξMax
    for m in -sf_mMax:sf_mMax
    for l in 0:sf_lMax
    for m_ in -l:l
        IntData[nξ+1, m+sf_mMax+1, l+1, m_+l+1] = real(Folds.mapreduce(i->conj(Ωνl_precomp(nξ,m,l,m_,i))*Vc_ψ0_dV[i], +, 1:N))
        next!(prog21); next!(prog22)
    end; end; end; end
    finish!(prog21); finish!(prog22); println()

    #* 3. Save the data

    function defaultFileName()
        Y,M,D = yearmonthday(now())
        h,m,s = hour(now()), minute(now()), second(now())
        return "WFAT_INTDATA-$(string(Y,pad=4))$(string(M,pad=2))$(string(D,pad=2))-$(string(h,pad=2))$(string(m,pad=2))$(string(s,pad=2)).h5"
    end
    if isfile(save_path)
        @warn "[PySCFMolecularCalculator] File \"$save_path\" already exists. Saving at \"$(defaultFileName())\"."
        save_path = defaultFileName()
    end
    if save_path == ""
        save_path = defaultFileName()
    end

    save_success = (true,)  # defining it as a tuple object so that code in the inner block can access it.
    file = nothing
    try
        file = h5open(save_path,"w")
        file["molInfo"] = exportMolAtomInfo(mol)
        file["basis"]   = molCalc.basis
        file["molName"] = mol.name
        file["nξMax"]   = sf_nξMax
        file["mMax"]    = sf_mMax
        file["lMax"]    = sf_lMax
        file["IntData"] = IntData
    catch
        @error "[PySCFMolecularCalculator] Failed to save file at \"$save_path\", trying to save at default path \"$(defaultFileName())\". Check if you have permission to write to the destination file."
        save_success = (false,)
    finally
        if ! isnothing(file)
            close(file)
        end
    end
    if save_success[1]
        @info "Output file saved at \"$save_path\"."
        return IntData
    end
    # if failed, try again
    try
        save_path = defaultFileName()
        file = h5open(save_path,"w")
        file["molInfo"] = exportMolAtomInfo(mol)
        file["basis"]   = molCalc.basis
        file["molName"] = mol.name
        file["nξMax"]   = sf_nξMax
        file["mMax"]    = sf_mMax
        file["lMax"]    = sf_lMax
        file["IntData"] = IntData
        save_success = (true,)
    catch
        @error "[PySCFMolecularCalculator] Failed to save at default path \"$(defaultFileName())\", the data would NOT be saved! Check if you have permission to write to the destination file."
        save_success = (false,)
    finally
        if ! isnothing(file)
            close(file)
        end
    end
    if save_success[1]
        @info "Output file saved at \"$save_path\"."
    end
    return IntData
end
