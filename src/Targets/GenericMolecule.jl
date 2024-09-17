
"""
    struct GenericMolecule <: MoleculeBase <: Target

Represents a generic molecule.
"""
mutable struct GenericMolecule <: MoleculeBase

    "Molecular calculator that calculates energy, structure factors ... of the molecule."
    mol_calc

    #* primitive properties

    atoms;
    atom_coords;
    charge::Integer;
    spin;
    name::String;

    #* properties that would be stored in data upon related calculation.

    # energy levels
    energy_data_available::Bool;
    energy_levels;
    orbit_occ;

    # WFAT IntData
    wfat_data_available::Bool;
    wfat_indices::Set;
    wfat_intdata::Dict;
    wfat_μ::Dict;

    # Asymptotic Coeff
    asymp_coeff_available::Bool;
    asymp_coeff_indices::Set;
    asymp_coeff::Dict;

    #* properties that would not be stored in data

    rot_α;
    rot_β;
    rot_γ;
end

#* fresh init without data.
"""
    GenericMolecule(atoms, atom_coords [,charge=0] [,spin=0] [,name] [,rot_α=0.0] [,rot_β=0.0] [,rot_γ=0.0])

Initializes a new `GenericMolecule` with given parameters.

## Parameters
- `atoms`       : Atoms in the molecule, stored as a `Vector` of `String`.
- `atom_coords` : Atoms' coordinates in the molecule (numerically in Å or a `Unitful.Quantity`), stored as a N×3 `Matrix`.
- `charge`      : Total charge of the molecule (ion) (*optional, default `0`*).
- `spin`        : Total spin of the molecule (*optional, default `0`*). Note that each unpaired electron contributes 1/2.
- `name`        : Name of the molecule (*optional*).
- `rot_α`,`rot_β`,`rot_γ`   : Euler angles (ZYZ convention) specifying the molecule's orientation (numerically in radian or a `Unitful.Quantity`) (*optional, default `0`*).

## Example
```
julia> m = GenericMolecule(atoms=["H","H"], atom_coords=[0.0 0.0 -0.375; 0.0 0.0 0.375], name="Hydrogen")
[GenericMolecule] Hydrogen

julia> using SemiclassicalSFI.Units

julia> m = GenericMolecule(atoms=["H","H"], atom_coords=[0.0 0.0 -0.375; 0.0 0.0 0.375]*Å, name="Hydrogen", rot_β=90°)
[GenericMolecule] Hydrogen, αβγ=(0.0°,90.0°,0.0°)
```
"""
function GenericMolecule(;atoms::Vector,atom_coords::Matrix,charge::Integer=0,spin=0,name::String="[NA]",rot_α=0.,rot_β=0.,rot_γ=0.)
    @assert eltype(atoms) <: String   "[GenericMolecule] Element type of `atoms` must be String."
    @assert atom_coords isa Matrix && ndims(atom_coords)==2 && size(atom_coords,2)==3 && size(atom_coords,1)==size(atoms,1)   "[GenericMolecule] `atom_coords` should be a Matrix of size N×3."
    @assert spin>=0 "[GenericMolecule] `spin` must be non-negative."
    # unit transformation
    (eltype(atom_coords)<:Quantity) && (atom_coords=map(q->uconvert(u"Å", q).val, atom_coords))
    (typeof(rot_α)<:Quantity) && (rot_α=uconvert(u"rad",rot_α).val)
    (typeof(rot_β)<:Quantity) && (rot_β=uconvert(u"rad",rot_β).val)
    (typeof(rot_γ)<:Quantity) && (rot_γ=uconvert(u"rad",rot_γ).val)
    mol = GenericMolecule(
        nothing,    # mol_calc
        atoms, atom_coords, charge, spin, name,
        false, Float64[], Float64[],        # energy_data
        false, Set(), Dict(), Dict(),   # wfat_data
        false, Set(), Dict(),           # asymp_coeff
        rot_α,rot_β,rot_γ)
    return mol
end

#* init from data.
"""
    LoadMolecule(ext_data_path; [rot_α=0.0] [,rot_β=0.0] [,rot_γ=0.0])

Initializes a new `GenericMolecule` with the data stored in `ext_data_path`.

## Parameters
- `ext_data_path`           : Path to the molecule's data stored externally.
- `rot_α`,`rot_β`,`rot_γ`   : Euler angles (ZYZ convention) specifying the molecule's orientation (numerically in radian or a `Unitful.Quantity`) (*optional, default `0`*).
"""
function LoadMolecule(ext_data_path::String; rot_α=0.,rot_β=0.,rot_γ=0.)
    file = jldopen(ext_data_path ,"r")
    @unpack atoms, atom_coords, charge, spin, name,
        energy_levels, orbit_occ,
        wfat_indices, wfat_intdata, wfat_μ,
        asymp_coeff_indices, asymp_coeff = file
    energy_data_available = !isempty(energy_levels)
    wfat_data_available = !isempty(wfat_indices)
    asymp_coeff_available = !isempty(asymp_coeff_indices)
    close(file)
    (typeof(rot_α)<:Quantity) && (rot_α=uconvert(u"rad",rot_α).val)
    (typeof(rot_β)<:Quantity) && (rot_β=uconvert(u"rad",rot_β).val)
    (typeof(rot_γ)<:Quantity) && (rot_γ=uconvert(u"rad",rot_γ).val)
    return GenericMolecule(
        nothing,
        atoms, atom_coords, charge, spin, name,
        energy_data_available, energy_levels, orbit_occ,
        wfat_data_available, wfat_indices, wfat_intdata, wfat_μ,
        asymp_coeff_available, asymp_coeff_indices, asymp_coeff,
        rot_α,rot_β,rot_γ)
end


#* Molecule's specific properties & methods

MolAtoms(mol::GenericMolecule) = mol.atoms
MolAtomCoords(mol::GenericMolecule) = mol.atom_coords
MolCharge(mol::GenericMolecule) = mol.charge
MolSpin(mol::GenericMolecule) = mol.spin

function MolEnergyDataAvailable(mol::GenericMolecule)
    return mol.energy_data_available
end

function MolEnergyLevels(mol::GenericMolecule, spin::Integer=0)
    if ! mol.energy_data_available
        error("[GenericMolecule] The energy data is not available, calculate first.")
    end
    if spin == 0 || mol.spin == 0
        return mol.energy_levels
    else
        if spin == 1
            mol.energy_levels[1,:]
        else # spin == 2
            mol.energy_levels[2,:]
        end
    end
end

function MolEnergyLevel(mol::GenericMolecule, orbit_ridx)
    if ! mol.energy_data_available
        error("[GenericMolecule] The energy data is not available, calculate first.")
    end
    if mol.spin != 0
        @assert orbit_ridx isa Tuple{Int,Int} && orbit_ridx[1] in (1,2) "[GenericMolecule] For an open-shell molecule, `orbit_ridx` should be a two-element tuple `(spin, ridx)`."
        mol.energy_levels[orbit_ridx[1],_get_HOMO_idx(mol.orbit_occ[orbit_ridx[1],:])+orbit_ridx[2]]
    else
        @assert isinteger(orbit_ridx) "[GenericMolecule] For a closed-shell molecule, `orbit_ridx` should be an integer."
        mol.energy_levels[_get_HOMO_idx(mol.orbit_occ)+orbit_ridx]
    end
end

function MolOrbitalOccupation(mol::GenericMolecule, spin::Integer=0)
    if ! mol.energy_data_available
        error("[GenericMolecule] The energy data is not available, calculate first.")
    end
    if spin == 0 || mol.spin == 0
        return mol.orbit_occ
    else
        if spin == 1
            mol.orbit_occ[1,:]
        else # spin == 2
            mol.orbit_occ[2,:]
        end
    end
end

function MolHOMOEnergy(mol::GenericMolecule, spin::Integer=0)
    if ! mol.energy_data_available
        error("[GenericMolecule] The energy data is not available, calculate first.")
    end
    if spin == 0
        return mol.energy_levels[_get_HOMO_idx(mol.orbit_occ)]
    else
        if spin == 1
            return mol.energy_levels[1,_get_HOMO_idx(mol.orbit_occ[1,:])]
        elseif spin == 2
            return mol.energy_levels[2,_get_HOMO_idx(mol.orbit_occ[2,:])]
        else
            return [mol.energy_levels[1,_get_HOMO_idx(mol.orbit_occ[1,:])], mol.energy_levels[2,_get_HOMO_idx(mol.orbit_occ[2,:])]]
        end
    end
end
function _get_HOMO_idx(orbit_occ)
    # gets the index of HOMO according to orbit_occ
    findlast(!iszero, orbit_occ)
end
function _get_LUMO_idx(orbit_occ)
    # gets the index of LUMO according to orbit_occ
    findfirst(iszero, orbit_occ)
end

function MolWFATAvailableIndices(mol::GenericMolecule)
    return if mol.wfat_data_available
        mol.wfat_indices
    else
        Set()
    end
end

function MolWFATData(mol::GenericMolecule, orbit_ridx)
    if ! mol.energy_data_available || ! (mol.wfat_data_available || (orbit_ridx in mol.wfat_indices))
        error("[GenericMolecule] The WFAT data is not available, calculate first.")
    end
    return mol.wfat_μ[orbit_ridx], mol.wfat_intdata[orbit_ridx]
end

function MolWFATStructureFactor_G(mol::GenericMolecule, orbit_ridx, nξ::Integer, m::Integer, β::Real, γ::Real)
    if ! mol.energy_data_available || ! mol.wfat_data_available || !(orbit_ridx in mol.wfat_indices)
        error("[GenericMolecule] The WFAT data is not available, calculate first.")
    end
    @assert nξ≥0 "[GenericMolecule] nξ should be non-negative."
    intdata = mol.wfat_intdata[orbit_ridx]
    nξMax = size(intdata,1) - 1
    mMax = round(Int,(size(intdata,2)-1)/2)
    lMax = size(intdata,3) - 1
    if nξ>nξMax
        @error "[GenericMolecule] The given nξ=$nξ is larger than the maximum value $nξMax, zero value would be returned."
        return 0.0
    end
    if abs(m)>mMax
        @error "[GenericMolecule] The given |m|=$(abs(m)) is larger than the maximum value $mMax, zero value would be returned."
        return 0.0
    end
    @inline μz(β,γ) = (RotZYZ(γ,β,0.0)*mol.wfat_μ[orbit_ridx])[3]   # the result is independent of α

    sum = zero(ComplexF64)
    for l in abs(m):lMax, m_ in -l:l
        sum += intdata[nξ+1,m+mMax+1,l+1,m_+l+1] * wignerdjmn(l,m,m_,β) * exp(-1im*m_*γ)
    end
    return sum * exp(-sqrt(2IonPotential(mol,orbit_ridx))*μz(β,γ))
end
function MolWFATStructureFactor_G(mol::GenericMolecule, orbit_ridx, nξ::Integer, m::Integer, β::AbstractVector{T} where T<:Real, γ::AbstractVector{T} where T<:Real)
    if ! mol.energy_data_available || ! mol.wfat_data_available || !(orbit_ridx in mol.wfat_indices)
        error("[GenericMolecule] The WFAT data is not available, calculate first.")
    end
    @assert size(β,1)==size(γ,1) "[GenericMolecule] Invalid input (β,γ), should be both `Real` values or two `Vector`s of `Real` and of same length."
    @assert nξ≥0 "[GenericMolecule] nξ should be non-negative."
    intdata = mol.wfat_intdata[orbit_ridx]
    nξMax = size(intdata,1) - 1
    mMax = round(Int,(size(intdata,2)-1)/2)
    lMax = size(intdata,3) - 1
    if nξ>nξMax
        @error "[GenericMolecule] The given nξ=$nξ is larger than the maximum value $nξMax, zero value would be returned."
        return 0.0
    end
    if abs(m)>mMax
        @error "[GenericMolecule] The given |m|=$(abs(m)) is larger than the maximum value $mMax, zero value would be returned."
        return 0.0
    end
    @inline μz(β,γ) = (RotZYZ(γ,β,0.0)*mol.wfat_μ[orbit_ridx])[3]   # the result is independent of α

    sum = zeros(ComplexF64,size(β))
    for l in abs(m):lMax, m_ in -l:l
        sum .+= intdata[nξ+1,m+mMax+1,l+1,m_+l+1] * @. wignerdjmn(l,m,m_,β) * exp(-1im*m_*γ)
    end
    κ = sqrt(2IonPotential(mol,orbit_ridx))
    return @. sum * exp(-κ*μz(β,γ))
end

function MolWFATMaxChannels(mol::GenericMolecule, orbit_ridx)
    μ, int_data = MolWFATData(mol, orbit_ridx)
    nξMax = size(int_data,1) - 1
    mMax = round(Int,(size(int_data,2)-1)/2)
    return (nξMax, mMax)
end

function MolAsympCoeffAvailableIndices(mol::GenericMolecule)
    return if mol.asymp_coeff_available
        mol.asymp_coeff_indices
    else
        Set()
    end
end

function MolAsympCoeff(mol::GenericMolecule, orbit_ridx)
    if ! mol.energy_data_available || ! mol.asymp_coeff_available || !(orbit_ridx in mol.asymp_coeff_indices)
        error("[GenericMolecule] The asymptotic coefficient is not available, calculate first.")
    end
    return mol.asymp_coeff[orbit_ridx]
end

function MolAsympCoeff_lMax(mol::GenericMolecule, orbit_ridx)
    return size(MolAsympCoeff(mol, orbit_ridx), 1) - 1
end

MolRotation(mol::GenericMolecule) = (mol.rot_α,mol.rot_β,mol.rot_γ)
function SetMolRotation!(mol::GenericMolecule, α,β,γ)
    (α isa Quantity) && (α=uconvert(u"rad",α).val)
    (β isa Quantity) && (β=uconvert(u"rad",β).val)
    (γ isa Quantity) && (γ=uconvert(u"rad",γ).val)
    mol.rot_α = α; mol.rot_β = β; mol.rot_γ = γ;
    return
end
function SetMolRotation!(mol::GenericMolecule, (α,β,γ))
    SetMolRotation!(mol, α,β,γ)
end

function MolExportAtomInfo(mol::GenericMolecule)
    atom2str(i_atm) = join([String(mol.atoms[i_atm]),mol.atom_coords[i_atm,1:3]], " ")
    return join(map(atom2str, eachindex(mol.atoms)),"; ")
end

# data calculation and operation

function MolInitCalculator!(mol::GenericMolecule, MCType::Type = PySCFMolecularCalculator; kwargs...)
    if !isnothing(mol.mol_calc)
        @warn "[GenericMolecule] Molecule's `MolecularCalculator` is present already, replacing."
    end
    if ! (MCType<:MolecularCalculatorBase)
        error("[GenericMolecule] `MCType`'s type $MCType mismatches `MolecularCalculatorBase`.")
    end
    mol.mol_calc = MCType(;mol=mol, kwargs...)
    mol.energy_data_available = true
    mol.energy_levels = EnergyLevels(mol.mol_calc)
    mol.orbit_occ = OrbitalOccupation(mol.mol_calc)
    return
end

function MolCalcWFATData!(mol::GenericMolecule, orbit_ridx; kwargs...)
    if isnothing(mol.mol_calc)
        error("[GenericMolecule] Molecule's `MolecularCalculator` is not initialized, call `MolInitCalculator!` first.")
    end
    if isnothing(mol.wfat_indices)
        mol.wfat_indices = Set()
        mol.wfat_intdata = Dict()
        mol.wfat_μ = Dict()
    end
    if mol.spin==0
        @assert isinteger(orbit_ridx) "[GenericMolecule] For a closed-shell molecule, `orbit_ridx` should be an integer."
    else
        @assert orbit_ridx isa Tuple{Int,Int} && orbit_ridx[1] in (1,2) "[GenericMolecule] For an open-shell molecule, `orbit_ridx` should be a two-element tuple `(spin, ridx)`."
    end
    mol.wfat_μ[orbit_ridx], mol.wfat_intdata[orbit_ridx] = calc_WFAT_data(;mc=mol.mol_calc, orbit_ridx=orbit_ridx, kwargs...)
    push!(mol.wfat_indices, orbit_ridx)
    mol.wfat_data_available = true
    return
end

function MolCalcAsympCoeff!(mol::GenericMolecule, orbit_ridx; kwargs...)
    if isnothing(mol.mol_calc)
        error("[GenericMolecule] Molecule's `MolecularCalculator` is not initialized, call `MolInitCalculator!` first.")
    end
    if isnothing(mol.asymp_coeff_available)
        mol.asymp_coeff_indices = Set()
        mol.asymp_coeff = Dict()
    end
    if mol.spin==0
        @assert isinteger(orbit_ridx) "[GenericMolecule] For a closed-shell molecule, `orbit_ridx` should be an integer."
    else
        @assert orbit_ridx isa Tuple{Int,Int} && orbit_ridx[1] in (1,2) "[GenericMolecule] For an open-shell molecule, `orbit_ridx` should be a two-element tuple `(spin, ridx)`."
    end
    mol.asymp_coeff[orbit_ridx] = calc_asymp_coeff(; mc=mol.mol_calc, orbit_ridx=orbit_ridx, kwargs...)
    push!(mol.asymp_coeff_indices, orbit_ridx)
    mol.asymp_coeff_available = true
    return
end

function MolSaveDataAs!(mol::GenericMolecule, data_path::String, overwrite::Bool=false)
    function default_filename()
        Y,M,D = yearmonthday(now())
        h,m,s = hour(now()), minute(now()), second(now())
        return "Molecule_$(mol.name)_$(string(Y,pad=4))$(string(M,pad=2))$(string(D,pad=2))-$(string(h,pad=2))$(string(m,pad=2))$(string(s,pad=2)).h5"
    end
    # if destination exists or not specified, saving as default file name.
    default_path = default_filename()
    if isfile(data_path)
        if !overwrite
            @warn "[GenericMolecule] Destination file `$data_path` already exists. Saving at `$default_path`."
            data_path = default_path
        else

        end
    elseif data_path==""
        @warn "[GenericMolecule] Destination file not specified. Saving at `$default_path`."
        data_path = default_path
    end
    atoms       = mol.atoms
    atom_coords = mol.atom_coords
    charge      = mol.charge
    spin        = mol.spin
    name        = mol.name
    energy_levels   = mol.energy_levels
    orbit_occ       = mol.orbit_occ
    wfat_indices    = mol.wfat_indices
    wfat_intdata    = mol.wfat_intdata
    wfat_μ      = mol.wfat_μ
    asymp_coeff_indices = mol.asymp_coeff_indices
    asymp_coeff = mol.asymp_coeff
    function write_jld2(path)
        file = jldopen(path, "w")
        @pack! file = atoms, atom_coords, charge, spin, name,
            energy_levels, orbit_occ,
            wfat_indices, wfat_intdata, wfat_μ,
            asymp_coeff_indices, asymp_coeff
        close(file)
    end
    try
        write_jld2(data_path)
        @info "[GenericMolecule] Data saved for molecule $(mol.name) at `$(data_path)`."
    catch
        @error "[GenericMolecule] Error writing to `$data_path`, trying to save at `$(default_path)`."
        write_jld2(default_path)
        @info "[GenericMolecule] Data saved for molecule $(mol.name) at `$(default_path)`."
    end
    return
end

#* Properties & methods that implement the supertype Target.

function IonPotential(mol::GenericMolecule)
    if ! mol.energy_data_available
        error("[GenericMolecule] The energy data is not available, calculate first.")
    end
    return -maximum(MolHOMOEnergy(mol))
end
function IonPotential(mol::GenericMolecule, orbit_ridx)
    if ! mol.energy_data_available
        error("[GenericMolecule] The energy data is not available, calculate first.")
    end
    return -MolEnergyLevel(mol,orbit_ridx)
end
AsympNuclCharge(mol::GenericMolecule) = mol.charge + 1
TargetName(mol::GenericMolecule) = mol.name
TargetPotential(mol::GenericMolecule) = (x,y,z) -> -(mol.charge+1)*(x^2+y^2+z^2+1.0)^(-0.5)
TargetForce(mol::GenericMolecule) = (x,y,z) -> -(mol.charge+1)*(x^2+y^2+z^2+1.0)^(-1.5) .* (x,y,z)
function TrajectoryFunction(mol::GenericMolecule, dimension::Integer, laserFx::Function, laserFy::Function, phase_method::Symbol; kwargs...)
    Z  = mol.charge+1
    Ip = IonPotential(mol)
    if dimension == 2
        if phase_method == :CTMC
            function traj_dipole_ctmc_2d(u,p,t)
                # tFx, tFy = targetF(u[1],u[2])
                tFx, tFy = -Z*(u[1]^2+u[2]^2+1.0)^(-1.5) .* (u[1],u[2])
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                @SVector [du1,du2,du3,du4]
            end
        elseif phase_method == :QTMC
            function traj_dipole_qtmc_2d(u,p,t)
                tFx, tFy = -Z*(u[1]^2+u[2]^2+1.0)^(-1.5) .* (u[1],u[2])
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]))
                du5 = -(Ip + (du1^2+du2^2)/2 - Z*(u[1]^2+u[2]^2+1.0)^(-0.5))
                @SVector [du1,du2,du3,du4,du5]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_2d(u,p,t)
                tFx, tFy = -Z*(u[1]^2+u[2]^2+1.0)^(-1.5) .* (u[1],u[2])
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]) + (u[1]*tFx+u[2]*tFy))
                du5 = -(Ip + (du1^2+du2^2)/2 - Z*(u[1]^2+u[2]^2+1.0)^(-0.5) + (u[1]*tFx+u[2]*tFy))
                @SVector [du1,du2,du3,du4,du5]
            end
        end
    else # dimension == 3
        if phase_method == :CTMC
            function traj_dipole_ctmc_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                @SVector [du1,du2,du3,du4,du5,du6]
            end
        elseif phase_method == :QTMC
            function traj_dipole_qtmc_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-0.5))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+1.0)^(-0.5) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    end
end

function Base.show(io::IO, mol::GenericMolecule)
    @printf(io, "[GenericMolecule] %s", mol.name)
    if !(mol.rot_α==0 && mol.rot_β==0 && mol.rot_γ==0)
        @printf(io, ", αβγ=(%.1f°,%.1f°,%.1f°)", mol.rot_α*180/π, mol.rot_β*180/π, mol.rot_γ*180/π)
    end
    if mol.asymp_coeff_available
        println(io)
        @printf(io, "Asymp coeff of %s available", join(_MOstring.(sort(mol.asymp_coeff_indices|>collect))," & "))
    end
    if mol.wfat_data_available
        println(io)
        @printf(io, "WFAT data of %s available", join(_MOstring.(sort(mol.wfat_indices|>collect))," & "))
    end
    if mol.energy_data_available
        println(io)
        if mol.spin==0
            spin11="-↿⇂-"
            spin00="----"
            HOMO_idx = _get_HOMO_idx(mol.orbit_occ)
            #        %-3s  %-7s    %7f__ %4s
            #        11122222223333333__4444
            @printf "#          E (Ha)  occp\n"
            @printf "⋮  ⋮         ⋮      ⋮⋮\n"
            @printf "%-3s%-7s%7.3f  %4s\n" HOMO_idx+2 "LUMO+1" mol.energy_levels[HOMO_idx+2] spin00
            @printf "%-3s%-7s%7.3f  %4s\n" HOMO_idx+1 "LUMO"   mol.energy_levels[HOMO_idx+1] spin00
            @printf "%-3s%-7s%7.3f  %4s"   HOMO_idx   "HOMO"   mol.energy_levels[HOMO_idx  ] spin11
            HOMO_idx>1 && (@printf "\n%-3s%-7s%7.3f  %4s" HOMO_idx-1 "HOMO-1" mol.energy_levels[HOMO_idx-1] spin11)
            HOMO_idx>2 && (@printf "\n%-3s%-7s%7.3f  %4s" HOMO_idx-2 "HOMO-2" mol.energy_levels[HOMO_idx-2] spin11)
            HOMO_idx>3 && (@printf "\n%-3s%-7s%7.3f  %4s" HOMO_idx-3 "HOMO-3" mol.energy_levels[HOMO_idx-3] spin11)
            HOMO_idx>4 && (@printf "\n⋮    ⋮        ⋮     ⋮⋮")
        else
            HOMO_alp_idx = _get_HOMO_idx(mol.orbit_occ[1,:])
            HOMO_bet_idx = _get_HOMO_idx(mol.orbit_occ[2,:])
            idx_max = max(HOMO_alp_idx,HOMO_bet_idx)+2
            idx_min = max(min(HOMO_alp_idx,HOMO_bet_idx)-3,1)
            #        %-3s  %-7s    %7f__%4s__    %6f_%-7s
            #        11122222223333333__4444__555555_6666666
            @printf "#          Eα(Ha)  occp  Eβ(Ha)\n"
            @printf "⋮    ⋮        ⋮     ⋮⋮      ⋮     ⋮"
            for i in idx_max:-1:idx_min
                @printf "\n%-3s%-7s%7.3f  %4s  %6.3f %-7s" i _MOstring(i-HOMO_alp_idx) mol.energy_levels[1,i] _MO_occ_string(mol.orbit_occ[1,i],mol.orbit_occ[2,i]) mol.energy_levels[2,i] _MOstring(i-HOMO_bet_idx)
            end
            idx_min>1 && (@printf "\n⋮    ⋮        ⋮     ⋮⋮      ⋮     ⋮")
        end
    end
end

function _MO_occ_string(α,β)
    spin11="-↿⇂-"
    spin10="-↿--"
    spin01="--⇂-"
    spin00="----"
    return if !iszero(α)
        if !iszero(β)
            spin11
        else
            spin10
        end
    else
        if !iszero(β)
            spin01
        else
            spin00
        end
    end
end

function _MOstring(orbit_ridx)
    if orbit_ridx isa Integer
        return if orbit_ridx == 0
            "HOMO"
        elseif orbit_ridx < 0
            "HOMO" * string(orbit_ridx)
        elseif orbit_ridx == 1
            "LUMO"
        else
            "LUMO+" * string(orbit_ridx-1)
        end
    else # open-shell (spin, idx)
        return (orbit_ridx[1]==1 ? "α-" : "β-") * _MOstring(orbit_ridx[2])
    end
end

function Serialize(t::GenericMolecule)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    atoms       = t.atoms
    atom_coords = t.atom_coords
    charge      = t.charge
    spin        = t.spin
    name        = t.name
    rot_alp     = t.rot_α
    rot_bet     = t.rot_β
    rot_gam     = t.rot_γ
    @pack! dict = (type, atoms, atom_coords, charge, spin, name, rot_alp, rot_bet, rot_gam)
    return dict
end