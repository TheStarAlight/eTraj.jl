
"Represents a Hydrogen-like atom."
struct HydrogenLikeAtom <: SAEAtomBase
    "Ionization potential of the atom."
    Ip;
    "Asymptotic charge of the inner nucleus."
    nucl_charge;
    "Angular quantum number l."
    l;
    "Magnetic quantum number m."
    m;
    "Asymptotic coefficient C_κl."
    asymp_coeff;
    "Orientation of the quantization axis θ."
    quan_ax_θ;
    "Orientation of the quantization axis ϕ."
    quan_ax_ϕ;
    "Soft core parameter of the atom."
    soft_core;
    "Name of the atom."
    name::String;
end

"""
    HydrogenLikeAtom(Ip, Z [,l=0] [,m=0] [,asymp_coeff=:hartree|<coeff>] [,quan_ax_θ=0.0] [,quan_ax_ϕ=0.0] [,soft_core=0.2] [,name]) <: SAEAtomBase

Initializes a new instance of `HydrogenLikeAtom`.

## Parameters
- `Ip`  : Ionization potential of the atom (numerically in **a.u.** or a `Unitful.Quantity`).
- `Z`   : Nuclear charge number.
- `l=0` : Angular quantum number (*optional, default 0*).
- `m=0` : Magnetic quantum number (*optional, default 0*).
- `asymp_coeff=:hartree` : Asymptotic coefficient related to the wavefunction's behavior when r→∞ (*`:hartree` or a positive number*). Passing `:hartree` (by default) indicates automatic calculation using the Hartree formula.
- `quan_ax_θ=0.0`   : Orientation angle θ of the quantization axis relative to the lab frame (*optional, default 0.0*).
- `quan_ax_ϕ=0.0`   : Orientation angle ϕ of the quantization axis relative to the lab frame (*optional, default 0.0*).
- `soft_core=0.2`   : Soft core parameter of the Coulomb potential (*optional, default 0.2*).
- `name::String`    : Name of the atom.
"""
function HydrogenLikeAtom(;Ip, Z::Integer, l::Integer=0, m::Integer=0, asymp_coeff=:hartree, quan_ax_θ::Real=0.0, quan_ax_ϕ::Real=0.0, soft_core::Real=0.2, name="[NA]")
    (Ip isa Quantity) && (Ip = (uconvert(eV,Ip) |> auconvert).val)
    @assert Ip>0 "[HydrogenLikeAtom] `Ip` should be positive."
    @assert l≥0 && m≥0 && l≥abs(m) "[HydrogenLikeAtom] Invalid (l,m)."
    @assert soft_core≥0 "[HydrogenLikeAtom] `soft_core` should be non-negative."
    @assert asymp_coeff in [:hartree] || asymp_coeff > 0 "[HydrogenLikeAtom] asymp_coeff should be either `:hartree` or a positive number."
    C = 0.0
    if asymp_coeff == :hartree
        C = hartree_asymp_coeff(Z,Ip,l)
    else
        C = asymp_coeff
    end
    HydrogenLikeAtom(Ip, Z, l, m, C, quan_ax_θ, quan_ax_ϕ, soft_core, name)
end

"Gets the ionization potential of the atom."
IonPotential(t::HydrogenLikeAtom) = t.Ip
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::HydrogenLikeAtom) = t.nucl_charge
"Gets the angular quantum number l of the atom."
AngularQuantumNumber(t::HydrogenLikeAtom) = t.l
"Gets the magnetic quantum number m of the atom."
MagneticQuantumNumber(t::HydrogenLikeAtom) = t.m
"Gets the soft core parameter of the atom."
SoftCore(t::HydrogenLikeAtom) = t.soft_core
"Gets the orientation of the quantization axis of the atom in spherical coordinates (θ,ϕ)."
QuantizationAxisOrientaion(t::HydrogenLikeAtom) = (t.quan_ax_θ, t.quan_ax_ϕ)

"Gets the asymptotic coefficient C_κl of the atom."
AsympCoeff(t::HydrogenLikeAtom) = t.asymp_coeff
"Gets the name of the atom."
TargetName(t::HydrogenLikeAtom) = t.name
"Gets the potential function of the atom."
TargetPotential(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-0.5)
"Gets the force exerted on the electron from the atom (which is the neg-grad of potential)."
TargetForce(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-1.5) .* (x,y,z)

"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::HydrogenLikeAtom, dimension::Integer, laserFx::Function, laserFy::Function, phase_method::Symbol; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    soft_core = t.soft_core
    if dimension == 2
        if phase_method == :CTMC
            function traj_dipole_ctmc_2d(u,p,t)
                # tFx, tFy = targetF(u[1],u[2])
                tFx, tFy = -Z*(u[1]^2+u[2]^2+soft_core)^(-1.5) .* (u[1],u[2])
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                @SVector [du1,du2,du3,du4]
            end
        elseif phase_method == :QTMC
            function traj_dipole_qtmc_2d(u,p,t)
                tFx, tFy = -Z*(u[1]^2+u[2]^2+soft_core)^(-1.5) .* (u[1],u[2])
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]))
                du5 = -(Ip + (du1^2+du2^2)/2 - Z*(u[1]^2+u[2]^2+soft_core)^(-0.5))
                @SVector [du1,du2,du3,du4,du5]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_2d(u,p,t)
                tFx, tFy = -Z*(u[1]^2+u[2]^2+soft_core)^(-1.5) .* (u[1],u[2])
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]) + (u[1]*tFx+u[2]*tFy))
                du5 = -(Ip + (du1^2+du2^2)/2 - Z*(u[1]^2+u[2]^2+soft_core)^(-0.5) + (u[1]*tFx+u[2]*tFy))
                @SVector [du1,du2,du3,du4,du5]
            end
        end
    else
        if phase_method == :CTMC
            function traj_dipole_ctmc_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
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
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                tFx, tFy, tFz = -Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-1.5) .* (u[1],u[2],u[3])
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    end
end

"Prints the information of the atom."
function Base.show(io::IO, t::HydrogenLikeAtom)
    @printf(io, "[HydrogenLikeAtom] Atom %s, Ip=%.4f, Z=%i, soft_core=%.4f", t.name, t.Ip, t.nucl_charge, t.soft_core)
end

"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(t::HydrogenLikeAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    l           = t.l
    m           = t.m
    asymp_coeff = t.asymp_coeff
    soft_core   = t.soft_core
    name        = t.name
    @pack! dict = (type, Ip, nucl_charge, l, m, asymp_coeff, soft_core, name)
    return dict
end