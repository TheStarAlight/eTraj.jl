
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
    "Orientation of the quantization axis θ."
    quan_ax_θ;
    "Orientation of the quantization axis ϕ."
    quan_ax_ϕ;
    "Soft core parameter of the atom."
    soft_core;
    "Name of the atom."
    name::String;
    "Initializes a new instance of `HydrogenLikeAtom`."
    function HydrogenLikeAtom(Ip, Z::Integer, l::Integer=0, m::Integer=0, quan_ax_orient_θ::Real=0.0, quan_ax_orient_ϕ::Real=0.0, soft_core::Real=0.2, name="[NA]")
        @assert Ip>0 "[HydrogenLikeAtom] Ip should be positive."
        @assert l≥0 && m≥0 && l≥abs(m) "[HydrogenLikeAtom] Invalid (l,m)."
        @assert soft_core≥0 "[HydrogenLikeAtom] Soft core should be non-negative."
        new(Ip, Z, l, m, quan_ax_orient_θ, quan_ax_orient_ϕ, soft_core, name)
    end
    HydrogenLikeAtom(;Ip, Z::Integer, l::Integer=0, m::Integer=0, quan_ax_orient_θ::Real=0.0, quan_ax_orient_ϕ::Real=0.0, soft_core::Real=0.2, name="[NA]") = HydrogenLikeAtom(Ip, Z, l, m, quan_ax_orient_θ, quan_ax_orient_ϕ, soft_core, name)
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

using SpecialFunctions
"Gets the asymptotic coefficient C_κl of the atom using the Hartree approximation formula."
function AsympCoeff(t::HydrogenLikeAtom)
    n = t.nucl_charge/sqrt(2*t.Ip)
    return 2^(n-1) / sqrt(n*gamma(n+t.l+1)*gamma(n-t.l))
end
"Gets the name of the atom."
TargetName(t::HydrogenLikeAtom) = t.name
"Gets the potential function of the atom."
TargetPotential(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-0.5)
"Gets the force exerted on the electron from the atom (which is the neg-grad of potential)."
TargetForce(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-1.5) .* (x,y,z)

using StaticArrays
"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::HydrogenLikeAtom, laserFx::Function, laserFy::Function, phase_method::Symbol; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    soft_core = t.soft_core
    if phase_method == :CTMC
        function traj_dipole_ctmc(u,p,t)
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
        function traj_dipole_qtmc(u,p,t)
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
        function traj_dipole_scts(u,p,t)
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

using Printf
"Prints the information of the atom."
function Base.show(io::IO, t::HydrogenLikeAtom)
    @printf(io, "[HydrogenLikeAtom] Atom %s, Ip=%.4f, Z=%i, soft_core=%.4f\n", t.name, t.Ip, t.nucl_charge, t.soft_core)
end

using Parameters, OrderedCollections
"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(t::HydrogenLikeAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    l           = t.l
    m           = t.m
    soft_core   = t.soft_core
    name        = t.name
    @pack! dict = (type, Ip, nucl_charge, l, m, soft_core, name)
    return dict
end