using StaticArrays

"Represents a Hydrogen-like atom."
struct HydrogenLikeAtom <: SAEAtomBase
    "Ionization potential of the atom."
    IonPotential;
    "Asymptotic charge of the inner nucleus."
    NuclCharge;
    "Soft core parameter of the atom."
    soft_core;
    "Name of the atom."
    name::String;
    "Initializes a new instance of HydrogenLikeAtom."
    HydrogenLikeAtom( Ip, Z, soft_core=1.0, name="[NA]") = new(Ip, Z, soft_core, name)
    HydrogenLikeAtom(;Ip, Z, soft_core=1.0, name="[NA]") = new(Ip, Z, soft_core, name)
end

"Gets the ionization potential of the atom."
IonPotential(t::HydrogenLikeAtom) = t.IonPotential
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::HydrogenLikeAtom) = t.NuclCharge
"Gets the soft core parameter of the atom."
SoftCore(t::HydrogenLikeAtom) = t.soft_core
"Gets the name of the atom."
TargetName(t::HydrogenLikeAtom) = t.name
"Gets the potential function of the atom."
TargetPotential(t::HydrogenLikeAtom) = (x,y,z) -> -t.NuclCharge*(x^2+y^2+z^2+t.soft_core)^(-0.5)
"Gets the force exerted on the electron from the atom (which is the neg-grad of potential)."
TargetForce(t::HydrogenLikeAtom) = (x,y,z) -> -t.NuclCharge*(x^2+y^2+z^2+t.soft_core)^(-1.5) .* (x,y,z)
"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::HydrogenLikeAtom, laserFx::Function, laserFy::Function, phaseMethod::Symbol, nonDipole::Bool; kwargs...)
    Z  = t.NuclCharge
    Ip = t.IonPotential
    soft_core = t.soft_core
    # including external function call is infeasible in GPU, thus the external targetF & targetP are replaced by pure Coulomb ones.
    return if ! nonDipole
        if phaseMethod == :CTMC
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
        elseif phaseMethod == :QTMC
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
        elseif phaseMethod == :SCTS
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
    else
        #TODO: add support for nondipole simulation.
    end
end
"""
Gets the exponential term of ADK rate which depends on
Field strength `F`,
Azimuthal angle of field `φ`,
momentum's transverse component `pd` (in xy plane),
and propagation-direction (which is Z axis) component `pz`.
"""
ADKRateExp(t::HydrogenLikeAtom) = (F,φ,pd,pz) -> exp(-2(pd^2+pz^2+2*t.IonPotential)^1.5/3F)

"Prints the information of the atom."
Base.show(io::IO, t::HydrogenLikeAtom) = print(io, "[HydrogenLikeAtom] Atom $(t.name), Ip=$(t.IonPotential), Z=$(t.NuclCharge), SoftCore=$(t.soft_core).")