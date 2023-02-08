using StaticArrays

"Represents an atom under single-active-electron (SAE) approximation."
struct SAEAtom <: SAEAtomBase
# should implement TargetPotential, TargetForce, TrajectoryFunction, ADKRateExp.
    "Ionization potential of the atom."
    IonPotential;
    "Asymptotic charge of the inner nucleus."
    NuclCharge;
    "Atomic parameters used to fit the atomic potential. See [J. Phys. B 38 2593 (2005)]"
    a1;
    b1;
    a2;
    b2;
    a3;
    b3;
    "Name of the atom."
    name::String;
    "Initializes a new instance of SAEAtom."
    SAEAtom( Ip, Z, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]") = new(Ip,Z,a1,b1,a2,b2,a3,b3,name)
    SAEAtom(;Ip, Z, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]") = new(Ip,Z,a1,b1,a2,b2,a3,b3,name)
end

"Gets the ionization potential of the atom."
IonPotential(t::SAEAtom) = t.IonPotential
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::SAEAtom) = t.NuclCharge
"Gets the name of the atom."
TargetName(t::SAEAtom) = t.name
"""
Gets the potential function of the atom.
Expression: V(r) = - [Z + a1*exp(-b1*r) + a2*r*exp(-b2*r) + a3*exp(-b3*r)] / r.
"""
function TargetPotential(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return - (t.NuclCharge + t.a1*exp(-t.b1*r) + t.a2*r*exp(-t.b2*r) + t.a3*exp(-t.b3*r)) / r
    end
end
"""
Gets the force exerted on the electron from the atom (which is the neg-grad of potential).
Expression: F(rvec) = - rvec / r^3 * [ Z + a1*(1+b1*r)*exp(-b1*r) + a3*(1+b3*r)*exp(-b3*r) ] - rvec * a2*b2/r * exp(-b2*r).
"""
function TargetForce(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return (-x,-y,-z) .* (r^(-3) * (t.NuclCharge + t.a1*(1+t.b1*r)*exp(-t.b1*r) + t.a3*(1+t.b3*r)*exp(-t.b3*r)) + t.a2*t.b2/r * exp(-t.b2*r))
    end
end
"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::SAEAtom, laserFx::Function, laserFy::Function, phaseMethod::Symbol, nonDipole::Bool; kwargs...)
    Z  = t.NuclCharge
    Ip = t.IonPotential
    a1,b1,a2,b2,a3,b3 = t.a1,t.b1,t.a2,t.b2,t.a3,t.b3
    return if ! nonDipole
        if phaseMethod == :CTMC
            function traj_dipole_ctmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
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
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r)
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phaseMethod == :SCTS
            function traj_dipole_scts(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
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
ADKRateExp(t::SAEAtom) = (F,φ,pd,pz) -> exp(-2(pd^2+pz^2+2*t.IonPotential)^1.5/3F)

"Prints the information of the atom."
Base.show(io::IO, t::SAEAtom) = print(io, "[SAEAtom] Atom $(t.name), Ip=$(t.IonPotential), Z=$(t.NuclCharge)")