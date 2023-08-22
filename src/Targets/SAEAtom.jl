
"Represents an atom under single-active-electron (SAE) approximation."
struct SAEAtom <: SAEAtomBase
# should implement TargetPotential, TargetForce, TrajectoryFunction.
    "Ionization potential of the atom."
    Ip;
    "Asymptotic charge of the inner nucleus."
    nucl_charge;
    "Angular quantum number l."
    l;
    "Magnetic quantum number m."
    m;
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
    function SAEAtom( Ip, Z::Integer, l::Integer=0, m::Integer=0, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]")
        @assert Ip>0 "[SAEAtom] Ip should be positive."
        @assert l≥0 && m≥0 && l≥abs(m) "[SAEAtom] Invalid (l,m)."
        @assert b1≥0 && b2≥0 && b3≥0 "[SAEAtom] b1,b2,b3 should be non-negative."
        new(Ip,Z,l,m,a1,b1,a2,b2,a3,b3,name)
    end
    SAEAtom(;Ip, Z::Integer, l::Integer=0, m::Integer=0, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]") = SAEAtom(Ip,Z,l,m,a1,b1,a2,b2,a3,b3,name)
end

"Gets the ionization potential of the atom."
IonPotential(t::SAEAtom) = t.Ip
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::SAEAtom) = t.nucl_charge
"Gets the angular quantum number l of the atom."
AngularQuantumNumber(t::SAEAtom) = t.l
"Gets the magnetic quantum number m of the atom."
MagneticQuantumNumber(t::SAEAtom) = t.m

using SpecialFunctions
"Gets the asymptotic coefficient C_κl of the atom using the Hartree approximation formula."
function AsympCoeff(t::SAEAtom)
    n = t.nucl_charge/sqrt(2*t.Ip)
    return 2^(n-1) / sqrt(n*gamma(n+t.l+1)*gamma(n-t.l))
end
"Gets the name of the atom."
TargetName(t::SAEAtom) = t.name
"""
Gets the potential function of the atom.
Expression: V(r) = - [Z + a1*exp(-b1*r) + a2*r*exp(-b2*r) + a3*exp(-b3*r)] / r.
"""
function TargetPotential(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return - (t.nucl_charge + t.a1*exp(-t.b1*r) + t.a2*r*exp(-t.b2*r) + t.a3*exp(-t.b3*r)) / r
    end
end
"""
Gets the force exerted on the electron from the atom (which is the neg-grad of potential).
Expression: F(rvec) = - rvec / r^3 * [ Z + a1*(1+b1*r)*exp(-b1*r) + a3*(1+b3*r)*exp(-b3*r) ] - rvec * a2*b2/r * exp(-b2*r).
"""
function TargetForce(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return (-x,-y,-z) .* (r^(-3) * (t.nucl_charge + t.a1*(1+t.b1*r)*exp(-t.b1*r) + t.a3*(1+t.b3*r)*exp(-t.b3*r)) + t.a2*t.b2/r * exp(-t.b2*r))
    end
end

using StaticArrays
"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::SAEAtom, laserFx::Function, laserFy::Function, phase_method::Symbol, non_dipole::Bool; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    a1,b1,a2,b2,a3,b3 = t.a1,t.b1,t.a2,t.b2,t.a3,t.b3
    return if ! non_dipole
        if phase_method == :CTMC
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
        elseif phase_method == :QTMC
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
        elseif phase_method == :SCTS
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
        if phase_method == :CTMC
            function traj_nondipole_ctmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                @SVector [du1,du2,du3,du4,du5,du6]
            end
        elseif phase_method == :QTMC
            function traj_nondipole_qtmc(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r)
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_nondipole_scts(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                c0 = 137.035999173
                du1 = u[4] + u[3]*laserFx(t)/c0
                du2 = u[5] + u[3]*laserFy(t)/c0
                du3 = u[6]
                du4 = tFx - laserFx(t)
                du5 = tFy - laserFy(t)
                du6 = tFz - (u[4]*laserFx(t)+u[5]*laserFy(t))/c0
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    end
end

using Printf
"Prints the information of the atom."
function Base.show(io::IO, t::SAEAtom)
    @printf(io, "[SAEAtom] Atom %s, Ip=%.4f, Z=%i\n", t.name, t.Ip, t.nucl_charge)
end

using Parameters, OrderedCollections
"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(t::SAEAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    l           = t.l
    m           = t.m
    name        = t.name
    a1 = t.a1
    b1 = t.b1
    a2 = t.a2
    b2 = t.b2
    a3 = t.a3
    b3 = t.b3
    @pack! dict = (type, Ip, nucl_charge, l, m, a1,b1,a2,b2,a3,b3, name)
    return dict
end