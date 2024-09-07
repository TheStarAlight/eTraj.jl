
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
    "Asymptotic coefficient C_κl."
    asymp_coeff;
    "Orientation of the quantization axis θ."
    quan_ax_θ;
    "Orientation of the quantization axis ϕ."
    quan_ax_ϕ;
    "Parameters used to fit the atomic potential. See [J. Phys. B 38 2593 (2005)]"
    a1;
    b1;
    a2;
    b2;
    a3;
    b3;
    "Name of the atom."
    name::String;
end

"""
    SAEAtom(Ip, Z [,l=0] [,m=0] [,asymp_coeff=:hartree|<coeff>] [,quan_ax_θ=0.0] [,quan_ax_ϕ=0.0] [,a1,b1,a2,b2,a3,b3] [,name]) <: SAEAtomBase

Initializes a new instance of `SAEAtom`.

## Parameters
- `Ip`  : Ionization potential of the atom (numerically in **a.u.** or a `Unitful.Quantity`).
- `Z`   : Nuclear charge number.
- `l=0` : Angular quantum number (*optional, default 0*).
- `m=0` : Magnetic quantum number (*optional, default 0*).
- `asymp_coeff=:hartree`: Asymptotic coefficient related to the wavefunction's behavior when r→∞ (*`:hartree` or a positive number*). Passing `:hartree` (by default) indicates automatic calculation using the Hartree formula.
- `quan_ax_θ=0.0`       : Orientation angle θ of the quantization axis relative to the lab frame (*optional, default 0.0*).
- `quan_ax_ϕ=0.0`       : Orientation angle ϕ of the quantization axis relative to the lab frame (*optional, default 0.0*).
- `a1,b1,a2,b2,a3,b3`   : Parameters used to fit the atomic potential. See [*J. Phys. B* **38**, 2593 (2005)]
- `name::String`        : Name of the atom.
"""
function SAEAtom(;Ip, Z::Integer, l::Integer=0, m::Integer=0, asymp_coeff=:hartree, quan_ax_θ::Real=0.0, quan_ax_ϕ::Real=0.0, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., name="[NA]")
    (Ip isa Quantity) && (Ip = (uconvert(eV,Ip) |> auconvert).val)
    @assert Ip>0 "[SAEAtom] Ip should be positive."
    @assert l≥0 && m≥0 && l≥abs(m) "[SAEAtom] Invalid (l,m)."
    @assert b1≥0 && b2≥0 && b3≥0 "[SAEAtom] b1,b2,b3 should be non-negative."
    @assert asymp_coeff in [:hartree] || asymp_coeff > 0 "[SAEAtom] asymp_coeff should be either `:hartree` or a positive number."
    C = 0.0
    if asymp_coeff == :hartree
        C = hartree_asymp_coeff(Z,Ip,l)
    else
        C = asymp_coeff
    end
    SAEAtom(Ip,Z,l,m,C,quan_ax_θ,quan_ax_ϕ,a1,b1,a2,b2,a3,b3,name)
end

"Gets the ionization potential of the atom."
IonPotential(t::SAEAtom) = t.Ip
"Gets the asymptotic nuclear charge of the atom."
AsympNuclCharge(t::SAEAtom) = t.nucl_charge
"Gets the angular quantum number l of the atom."
AngularQuantumNumber(t::SAEAtom) = t.l
"Gets the magnetic quantum number m of the atom."
MagneticQuantumNumber(t::SAEAtom) = t.m

"Gets the orientation of the quantization axis of the atom in spherical coordinates (θ,ϕ)."
QuantizationAxisOrientaion(t::SAEAtom) = (t.quan_ax_θ, t.quan_ax_ϕ)

"Gets the asymptotic coefficient C_κl of the atom."
AsympCoeff(t::SAEAtom) = t.asymp_coeff

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

"Gets the trajectory function according to given parameter."
function TrajectoryFunction(t::SAEAtom, dimension::Integer, laserFx::Function, laserFy::Function, phase_method::Symbol; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    a1,b1,a2,b2,a3,b3 = t.a1,t.b1,t.a2,t.b2,t.a3,t.b3
    if dimension == 2
        if phase_method == :CTMC
            function traj_dipole_ctmc_2d(u,p,t)
                # tFx, tFy = targetF(u[1],u[2])
                r = sqrt(u[1]^2+u[2]^2)
                tFx, tFy = (u[1],u[2]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                @SVector [du1,du2,du3,du4]
            end
        elseif phase_method == :QTMC
            function traj_dipole_qtmc_2d(u,p,t)
                r = sqrt(u[1]^2+u[2]^2)
                tFx, tFy = (u[1],u[2]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]))
                du5 = -(Ip + (du1^2+du2^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r)
                @SVector [du1,du2,du3,du4,du5]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_2d(u,p,t)
                r = sqrt(u[1]^2+u[2]^2)
                tFx, tFy = (u[1],u[2]) .* -(r^(-3)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/r*exp(-b2*r))
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]) + (u[1]*tFx+u[2]*tFy))
                du5 = -(Ip + (du1^2+du2^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/r + (u[1]*tFx+u[2]*tFy))
                @SVector [du1,du2,du3,du4,du5]
            end
        end
    else
        if phase_method == :CTMC
            function traj_dipole_ctmc_3d(u,p,t)
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
            function traj_dipole_qtmc_3d(u,p,t)
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
            function traj_dipole_scts_3d(u,p,t)
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
    end
end

"Prints the information of the atom."
function Base.show(io::IO, t::SAEAtom)
    @printf(io, "[SAEAtom] Atom %s, Ip=%.4f, Z=%i", t.name, t.Ip, t.nucl_charge)
end

"Returns a `Dict{Symbol,Any}` containing properties of the object."
function Serialize(t::SAEAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    l           = t.l
    m           = t.m
    asymp_coeff = t.asymp_coeff
    name        = t.name
    a1 = t.a1
    b1 = t.b1
    a2 = t.a2
    b2 = t.b2
    a3 = t.a3
    b3 = t.b3
    @pack! dict = (type, Ip, nucl_charge, l, m, asymp_coeff, a1,b1,a2,b2,a3,b3, name)
    return dict
end