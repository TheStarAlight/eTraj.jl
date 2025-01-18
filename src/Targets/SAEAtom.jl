
"""
    struct SAEAtom <: SAEAtomBase <: Targets

Represents an atom under single-active-electron (SAE) approximation.
"""
struct SAEAtom <: SAEAtomBase
    Ip;
    nucl_charge;
    l;
    m;
    asymp_coeff;
    quan_ax_θ;
    quan_ax_ϕ;
    a1;
    b1;
    a2;
    b2;
    a3;
    b3;
    soft_core;
    name::String;
end

"""
    SAEAtom(Ip, Z [,l=0] [,m=0] [,asymp_coeff=:hartree|<coeff>] [,quan_ax_θ=0.0] [,quan_ax_ϕ=0.0] [,a1,b1,a2,b2,a3,b3] [,soft_core=1e-10] [,name]) <: SAEAtomBase

Initializes a new `SAEAtom`.

## Parameters
- `Ip`  : Ionization potential of the atom (numerically in **a.u.** or a `Unitful.Quantity`).
- `Z`   : Nuclear charge number.
- `l=0` : Angular quantum number (*optional, default 0*).
- `m=0` : Magnetic quantum number (*optional, default 0*).
- `asymp_coeff=:hartree`: Asymptotic coefficient related to the wavefunction's behavior when r→∞ (*`:hartree` or a positive number*). Passing `:hartree` (by default) indicates automatic calculation using the Hartree formula.
- `quan_ax_θ=0.0`       : Orientation angle θ of the quantization axis relative to the lab frame (numerically in radian or a `Unitful.Quantity`) (*optional, default 0.0*).
- `quan_ax_ϕ=0.0`       : Orientation angle ϕ of the quantization axis relative to the lab frame (numerically in radian or a `Unitful.Quantity`) (*optional, default 0.0*).
- `a1,b1,a2,b2,a3,b3`   : Parameters used to fit the atomic potential. See [*J. Phys. B* **38**, 2593 (2005)]
- `soft_core=1e-10`     : Soft core parameter to avoid singularity in the potential (*optional, default 1e-10*).
- `name::String`        : Name of the atom.

## Examples
```jldoctest
julia> t = SAEAtom(Ip=0.9035, Z=1, asymp_coeff=:hartree, a1=1.230723, b1=0.6620055, a2=-1.325040, b2=1.236224, a3=-0.2307230, b3=0.4804286, name="He")
[SAEAtom] Atom He, Ip=0.9035 (24.59 eV), Z=1

julia> using eTraj.Units

julia> t = SAEAtom(Ip=12.13eV, Z=1, l=1, a1=51.35554, b1=2.111554, a2=-99.92747, b2=3.737221, a3=1.644457, b3=0.4306465, asymp_coeff=1.3, name="Xe")
[SAEAtom] Atom Xe (p orbital, m=0), Ip=0.4458 (12.13 eV), Z=1
```

## See Also
The [`get_atom`](@ref) method provides some atom presets for use.
"""
function SAEAtom(;Ip, Z::Integer, l::Integer=0, m::Integer=0, asymp_coeff=:hartree, quan_ax_θ=0.0, quan_ax_ϕ=0.0, a1=0., b1=0., a2=0., b2=0., a3=0., b3=0., soft_core=1e-10, name="[NA]")
    (Ip isa Quantity) && (Ip = (uconvert(eV,Ip) |> auconvert).val)
    (quan_ax_θ isa Quantity) && (quan_ax_θ=uconvert(u"rad",quan_ax_θ).val)
    (quan_ax_ϕ isa Quantity) && (quan_ax_ϕ=uconvert(u"rad",quan_ax_ϕ).val)
    @assert Ip>0 "[SAEAtom] Ip should be positive."
    @assert l≥0 && l≥abs(m) "[SAEAtom] Invalid (l,m)."
    @assert b1≥0 && b2≥0 && b3≥0 "[SAEAtom] b1,b2,b3 should be non-negative."
    @assert asymp_coeff in [:hartree] || asymp_coeff > 0 "[SAEAtom] asymp_coeff should be either `:hartree` or a positive number."
    @assert soft_core>0 "[SAEAtom] soft_core should be positive."
    C = 0.0
    if asymp_coeff == :hartree
        C = hartree_asymp_coeff(Z,Ip,l)
    else
        C = asymp_coeff
    end
    SAEAtom(Ip,Z,l,m,C,quan_ax_θ,quan_ax_ϕ,a1,b1,a2,b2,a3,b3,soft_core,name)
end


IonPotential(t::SAEAtom) = t.Ip
AsympNuclCharge(t::SAEAtom) = t.nucl_charge
AngularQuantumNumber(t::SAEAtom) = t.l
MagneticQuantumNumber(t::SAEAtom) = t.m
QuantizationAxisOrientaion(t::SAEAtom) = (t.quan_ax_θ, t.quan_ax_ϕ)
AsympCoeff(t::SAEAtom) = t.asymp_coeff
SoftCore(t::SAEAtom) = t.soft_core
TargetName(t::SAEAtom) = t.name

function TargetPotential(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return - (t.nucl_charge + t.a1*exp(-t.b1*r) + t.a2*r*exp(-t.b2*r) + t.a3*exp(-t.b3*r)) / sqrt(r^2+t.soft_core)
    end
end

function TargetForce(t::SAEAtom)
    return function(x,y,z)
        r = sqrt(x^2+y^2+z^2)
        return (-x,-y,-z) .* ((r^2+t.soft_core)^(-1.5) * (t.nucl_charge + t.a1*(1+t.b1*r)*exp(-t.b1*r) + t.a3*(1+t.b3*r)*exp(-t.b3*r)) + t.a2*t.b2/sqrt(r^2+t.soft_core) * exp(-t.b2*r))
    end
end

function TrajectoryFunction(t::SAEAtom, dimension::Integer, laserFx::Function, laserFy::Function, phase_method::Symbol; kwargs...)
    Z  = t.nucl_charge
    Ip = t.Ip
    sc = t.soft_core
    a1,b1,a2,b2,a3,b3 = t.a1,t.b1,t.a2,t.b2,t.a3,t.b3
    if dimension == 2
        if phase_method == :CTMC
            function traj_dipole_ctmc_2d(u,p,t)
                # tFx, tFy = targetF(u[1],u[2])
                r = sqrt(u[1]^2+u[2]^2)
                tFx, tFy = (u[1],u[2]) .* -((r^2+sc)^(-1.5)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/sqrt(r^2+sc)*exp(-b2*r))
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                @SVector [du1,du2,du3,du4]
            end
        elseif phase_method == :QTMC
            function traj_dipole_qtmc_2d(u,p,t)
                r = sqrt(u[1]^2+u[2]^2)
                tFx, tFy = (u[1],u[2]) .* -((r^2+sc)^(-1.5)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/sqrt(r^2+sc)*exp(-b2*r))
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -(Ip + (du1^2+du2^2)/2 + targetP(u[1],u[2]))
                du5 = -(Ip + (du1^2+du2^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/sqrt(r^2+sc))
                @SVector [du1,du2,du3,du4,du5]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_2d(u,p,t)
                r = sqrt(u[1]^2+u[2]^2)
                tFx, tFy = (u[1],u[2]) .* -((r^2+sc)^(-1.5)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/sqrt(r^2+sc)*exp(-b2*r))
                du1 = u[3]
                du2 = u[4]
                du3 = tFx-laserFx(t)
                du4 = tFy-laserFy(t)
                # du5 = -((du1^2+du2^2)/2 + targetP(u[1],u[2]) + (u[1]*tFx+u[2]*tFy))
                du5 = -((du1^2+du2^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/(r^2+sc) + (u[1]*tFx+u[2]*tFy))
                @SVector [du1,du2,du3,du4,du5]
            end
        end
    else
        if phase_method == :CTMC
            function traj_dipole_ctmc_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -((r^2+sc)^(-1.5)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/sqrt(r^2+sc)*exp(-b2*r))
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
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -((r^2+sc)^(-1.5)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/sqrt(r^2+sc)*exp(-b2*r))
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -(Ip + (du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]))
                du7 = -(Ip + (du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/sqrt(r^2+sc))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        elseif phase_method == :SCTS
            function traj_dipole_scts_3d(u,p,t)
                # tFx, tFy, tFz = targetF(u[1],u[2],u[3])
                r = sqrt(u[1]^2+u[2]^2+u[3]^2)
                tFx, tFy, tFz = (u[1],u[2],u[3]) .* -((r^2+sc)^(-1.5)*(Z+a1*(1+b1*r)*exp(-b1*r)+a3*(1+b3*r)*exp(-b3*r)) + a2*b2/sqrt(r^2+sc)*exp(-b2*r))
                du1 = u[4]
                du2 = u[5]
                du3 = u[6]
                du4 = tFx-laserFx(t)
                du5 = tFy-laserFy(t)
                du6 = tFz
                # du7 = -((du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -((du1^2+du2^2+du3^2)/2 - (Z+a1*exp(-b1*r)+a2*r*exp(-b2*r)+a3*exp(-b3*r))/sqrt(r^2+sc) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    end
end

function Base.show(io::IO, t::SAEAtom)
    @printf(io, "[SAEAtom] Atom %s%s, Ip=%.4f (%.2f eV), Z=%i", t.name, (t.l==0 ? "" : " ($(l_info(t.l)), m=$(t.m))"), t.Ip, auconvert(eV, t.Ip).val, t.nucl_charge)
    (t.quan_ax_θ!=0 || t.quan_ax_ϕ!=0) && @printf(io, ", θϕ=(%.1f°,%.1f°)", t.quan_ax_θ*180/π, t.quan_ax_ϕ*180/π)
end

function Serialize(t::SAEAtom)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(t)
    Ip          = t.Ip
    nucl_charge = t.nucl_charge
    l           = t.l
    m           = t.m
    asymp_coeff = t.asymp_coeff
    name        = t.name
    soft_core   = t.soft_core
    quan_ax_θ   = t.quan_ax_θ
    quan_ax_ϕ   = t.quan_ax_ϕ
    a1 = t.a1
    b1 = t.b1
    a2 = t.a2
    b2 = t.b2
    a3 = t.a3
    b3 = t.b3
    @pack! dict = (type, Ip, nucl_charge, l, m, asymp_coeff, quan_ax_θ, quan_ax_ϕ, a1,b1,a2,b2,a3,b3, soft_core, name)
    return dict
end