
"""
    struct HydrogenLikeAtom <: SAEAtomBase <: Target

Represents a Hydrogen-like atom.
"""
struct HydrogenLikeAtom <: SAEAtomBase
    Ip;
    nucl_charge;
    l;
    m;
    asymp_coeff;
    quan_ax_θ;
    quan_ax_ϕ;
    soft_core;
    name::String;
end

"""
    HydrogenLikeAtom(Ip, Z [,l=0] [,m=0] [,asymp_coeff=:hartree|<coeff>] [,quan_ax_θ=0.0] [,quan_ax_ϕ=0.0] [,soft_core=1e-10] [,name]) <: SAEAtomBase

Initializes a new `HydrogenLikeAtom`.

## Parameters
- `Ip`  : Ionization potential of the atom (numerically in **a.u.** or a `Unitful.Quantity`).
- `Z`   : Nuclear charge number.
- `l=0` : Angular quantum number (*optional, default 0*).
- `m=0` : Magnetic quantum number (*optional, default 0*).
- `asymp_coeff=:hartree` : Asymptotic coefficient related to the wavefunction's behavior when r→∞ (*`:hartree` or a positive number*). Passing `:hartree` (by default) indicates automatic calculation using the Hartree formula.
- `quan_ax_θ=0.0`   : Orientation angle θ of the quantization axis relative to the lab frame (numerically in radian or a `Unitful.Quantity`) (*optional, default 0.0*).
- `quan_ax_ϕ=0.0`   : Orientation angle ϕ of the quantization axis relative to the lab frame (numerically in radian or a `Unitful.Quantity`) (*optional, default 0.0*).
- `soft_core=1e-10` : Soft core parameter of the Coulomb potential (*optional, default 1e-10*).
- `name::String`    : Name of the atom.

## Examples
```jldoctest
julia> t = HydrogenLikeAtom(Ip=0.5, Z=1, name="H")
[HydrogenLikeAtom] Atom H, Ip=0.5000 (13.61 eV), Z=1

julia> using eTraj.Units

julia> t = HydrogenLikeAtom(Ip=3.4eV, Z=1, l=1, name="H")
[HydrogenLikeAtom] Atom H (p orbital, m=0), Ip=0.1249 (3.40 eV), Z=1
```

## See Also
The [`get_atom`](@ref) method provides some atom presets for use.
"""
function HydrogenLikeAtom(;Ip, Z::Integer, l::Integer=0, m::Integer=0, asymp_coeff=:hartree, quan_ax_θ=0.0, quan_ax_ϕ=0.0, soft_core=1e-10, name="[NA]")
    (Ip isa Quantity) && (Ip = (uconvert(eV,Ip) |> auconvert).val)
    (quan_ax_θ isa Quantity) && (quan_ax_θ=uconvert(u"rad",quan_ax_θ).val)
    (quan_ax_ϕ isa Quantity) && (quan_ax_ϕ=uconvert(u"rad",quan_ax_ϕ).val)
    @assert Ip>0 "[HydrogenLikeAtom] `Ip` should be positive."
    @assert l≥0 && l≥abs(m) "[HydrogenLikeAtom] Invalid (l,m)."
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

IonPotential(t::HydrogenLikeAtom) = t.Ip
AsympNuclCharge(t::HydrogenLikeAtom) = t.nucl_charge
AngularQuantumNumber(t::HydrogenLikeAtom) = t.l
MagneticQuantumNumber(t::HydrogenLikeAtom) = t.m
SoftCore(t::HydrogenLikeAtom) = t.soft_core
QuantizationAxisOrientaion(t::HydrogenLikeAtom) = (t.quan_ax_θ, t.quan_ax_ϕ)

AsympCoeff(t::HydrogenLikeAtom) = t.asymp_coeff
TargetName(t::HydrogenLikeAtom) = t.name
TargetPotential(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-0.5)
TargetForce(t::HydrogenLikeAtom) = (x,y,z) -> -t.nucl_charge*(x^2+y^2+z^2+t.soft_core)^(-1.5) .* (x,y,z)

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
                # du5 = -((du1^2+du2^2)/2 + targetP(u[1],u[2]) + (u[1]*tFx+u[2]*tFy))
                du5 = -((du1^2+du2^2)/2 - Z*(u[1]^2+u[2]^2+soft_core)^(-0.5) + (u[1]*tFx+u[2]*tFy))
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
                # du7 = -((du1^2+du2^2+du3^2)/2 + targetP(u[1],u[2],u[3]) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                du7 = -((du1^2+du2^2+du3^2)/2 - Z*(u[1]^2+u[2]^2+u[3]^2+soft_core)^(-0.5) + (u[1]*tFx+u[2]*tFy+u[3]*tFz))
                @SVector [du1,du2,du3,du4,du5,du6,du7]
            end
        end
    end
end

function Base.show(io::IO, t::HydrogenLikeAtom)
    @printf(io, "[HydrogenLikeAtom] Atom %s%s, Ip=%.4f (%.2f eV), Z=%i", t.name, (t.l==0 ? "" : " ($(l_info(t.l)), m=$(t.m))"), t.Ip, auconvert(eV, t.Ip).val, t.nucl_charge)
    (t.quan_ax_θ!=0 || t.quan_ax_ϕ!=0) && @printf(io, ", θϕ=(%.1f°,%.1f°)", t.quan_ax_θ*180/π, t.quan_ax_ϕ*180/π)
end

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
