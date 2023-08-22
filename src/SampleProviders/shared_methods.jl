# some shared methods

"Generates a 2D point inside a circle of radius `r0` using the random generator `rng`."
function gen_rand_pt_circ(rng, r0)
    x,y = (rand(rng)-0.5)*2r0, (rand(rng)-0.5)*2r0
    while x^2+y^2>r0^2
        x,y = (rand(rng)-0.5)*2r0, (rand(rng)-0.5)*2r0
    end
    return x,y
end

using StaticArrays
using LinearAlgebra
using SphericalHarmonics
"Computes the Y_lm(̂k(ts))."
function sph_harm_lm_khat(l,m,(kxts,kyts,kzts),(Fxtr,Fytr))
    cos_theta = (kxts*Fxtr+kyts*Fytr)/sqrt(kxts^2+kyts^2+kzts^2)/sqrt(Fxtr^2+Fytr^2)
    return if m == 0
        SphericalHarmonics.sphericalharmonic(acos(complex(cos_theta)),0.0; l=l,m=m)
    else
        # determine ϕ: F is set as the new z axis, Z is set as the new x axis.
        y_unit = normalize((@SVector [0.0,0.0,1.0]) × (@SVector [Fxtr,Fytr,0.0]))
        ϕ = atan(real(y_unit ⋅ @SVector [kxts,kyts,kzts]), real(kzts))
        SphericalHarmonics.sphericalharmonic(acos(complex(cos_theta)),ϕ; l=l,m=m)
    end
end

using Rotations
"""
Obtains the rotational Euler angles `(α,β,γ)` from LF to the MF.
The LF is defined as the frame where the field F is aligned with the z axis.
"""
function obtain_Euler(mol_rot, F)
    F_unit = normalize(F)
    mα, mβ, mγ = mol_rot
    RotMol = RotZYZ(mγ, mβ, mα)
    RotLaser = RotMatrix3([ 0 -F_unit[2] F_unit[1];
                            0  F_unit[1] F_unit[2];
                           -1  0         0         ])
    α,β,γ = Rotations.params(RotZYZ(inv(RotLaser)*RotMol))
    return α,β,γ
end