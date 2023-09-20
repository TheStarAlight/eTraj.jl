# some shared methods

"Generates a 2D point inside a given area x∈[-x0,x0], y∈[-y0,y0] using the random generator `rng`."
function gen_rand_pt(rng, x0,y0)
    x,y = (rand(rng)-0.5)*2*x0, (rand(rng)-0.5)*2*y0
    return x,y
end

using Printf
using StaticArrays
using LinearAlgebra
using SphericalHarmonics
"Computes the Y_lm(̂k(ts))."
@inline function sph_harm_lm_khat(l,m, kxts,kyts,kzts, x_unit,y_unit,z_unit) # for atoms the func is invoked once per sample, use the more precise version
    cosθ = (kxts*z_unit[1]+kyts*z_unit[2]+kzts*z_unit[3])/sqrt(kxts^2+kyts^2+kzts^2) |> abs
    (cosθ≤1.0) && (cosθ=1.0)
    return if m == 0
        SphericalHarmonics.sphericalharmonic(acos(complex(abs(cosθ))),0.0; l=l,m=m)
    else
        # determine ϕ: F is set as the new z axis, and the new x and y axis are determined using the rotation angles
        ϕ = atan(y_unit ⋅ real([kxts,kyts,kzts]), x_unit ⋅ real([kxts,kyts,kzts]))
        SphericalHarmonics.sphericalharmonic(acos(complex(abs(cosθ))),ϕ; l=l,m=m)
    end
end
"Computes the Y_lm(̂k(ts)) using the leading-order expansion around cosθ=1."
@inline function sph_harm_lm_khat_approx(l,m, kxts,kyts,kzts, x_unit,y_unit,z_unit) # for molecules the func would be invoked hundreds of times per sample, which greatly affects the efficiency, use the approximated version
    cosθ = (kxts*z_unit[1]+kyts*z_unit[2]+kzts*z_unit[3])/sqrt(kxts^2+kyts^2+kzts^2) |> abs # cosθ is real, but the result can be near +1 or -1 because sqrt is bi-valued, we can give it an abs to make it near +1.
    if cosθ≤1.0 # for SFA, some solutions are not accurate enough, which results in |cosθ| ≲ 1.0, simply set to 1.0
        # @printf "kxts=%.4f+%.6fim; kyts=%.4f+%.6fim; kzts=%.4f+%.6fim, k²(ts)=%.4f+%.6fim, cosθ=%.4f+%.4fim\n" real(kxts) imag(kxts) real(kyts) imag(kyts) real(kzts) imag(kzts) real(kxts^2+kyts^2+kzts^2) imag(kxts^2+kyts^2+kzts^2) real(cosθ) imag(cosθ)
        cosθ = 1.0
    end
    Q_lm = (-1)^m * sqrt((2l+1)/2*gamma(l+abs(m)+1)/gamma(l-abs(m)+1)) # use gamma function instead of factorial to avoid overflow
    return if m == 0
        Q_lm / sqrt(2π)
    else
        # determine ϕ: F is set as the new z axis, and the new x and y axis are determined using the rotation angles
        ϕ = atan(y_unit ⋅ real([kxts,kyts,kzts]), x_unit ⋅ real([kxts,kyts,kzts]))
        Q_lm / gamma(abs(m)+1) * (abs(abs(cosθ)-1)/2)^(abs(m)/2) / sqrt(2π) * cis(m*ϕ)
    end
end

using Rotations
"""
Obtains the rotational Euler angles `(α,β,γ)` from FF to the MF.
The FF is defined as the frame where the field F is aligned with the z axis.
"""
function obtain_FF_MF_Euler(mol_rot, F)
    mα, mβ, mγ = mol_rot
    # note: in the definition of RotZYZ(t1,t2,t3), the rot consists of Z(t3)->Y(t2)->Z(t1), while is extrinsic, but in our definition the rot is intrinsic, thus we reversed the order of angles
    RotMFLF = inv(RotZYZ(mα, mβ, mγ)) # MF → LF, inversed because it is initially an active rotation (rot of vector) but a passive rot is needed
    RotLFFF = RotZYZ(atan(F[2],F[1]), π/2, 0) # LF → FF
    α,β,γ = Rotations.params(RotZYZ(inv(RotLFFF*RotMFLF)))
    return α,β,γ
end
"""
Obtains the rotated axes `(x_axis,y_axis,z_axis)` in FF represented in LF.
"""
function obtain_xyz_FF_LF(Fx, Fy)
    LFFFRotMat = RotZYZ(atan(Fy,Fx), π/2, 0.0) # active rot from Lf to FF
    new_x_axis = LFFFRotMat * @SVector [1.0,0.0,0.0]
    new_y_axis = LFFFRotMat * @SVector [0.0,1.0,0.0]
    new_z_axis = LFFFRotMat * @SVector [0.0,0.0,1.0]
    return new_x_axis, new_y_axis, new_z_axis
end