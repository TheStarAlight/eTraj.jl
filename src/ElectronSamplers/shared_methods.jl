# some shared methods

"Generates a 2D point inside a given area x∈[-x0,x0], y∈[-y0,y0] using the random generator `rng`."
function gen_rand_pt_2dsq(rng, x0,y0)
    x,y = (rand(rng)-0.5)*2*x0, (rand(rng)-0.5)*2*y0
    return x,y
end

"Generates a 1D point inside x∈[-x0,x0] using the random generator `rng`."
function gen_rand_pt_1d(rng, x0)
    return (rand(rng)-0.5)*2*x0
end

"Generates a set of spherical harmonics Y_lm(kx,ky,kz) that evaluates in the Cartesian coordinate where k=(kx,ky,kz) satisfies the saddle-point equation k²=-κ²."
@inline function gen_sph_harm_funcs(lmax, κ)
    fact = factorial    # shortcut of factorial
    # R(l,m) = fact(l+m)*mapreduce(k->(-1)^k/(2^(2k+m)*fact(k+m)*fact(k)*fact(l-m-2k)) * (x+1im*y)^(k+m) * (x-1im*y)^k * z^(l-m-2k), +, max(0,-m):round(Int, (l-m)/2, RoundDown)) # solid harmonics
    # Y(l,m) = (-1)^m * sqrt((2l+1)/4π*fact(l-m)/fact(l+m)) / r^l * R(l,m)
    SHfunc = Matrix(undef, lmax+1, 2lmax+1)
    for l in 0:lmax for m in -l:l
        SHfunc[l+1, l+m+1] = (x,y,z) -> (-1)^m * sqrt((2l+1)/4π*fact(l-m)/fact(l+m)) / (-1im*κ)^l * fact(l+m)*mapreduce(k->(-1)^k/(2^(2k+m)*fact(k+m)*fact(k)*fact(l-m-2k)) * (x+1im*y)^(k+m) * (x-1im*y)^k * z^(l-m-2k), +, max(0,-m):round(Int, (l-m)/2, RoundDown))
    end; end
    return SHfunc
end

using StaticArrays
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