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

"""
Obtains the rotational Euler angles `(α,β,γ)` from FF to the MF.
The FF is defined as the frame where the field F_vec is aligned with the z axis.
This transform is only valid under the case when F_vec is in xy plane in LF.
"""
function obtain_FF_MF_Euler(mol_rot, F_vec)
    mα, mβ, mγ = mol_rot
    Fx = F_vec[1]; Fy = F_vec[2]; F = hypot(Fx,Fy)
    FF_LFrep = [0 Fy/F Fx/F; 0 -Fx/F Fy/F; 1 0 0]
    MF_LFrep = RotZYZ(mα,mβ,mγ)
    α,β,γ = params(RotZYZ(MF_LFrep\FF_LFrep)) # basis transformation is a column transformation, the transition matrix is multiplied on the right side.
    return α,β,γ
end
"""
Obtains the rotated axes `(x_axis,y_axis,z_axis)` in FF represented in LF.
"""
function obtain_xyz_FF_LF(Fx, Fy)
    F = hypot(Fx,Fy)
    new_x_axis = @SVector [ 0.0 ,  0.0, 1.0]
    new_y_axis = @SVector [ Fy/F,-Fx/F, 0.0]
    new_z_axis = @SVector [ Fx/F, Fy/F, 0.0]
    return new_x_axis, new_y_axis, new_z_axis
end