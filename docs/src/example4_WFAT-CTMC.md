# Example: WFAT-CTMC for Molecules

The WFAT provides a precise means of calculating the probability of tunneling ionization of molecules, especially for complex molecular targets.
In this section we present an example of using the WFAT-CTMC simulation scheme for molecular targets.

```julia
# examples/test_Molecules.jl

using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos2Laser(peak_int=4e14W/cm^2, wave_len=800.0nm, cyc_num=6, ellip=1.0)
t = [get_mol("Hydrogen"; rot_β=90°),
     get_mol("Carbon Monoxide"; rot_β=90°),
     get_mol("Oxygen"; rot_β=90°),
     get_mol("Oxygen"; rot_β=90°),
     get_mol("Benzene"; rot_β=90°),
     get_mol("Benzene"; rot_β=90°)]
orbit_ridx = [0, 0, (1,0), (1,-1), 0, -1]
path = [
    "WFAT-CTMC_Hydrogen_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_CarbonMonoxide_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Oxygen_α-HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Oxygen_α-HOMO-1_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Benzene_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Benzene_HOMO-1_4e14_800nm_6cyc_CP.jld2"
]
for i in eachindex(t)
    perform_traj_simulation(
        init_cond_method    = :WFAT,
        laser               = l,
        target              = t[i],
        dimension           = 2,
        sample_t_intv       = (-300,300),
        sample_t_num        = 10000,
        traj_t_final        = 350,
        final_p_max         = (2.0,2.0),
        final_p_num         = (500,500),
        ss_kd_max           = 2.0,
        ss_kd_num           = 5000,
        output_path         = path[i],
        traj_phase_method   = :CTMC,    # WFAT supports CTMC only
        mol_orbit_ridx      = orbit_ridx[i]
    )
end
```

In this code, we set the molecules' orientation such that the ``z``-axis in the MF is parallel to the ``x``-axis in the LF, which can be realized by setting `rot_β=90°` (equivalent to a 90° counterclockwise rotation of the molecule around the ``y``-axis).
The laser is a circularly polarized 800-nm laser pulse, whose duration is set to 6 cycles to ensure the symmetry of the PMD for an isotropic atom, which facilitates the imaging of the molecules' orbitals through the PMD: a node or a dark curve would show up in the PMD if the laser's electric field vector crosses or scans inside a nodal plane of the orbital's wavefunction.

The PMDs obtained from different molecule's HOMO orbitals are shown below.
The structures of the PMDs reflect the geometries of the orbitals, which are also mirrored in the orientation-dependent structure factors, which are also shown in the figure.

PMD Fig. (a) shows the PMD of the hydrogen molecule's HOMO orbital 1s``\sigma_{\rm{g}}``, which is mainly contributed by the in-phase combination of two 1s orbitals from the two hydrogen atoms.
The overall shape of the orbital is similar to a sphere, which is consistent with the lowest order [``\nu=(n_\xi,m)=(0,0)``] squared structure factor ``\abs{G_{00}}^2`` shown in the structure factor Fig. (a), and finally leads to the evenly-distributed ring-like shape of the PMD.

The carbon monoxide (CO) molecule is a heteronuclear diatomic molecule, whose HOMO orbital is named after 3``\sigma_{\rm{g}}``.
The 3``\sigma_{\rm{g}}`` orbital of CO molecule is different from a conventional ``\sigma_{\rm{g}}`` orbital of a homonuclear molecule in that the electron density is concentrated more on the carbon atom [see the structure factor Fig. (b)], which results in a dramatic increase in the squared structure factor and ionization probability when the negative electric field points towards the carbon atom.
Such a localized peak in the ring structure is observed in the PMD of CO molecule in PMD Fig. (b).

The HOMOs of the oxygen (``\rm{O_2}``) and benzene (``\rm{C_6H_6}``) molecules (in MF) are shown in the structure factor Figs. (c) and (d), which are degenerate orbitals of ``\pi`` symmetry.
The oxygen's "α-HOMO" (one of 2p``\pi_{\rm{u}}``) and the benzene's "HOMO" (2p``\pi_3``), after the given rotation (`rot_β=90°`), have ``x-z`` and ``y-z`` as their nodal planes, hence the rotating electric field in the ``x-y`` plane "sees" a four-lobe structure, which is revealed in the PMDs in Figs. (c1) and (d1).
Whereas the oxygen's "α-HOMO-1" (another one of 2p``\pi_{\rm{u}}``) and the benzene's "HOMO-1" (2p``\pi_2``), each have a nodal plane on the ``x-y`` plane, which indicates that in the zeroth order (``m=0``), the outgoing electron waves that are contributed by the "``+``" parts and the "``-``" parts of the orbital cancel out each other, leading to a net-zero ionization probability.
For non-zero-``m`` channels we have ``\mathcal{W}_\nu(F,\kt=0)\equiv 0`` [see Eq. (57) in [WFAT](@ref WFAT)], hence a nodal ring which corresponds to zero initial momenta ``\kt=0`` would appear in the PMDs of Figs. (c2) and (d2).
Under such circumstance, the ionization gets suppressed, and the first-order channels (``m=\pm 1``) would take dominance in the contribution of ionization probability.

![fig:example_Molecules](assets/figure_Molecules.png)
![fig:example_Molecules_G](assets/figure_Molecules_G_orbital.png)
