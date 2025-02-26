# Example: WFAT-CTMC for Molecules

The WFAT framework offers a rigorous approach for computing tunneling ionization probabilities in molecular systems, particularly advantageous for complex molecular structures.
This section demonstrates the application of the WFAT-CTMC simulation scheme to various molecular targets.

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
        traj_phase_method   = :CTMC,    # WFAT supports CTMC only
        laser               = l,
        target              = t[i],
        mol_orbit_ridx      = orbit_ridx[i],
        dimension           = 2,
        sample_t_intv       = (-300,300),
        sample_t_num        = 10000,
        traj_t_final        = 350,
        final_p_max         = (2.0,2.0),
        final_p_num         = (500,500),
        ss_kd_max           = 2.0,
        ss_kd_num           = 5000,
        output_path         = path[i]
    )
end
```

The simulation employs a 6-cycle, circularly polarized 800-nm laser pulse. This duration ensures PMD symmetry for isotropic atoms while enabling effective orbital imaging through the PMD structure: nodal planes in the molecular orbital manifest as nodes or dark regions in the PMD when intersected by the laser's electric field vector.
The molecular orientation is configured to align the molecular frame (MF) ``z``-axis parallel to the laboratory frame (LF) ``x``-axis, achieved through a 90° counterclockwise rotation around the ``y``-axis (specified by `rot_β=90°`).

The PMDs obtained from the HOMOs of different molecules are presented in the figure below.
The structures of the PMDs reflect the geometries of the corresponding orbitals, which are further mirrored in the orientation-dependent structure factors, as illustrated in another figure below.

Fig. (a) displays the PMD of the hydrogen molecule's HOMO orbital, ``1\rm{s}\sigma_{\rm{g}}``, which primarily arises from the in-phase combination of the two 1s orbitals of the hydrogen atoms.
The overall spherical shape of the orbital is consistent with the lowest-order [``\nu=(n_\xi,m)=(0,0)``] squared structure factor ``\abs{G_{00}}^2`` shown in sub figure (a), resulting in the evenly distributed ring-like structure of the PMD.

The carbon monoxide (CO) molecule, a heteronuclear diatomic molecule, has a HOMO orbital designated as ``3\sigma_{\rm{g}}``.
The ``3\sigma_{\rm{g}}`` orbital of the CO molecule differs from a conventional ``\sigma_{\rm{g}}`` orbital of a homonuclear molecule, as the electron density is predominantly localized on the carbon atom. This localization leads to a significant increase in the squared structure factor and ionization probability when the negative electric field is oriented toward the carbon atom.
This localized peak within the ring structure is evident in the PMD of the CO molecule, as shown in sub figure (b).

The HOMOs of the oxygen (``\rm{O_2}``) and benzene (``\rm{C_6H_6}``) molecules, depicted in sub figures (c) and (d), respectively, are degenerate orbitals with ``\pi`` symmetry.
The oxygen molecule's "α-HOMO" (one of the 2p``\pi_{\rm{u}}`` orbitals) and the benzene molecule's "HOMO" (``2\rm{p}\pi_3``), after a rotation of `rot_β=90°`, exhibit nodal planes along the ``x-z`` and ``y-z`` directions, respectively. Consequently, the rotating electric field in the ``x-y`` plane interacts with a four-lobe structure, as revealed in the PMDs shown in sub figures (c1) and (d1).
In contrast, the oxygen molecule's "α-HOMO-1" (another ``2\rm{p}\pi_{\rm{u}}`` orbital) and the benzene molecule's "HOMO-1" (``2\rm{p}\pi_2``) each possess a nodal plane in the ``x-y`` plane. This implies that, in the zeroth order (``m=0``), the outgoing electron waves originating from the "``+``" and "``-``" regions of the orbital interfere destructively, resulting in a net-zero ionization probability.
For non-zero-``m`` channels, ``\mathcal{W}_\nu(F,\kt=0)\equiv 0``, leading to the appearance of a nodal ring corresponding to zero initial momenta (``\kt=0``) in the PMDs of sub figures (c2) and (d2).
Under these conditions, ionization is suppressed, and the first-order channels (``m=\pm 1``) dominate the contribution to the ionization probability.

![fig:example_Molecules](assets/figure_Molecules.png)
![fig:example_Molecules_G](assets/figure_Molecules_G_orbital.png)
