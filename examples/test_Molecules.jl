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
