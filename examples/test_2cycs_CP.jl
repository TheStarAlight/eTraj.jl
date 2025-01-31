using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos4Laser(peak_int=0.4PW/cm^2, wave_len=800.0nm, cyc_num=2, ellip=1.0)
t = get_atom("H")

for init_cond in [:ADK, :SPANE, :SPA]
    perform_traj_simulation(
        init_cond_method    = init_cond,
        traj_phase_method   = :CTMC,
        laser               = l,
        target              = t,
        dimension           = 2,            # two-dimensional simulation, xy plane only
        sample_t_intv       = (-100,100),   # equivalent to `(-2.42fs, 2.42fs)`
        sample_t_num        = 20000,        # will sample 20000 equidistant time points between -100 and 100 a.u.
        traj_t_final        = 120,          # the traj end at 120 a.u., equivalent to `2.90fs`
        final_p_max         = (2.5,2.5),    # the momentum spec collection grid's border (-2.5 to +2.5 a.u.)
        final_p_num         = (500,500),    # the momentum spec collection grid's size (500x500)
        ss_kd_max           = 2.0,
        ss_kd_num           = 10000,        # will sample 10000 equidistant k⟂ points between -2 to +2 a.u.
        output_path         = "$(init_cond)-CTMC_4e14_800nm_cos4_2cyc_CP.jld2"
    )
end