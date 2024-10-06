using eTraj
using eTraj.Targets
using eTraj.Lasers
using eTraj.Units

l = Cos4Laser(peak_int=0.4PW/cm^2, wave_len=800.0nm, cyc_num=2, ellip=1.0)
t = get_atom("H")

for init_cond in [:ADK, :SPANE, :SPA]
    perform_traj_simulation(
            init_cond_method    = init_cond,
            laser               = l,
            target              = t,
            dimension           = 2,
            sample_t_intv       = (-100,100),
            sample_t_num        = 20000,
            traj_t_final        = 120,
            final_p_max         = (2.5,2.5),
            final_p_num         = (500,500),
            ss_kd_max           = 2.0,
            ss_kd_num           = 10000,
            output_path         = "$(init_cond)-CTMC_4e14_800nm_cos4_2cyc_CP.jld2",
            traj_phase_method   = :CTMC,
            traj_rtol           = 1e-6,
            rate_prefix         = :Full
        )
end