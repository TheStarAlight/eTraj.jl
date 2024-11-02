using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos2Laser(peak_int=90.0TW/cm^2, wave_len=800.0nm, cyc_num=8, ellip=0.0)
t = get_atom("H")

for phase_method in [:QTMC, :SCTS]
    perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 2,
        sample_t_intv       = (-350,350),
        sample_t_num        = 50000,
        traj_t_final        = 500,
        final_p_max         = (1.0,1.0),
        final_p_num         = (500,500),
        ss_kd_max           = 1.0,
        ss_kd_num           = 20000,
        output_path         = "ADK-$(phase_method)_9e13_800nm_8cyc_LP_ExpRate.jld2",
        traj_phase_method   = phase_method,
        rate_prefix         = :Exp
    )
end