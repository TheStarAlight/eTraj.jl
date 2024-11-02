using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos2Laser(peak_int=90.0TW/cm^2, wave_len=800.0nm, cyc_num=1, cep=π/2, ellip=0.0)
t = get_atom("H"; soft_core=1e-12)

for phase_method in [:QTMC, :SCTS]
    perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 2,
        sample_t_intv       = (-50,50),
        sample_t_num        = 30000,
        traj_t_final        = 100,
        final_p_max         = (1.0,1.0),
        final_p_num         = (500,500),
        ss_kd_max           = 1.5,
        ss_kd_num           = 10000,
        output_path         = "ADK-$(phase_method)_9e13_800nm_1cyc_LP_ExpRate.jld2",
        traj_phase_method   = phase_method,
        rate_prefix         = :Exp
    )
end