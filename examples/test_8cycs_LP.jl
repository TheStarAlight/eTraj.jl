using eTraj
using eTraj.Targets
using eTraj.Lasers
using eTraj.Units

l = Cos2Laser(peak_int=90.0TW/cm^2, wave_len=800.0nm, cyc_num=8, ellip=0.0)
t = get_atom("H"; soft_core=1e-10)

perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 3,
        sample_t_intv       = (-350,350),
        sample_t_num        = 50000,
        traj_t_final        = 500,
        final_p_max         = (1.0,0.0,1.0),
        final_p_num         = (500,1,500),
        ss_kd_max           = 0,
        ss_kd_num           = 1,
        ss_kz_max           = 1.0,
        ss_kz_num           = 20000,
        output_path         = "ADK-QTMC_9e13_800nm_8cyc_LP_ExpRate.jld2",
        traj_phase_method   = :QTMC,
        traj_rtol           = 1e-6,
        rate_prefix         = :Exp
    )
# ✓ Launching electrons and collecting ... [batch #50000/50000, 311663176 electrons collected]    Time: 0:35:37
# Progress: 100%[●●●●●●●●●●●●●●●●●●●●●●●●●] Time: 0:35:37 (42.75 ms/it)