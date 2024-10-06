using eTraj
using eTraj.Targets
using eTraj.Lasers
using eTraj.Units

l = Cos2Laser(peak_int=90.0TW/cm^2, wave_len=800.0nm, cyc_num=1, cep=π/2, ellip=0.0)
t = get_atom("H"; soft_core=1e-12)

perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 3,
        sample_t_intv       = (-50,50),
        sample_t_num        = 30000,
        traj_t_final        = 100,
        final_p_max         = (1.0,0.0,1.0),
        final_p_num         = (500,1,500),
        ss_kd_max           = 0,
        ss_kd_num           = 1,
        ss_kz_max           = 1.5,
        ss_kz_num           = 10000,
        output_path         = "ADK-SCTS_9e13_800nm_1cyc_LP_ExpRate.jld2",
        traj_phase_method   = :SCTS,
        traj_rtol           = 1e-6,
        rate_prefix         = :Exp
    )
# ✓ Launching electrons and collecting ... [batch #30000/30000, 74901216 electrons collected]    Time: 0:06:17
# Progress: 100%[●●●●●●●●●●●●●●●●●●●●●●●●●] Time: 0:06:17 (12.58 ms/it)