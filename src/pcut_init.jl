function pcut_init()
    # Initialize all of the *_sav arrays to help prevent bleeding over between pcuts or ion species
    l_save .= false  # Whole array must be initialized in case number of particles changes from pcut to pcut
    weight_sav      .= 0.0
    ptot_pf_sav       .= 0.0
    pb_pf_sav       .= 0.0
    x_PT_cm_sav     .= 0.0
    grid_sav        .= 0
    DwS_sav         .= false
    inj_sav         .= false
    xn_per_sav      .= 0.0
    prp_x_cm_sav    .= 0.0
    acctime_sec_sav .= 0.0
    φ_rad_sav       .= 0.0
    tcut_sav        .= 0

    # A separate variable tracks the number of finished particles, so
    # that race conditions can be avoided in OMP mode
    i_fin = 0

    # For high-energy electrons in a strong magnetic field, need to know
    # previous cutoff momentum for calculating new PRP downstream
    pcut_prev = i_cut > 1 ? pcuts_use[i_cut-1] : 0.0
end
