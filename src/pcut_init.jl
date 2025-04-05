function pcut_init(i_cut)
    # Initialize all of the *_sav arrays to help prevent bleeding over between pcuts or ion species
    l_save .= false  # Whole array must be initialized in case number of particles changes from pcut to pcut
    for arr in (:weight_sav, :ptot_pf_sav, :pb_pf_sav, :x_PT_cm_sav, :grid_sav, :DwS_sav,
                :inj_sav, :xn_per_sav, :prp_x_cm_sav, :acctime_sec_sav, :φ_rad_sav, :tcut_sav)
        # zero out each array in a type stable manner
        @eval $arr .= zero(eltype($arr))
    end
    #weight_sav      .= 0.0
    #ptot_pf_sav     .= zero(eltype(ptot_pf_save))
    #pb_pf_sav       .= zero(eltype(pb_pf_save))
    #x_PT_cm_sav     .= zero(eltype(x_PT_cm_sav))
    #grid_sav        .= zero(eltype(grid_sav))
    #DwS_sav         .= false
    #inj_sav         .= false
    #xn_per_sav      .= 0.0
    #prp_x_cm_sav    .= 0.0
    #acctime_sec_sav .= 0.0
    #φ_rad_sav       .= 0.0
    #tcut_sav        .= 0

    # A separate variable tracks the number of finished particles, so
    # that race conditions can be avoided in OMP mode
    i_fin = 0

    # For high-energy electrons in a strong magnetic field, need to know
    # previous cutoff momentum for calculating new PRP downstream
    pcut_prev = i_cut > 1 ? pcuts_use[i_cut-1] : 0.0
end
