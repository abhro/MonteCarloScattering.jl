function pcut_finalize()
    break_pcut = false
    n_saved = count(l_save)

    t_end = now()
    run_time = t_end - t_start

    @info("", i_iter, i_ion, i_cut,
          pcuts_in[i_cut], pcuts_use[i_cut]/(mp_cgs*c_cgs),
          n_saved, n_pts_use, weight_running, run_time)
    println(outfile, " itr=", i_iter, " ion=", i_ion, " icut=", i_cut,
            pcuts_in[i_cut], pcuts_use[i_cut]/(mp_cgs*c_cgs),
            "  n_sav=", n_saved, "/", n_pts_use, weight_running, run_time)


    # If no particles saved, don't bother with remaining pcuts
    break_pcut = true
    if n_saved == 0
        return break_pcut
    end

    # Prepare population for next pcut
    n_pts_target = pcuts_use[i_cut] < p_pcut_hi ? n_pts_pcut : n_pts_pcut_hi

    global (
            grid_new, tcut_new, DwS_new, inj_new, weight_new, ptot_pf_new, pb_pf_new,
            x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new,
            n_pts_use, weight_running
           ) = new_pcut(
                        n_pts_target, n_saved, l_save, grid_sav, DwS_sav, inj_sav,
                        weight_sav, ptot_pf_sav, pb_pf_sav, x_PT_cm_sav, xn_per_sav,
                        prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav,
                          n_pts_use, weight_running)
    return break_pcut
end
