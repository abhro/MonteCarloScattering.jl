function iter_finalize()
    # Compute the escaping flux for this iteration
    pₓ_esc_flux_UpS[i_iter] = pₓ_esc_UpS / flux_px_UpS
    energy_esc_flux_UpS[i_iter] = energy_esc_UpS / flux_energy_UpS

    # Compute the adiabatic index everywhere on grid now that all particles are accounted
    # for. First store value from previous iteration, since both the pre- and post-iteration
    # adiabatic indices will be used in the profile smoothing subroutine smooth_grid
    if i_iter == 1
        index_mask = (x_grid_cm .≤ 0)
        Γ_grid[  index_mask,1] .= 5//3         #    #assumecold
        Γ_grid[.!index_mask,1] .= Γ₂_RH
    else
        Γ_grid[:,1] .= Γ_grid[:,2]
    end

    Γ_grid[:,2] .= @. 1 + (pressure_psd_par + pressure_psd_perp) / energy_density_psd

    index_mask = (energy_density_psd .== 1e-99)
    Γ_grid[index_mask,2] .= 1e-99

    # Also compute the adiabatic index of particles that were lost DwS
    Γ_DwS[i_iter] = 1 + ∑P_DwS/∑KEdensity_DwS

    # Calculate expected escaping fluxes, now that adiabatic index is known
    # far DwS. Also average them so that the smoothing subroutine treats
    # calculated and actual escaping fluxes identically
    q_esc_cal_pₓ[i_iter], q_esc_cal_energy[i_iter] = q_esc_calcs(Γ_DwS[i_iter],)
    n_avg = min(i_iter, 4)
    q_esc_cal_pₓ_avg = mean(q_esc_cal_pₓ[i_iter-n_avg+1:i_iter])
    q_esc_cal_energy_avg = mean(q_esc_cal_energy[i_iter-n_avg+1:i_iter])

    # In runs with OpenMP, race conditions may cause roundoff error at the
    # least significant decimal place of the **_flux variables. This causes
    # slight variations in the smoothing algorithm across different runs.
    # Correct for that here by rounding off the last two decimal places
    # (the second place is a fudge factor).
    pxx_flux[1:n_grid] .= round(pxx_flux[1:n_grid], digits=13)
    # Commented out because it's not necessary for smoothing parallel shocks
    #pxz_flux[1:n_grid] = round(pxz_flux[1:n_grid], digits=13)
    energy_flux[1:n_grid] = round(energy_flux[1:n_grid], digits=13)

    # Output grid data for this iteration, and smooth the grid for the next iteration
    smooth_grid_par(i_iter, i_shock, n_grid, x_grid_rg, x_grid_cm,
                    Γ_grid, uz_sk_grid, θ_grid,
                    pressure_psd_par, pressure_psd_perp, flux_px_UpS,
                    flux_energy_UpS, Γ₂_RH, q_esc_cal_pₓ_avg,
                    q_esc_cal_energy_avg, pxx_flux, energy_flux, uₓ_sk_grid,
                    γ_sf_grid, btot_grid, utot_grid, γ_ef_grid,
                    β_ef_grid, εB_grid)


    # Compute average escaping flux over last four iterations and write to
    # file; average over all iterations if i_iter < 4.
    n_avg          = min(i_iter, 4)
    pₓ_esc_avg     = mean(pₓ_esc_flux_UpS[1:n_avg])
    energy_esc_avg = mean(energy_esc_flux_UpS[1:n_avg])

    println(outfile, " Parallel shock q_esc from Double et al (2004) equations:")
    println(outfile, "     Esc. energy flux/UpS    = ", q_esc_cal_energy_avg)
    println(outfile, "     Esc. momentum flux/UpS  = ", q_esc_cal_pₓ_avg)
    energy_esc_flux_UpS[i_iter] = max(energy_esc_flux_UpS[i_iter], 1e-99)
    energy_esc_avg             = max(energy_esc_avg,             1e-99)
    pₓ_esc_flux_UpS[i_iter] = max(pₓ_esc_flux_UpS[i_iter], 1e-99)
    pₓ_esc_avg             = max(pₓ_esc_avg,             1e-99)
    println(outfile,
            " Esc. en flux FEB/UpS  for i_iter = ", i_iter, ":   en esc = ",
            energy_esc_flux_UpS[i_iter], "   Avg. esc en  = ", energy_esc_avg)
    println(outfile,
            " Esc. pxx flux FEB/UpS for i_iter = ", i_iter, ":  pxx esc = ",
            pₓ_esc_flux_UpS[i_iter], "   Avg. esc pxx = ", pₓ_esc_avg)

    if iszero(q_esc_cal_pₓ_avg)
        println(outfile, " Avg q_pₓ_MC/q_pₓ_cal N/A, because q_pₓ_cal = 0")
    else
        println(outfile, " Avg q_pₓ_MC/q_pₓ_cal. = ", pₓ_esc_avg/q_esc_cal_pₓ_avg)
    end
    if iszero(q_esc_cal_energy_avg)
        println(outfile, " Avg q_energy_MC/q_energy_cal N/A, because q_energy_cal = 0")
    else
        println(outfile, " Avg q_energy_MC/q_energy_cal. = ", energy_esc_avg/q_esc_cal_energy_avg)
    end
    println(outfile)


    # Compute various adiabatic indices and write them out
    n_avg     = min(i_iter, 4)
    Γ_DwS_esc = mean(Γ_DwS[i_iter-n_avg+1:i_iter])
    Γ_UpS     = 5//3 # #assumecold

    println(outfile, " Iteration #", i_iter)
    println(outfile, "   r_comp = ", r_comp, "      r_RH = ", r_RH)
    println(outfile, "   Adiab index for far UpS particles = ", Γ_UpS)
    println(outfile, "   Adiab index for DwS PRP particles = ", Γ_DwS_esc)
    println(outfile, "   Adiab index from R-H relations    = ", Γ₂_RH)
    println(outfile)


    # If tcut tracking was enabled, print out the results here
    if do_tcuts
        tcut_print(weights_file, spectra_file, n_tcuts, tcuts, n_ions,
                   num_psd_mom_bins, psd_mom_bounds,
                   i_iter, weight_coupled, spectra_coupled)
    end
end
