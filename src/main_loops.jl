"""
--------------------------------------------------------------------------
  Main computational loops
  Loops are set in the following sequence of nesting
   1. loop_itr:   i_iter   Iteration number
   2. loop_ion:   i_ion    Particle species number
   3. loop_pcut:  i_cut    Particle splitting (pcut) number
   4. loop_pt:    i_prt    Individual particle number
"""
function main_loops(
        n_itrs, n_ions, n_pcuts, n_grid, n_pts_inj, n_tcuts, species,
        (u₀, β₀, γ₀), (u₂, β₂, γ₂),
        Emax, Emax_per_aa, energy_pcut_hi, pmax,
        pxx_flux, pxz_flux, energy_flux,
        pressure_psd_par, pressure_psd_perp, energy_density_psd,
        esc_spectra_feb_upstream, esc_spectra_feb_downstream,
        weight_coupled,
        ε_target, Γ₂_RH, εB_grid, Γ_grid, γ_sf_grid, uₓ_sk_grid, uz_sk_grid, utot_grid, energy_transfer_frac,
        energy_transfer_pool, energy_recv_pool, energy_density, therm_energy_density,
        num_crossings, therm_grid, therm_pₓ_sk, therm_ptot_sk,
        therm_weight, psd, esc_psd_feb_upstream, esc_psd_feb_downstream,
        do_fast_push, inp_distr, energy_inj, inj_weight, ptot_inj, weight_inj,
        x_grid_start, x_grid_stop, rg₀, η_mfp, x_fast_stop_rg, x_grid_rg, x_grid_cm, n_pts_MB,
        xn_per_fine, xn_per_coarse, xn_per_sav, xn_per_new,
        weight_sav, weight_new, ptot_pf_sav, ptot_pf_new, pb_pf_sav, pb_pf_new,
        x_PT_cm_sav, x_PT_cm_new, grid_sav, grid_new, inj_sav, inj_new, downstream_sav, downstream_new,
        prp_x_cm_sav, prp_x_cm_new, acctime_sec_sav, acctime_sec_new, tcut_sav, tcut_new, φ_rad_sav, φ_rad_new,
        pcuts_in, pcuts_use, l_save, i_grid_feb, i_shock,
        n_pts_max, n_xspec, n_print_pt, num_psd_mom_bins, num_psd_θ_bins,
        psd_lin_cos_bins, psd_cos_fine, Δcos, psd_mom_bounds, psd_θ_bounds, psd_θ_min,
        psd_mom_min, psd_bins_per_dec_mom, psd_bins_per_dec_θ,
        bmag₂, zone_vol, pₑ_crit, γₑ_crit,
        x_spec, feb_upstream, feb_downstream, B_CMBz, use_custom_εB,
        γ_ef_grid, β_ef_grid, btot_grid, θ_grid, pₓ_esc_feb, energy_esc_feb,
        do_rad_losses, do_retro, do_tcuts, dont_DSA, dont_scatter, use_custom_frg,
        inj_fracs, spectra_pf, spectra_sf, tcuts, age_max, spectra_coupled,
        esc_energy_eff, esc_num_eff, esc_flux, electron_weight_fac,
        n_pts_pcut, n_pts_pcut_hi, t_start, weights_file, spectra_file, outfile,
        pₓ_esc_flux_upstream, flux_px_upstream, flux_energy_upstream, energy_esc_flux_upstream,
        Γ_downstream, q_esc_cal_pₓ, q_esc_cal_energy,
        do_multi_dNdps, do_photons,
        jet_rad_pc, jet_sph_frac, m_ion, aa_ion, zz_ion, T₀_ion, n₀_ion,
        r_comp, r_RH, n_shell_endpoints,
    )
    # Start of loop over iterations
    for i_iter in 1:n_itrs # loop_itr

        @info("Starting loop", i_iter)

        # Zero out numerous quantities that will be modified over the course of this iteration.
        # Minimally positive number is used to prevent errors when taking logarithms later
        # XXX uses the type-pirated version of Base.fill!
        fill!((pxx_flux, pxz_flux), 1e-99erg/cm^3)
        fill!((energy_flux, pressure_psd_par, pressure_psd_perp, energy_density_psd), 1e-99)
        fill!((esc_spectra_feb_upstream, esc_spectra_feb_downstream), 1e-99)
        fill!(weight_coupled, 1e-99)


        # Additionally, set/reset scalar quantities that will change
        ∑P_downstream  = 1e-99erg/cm^3         # total downstream pressure
        ∑KEdensity_downstream = 1e-99erg/cm^3  # total downstream kinetic energy density

        energy_esc_upstream    = 1e-99
        pₓ_esc_upstream        = 1e-99


        # To facilitate energy transfer from ions to electrons, calculate here the target energy
        # density fraction for electrons at each grid zone, and zero out the pool of
        # plasma-frame energy that will be taken from ions and donated to electrons
        # Per Ardaneh+ (2015) [10.1088/0004-637X/811/1/57], ε_electron is proportional to √ε_b.
        # ε_b is itself roughly proportional to density² -- B² is proportional to z² (z being the
        # density compression factor) -- so ε_electron should vary roughly linearly with density.
        z_max = γ₀ * β₀ / (γ₂ * β₂)
        populate_ε_target!(ε_target, z_max, γ_sf_grid, uₓ_sk_grid, u₀, γ₀, energy_transfer_frac)

        fill!((energy_transfer_pool, energy_recv_pool), 0.0erg)
        fill!((energy_density, therm_energy_density), 0.0)

        #------------------------------------------------------------------------
        #  Start of loop over particle species
        #
        #  First species always protons. If electrons present they MUST be last species.
        #  Each species has mass "aa" in units of proton mass and charge "zz" in units of
        #  elementary charge.
        #------------------------------------------------------------------------
        for i_ion in 1:n_ions # loop_ion

            zz = charge(species[i_ion])
            m = mass(species[i_ion])
            aa = m/mp |> NoUnits
            mc = m*c

            # At the start of each ion, print a glyph to the screen
            @info("Starting species iteration", i_iter, i_ion)

            pmax_cutoff = get_pmax_cutoff(Emax, Emax_per_aa, pmax, aa)

            # Zero out the phase space distributions and set variables related to
            # tracking thermal particles
            n_cr_count = clear_psd!(num_crossings, therm_grid, therm_pₓ_sk, therm_ptot_sk,
                                    therm_weight, psd, esc_psd_feb_upstream, esc_psd_feb_downstream)

            # In addition to initializing the phase space distribution, open the scratch (i.e.
            # temporary) file to which we will write information about thermal particle grid crossings
            nc_unit = open("mc_crossings.dat", "a")

            # To maintain identical results between OpenMP and serial versions,
            # set RNG seed based on current iteration/ion/pcut/particle number
            iseed_mod = (i_iter - 1)*n_ions + (i_ion - 1)
            Random.seed!(iseed_mod)


            # Initialize the particle populations that will be propagated through
            # the shock structure
            (n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in, x_PT_cm_in,
             pxx_flux, pxz_flux, energy_flux) = init_pop(
                do_fast_push, inp_distr, i_ion, m,
                temperature.(species), energy_inj, inj_weight, n_pts_inj,
                density.(species), x_grid_start, rg₀, η_mfp, x_fast_stop_rg,
                β₀, γ₀, u₀, n_ions, mass.(species),
                n_grid, x_grid_rg, uₓ_sk_grid, γ_sf_grid,
                ptot_inj, weight_inj, n_pts_MB,
            )
            @info("Finished init_pop on",
                  i_iter, i_ion,
                  n_pts_use, i_grid_in, weight_in, ptot_pf_in,
                  pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux)

            # Assign the various particle properties to the population
            assign_particle_properties_to_population!(
                n_pts_use, xn_per_fine, x_grid_stop,
                weight_new, weight_in, ptot_pf_new, ptot_pf_in,
                pb_pf_new, pb_pf_in, x_PT_cm_new, x_PT_cm_in, grid_new, i_grid_in,
                downstream_new, inj_new, xn_per_new, prp_x_cm_new, acctime_sec_new, tcut_new, φ_rad_new)

            # Weight of remaining particles, printed after each pcut; note that this will not be
            # correct for all particles if they were originally created so each thermal bin
            # would have equal weight
            weight_running = weight_in[1]

            # When using OpenMP, the array energy_transfer_pool can't be conditionally assigned
            # shared or reduction status, so it can't be used for both the ion and electron loops.
            # To get around this, use one array to hold the donated energy, and another to hold
            # the received energy.
            energy_recv_pool .= energy_transfer_pool


            # The array of pcuts read in by data_input has units momentum/mc.
            # Convert to momentum for this species
            pcuts_use[1:n_pcuts] .= pcuts_in[1:n_pcuts] * aa*mp*c

            # Determine the pcut at which to switch from low-E particle counts to high-E
            # particle counts. Recall that energy_pcut_hi has units of keV per aa, so when
            # dividing by particle mass the factor of aa is already present in the denominator.
            # Also set the maximum momentum cutoff based on the values given in keyword
            # "maximum-energy"
            p_pcut_hi = pcut_hi(energy_pcut_hi, E_rel_pt, mass(species[i_ion]))

            #----------------------------------------------------------------------
            #  Start of loop over pcuts
            #
            #  Particle splitting and resetting of n_pts_use is handled by the call
            #  to "new_pcut" at the end of each iteration of loop_pcut.
            #----------------------------------------------------------------------
            for i_cut in 1:n_pcuts # loop_pcut

                # Initialize all of the *_sav arrays to help prevent bleeding over between pcuts or ion species
                l_save .= false  # Whole array must be initialized in case number of particles changes from pcut to pcut

                zero!(weight_sav)
                zero!(ptot_pf_sav)
                zero!(pb_pf_sav)
                zero!(x_PT_cm_sav)
                zero!(grid_sav)
                zero!(downstream_sav)
                zero!(inj_sav)
                zero!(xn_per_sav)
                zero!(prp_x_cm_sav)
                zero!(acctime_sec_sav)
                zero!(φ_rad_sav)
                zero!(tcut_sav)

                # A separate variable tracks the number of finished particles,
                # so that race conditions can be avoided in OMP mode
                i_fin = 0

                # For high-energy electrons in a strong magnetic field, need to know
                # previous cutoff momentum for calculating new PRP downstream
                pcut_prev = i_cut > 1 ? pcuts_use[i_cut-1] : 0.0

                #--------------------------------------------------------------------
                #  Start of loop over particles
                #
                #  Quick explanation of particle properties. All units cgs where a unit exists.
                #
                #  weight: weighting value (used in momentum splitting)
                #  ptot_pf: total plasma frame momentum
                #  pb_pf: momentum along B field in plasma frame. NOT along x-axis unless
                #         upstream orientation is parallel
                #  r_PT_cm: current particle position
                #  i_grid: current grid zone number
                #  l_downstream: whether particle has been downstream
                #  inj: whether particle has been back upstream (i.e. is a CR)
                #  xn_per: time steps per gyro period, Δt = T_g/xn_per
                #  prp_x_cm: downstream position of PRP; adjusted to allow all particles,
                #            regardless of momentum, to isotropize before reaching
                #  acctime_sec: time since crossing shock for first time; not started until l_downstream = true
                #  φ_rad: phase angle of gyration relative to z axis
                #  tcut_curr: next time at which particle tracking takes place
                #--------------------------------------------------------------------
                #$omp parallel for default(none), schedule(dynamic,1), num_threads(6)
                for i_prt in 1:n_pts_use # loop_pt

                    (i_fin, i_reason, pb_pf, p_perp_b_pf, γₚ_pf, γᵤ_sf, b_cosθ, b_sinθ, weight,
                     uₓ_sk, uz_sk, utot, φ_rad,
                     ∑P_downstream, ∑KEdensity_downstream) = particle_loop(
                        i_iter, i_ion, i_cut, i_prt, i_grid_feb, i_shock,
                        n_ions, n_pcuts, n_pts_max, n_xspec, n_grid, n_print_pt, num_psd_mom_bins, num_psd_θ_bins,
                        psd_cos_fine, Δcos, psd_θ_min,
                        species,
                        γ₀, β₀, u₀, u₂, bmag₂,
                        pₑ_crit, γₑ_crit, η_mfp,
                        psd_mom_min, psd_bins_per_dec_mom, psd_bins_per_dec_θ,
                        energy_transfer_frac, energy_recv_pool,
                        psd, num_crossings, x_spec, feb_upstream, feb_downstream,
                        energy_esc_upstream, pₓ_esc_upstream,
                        pcut_prev, i_fin, ∑P_downstream, ∑KEdensity_downstream,
                        aa, zz, m, mc,
                        nc_unit, n_cr_count, pmax_cutoff,
                        B_CMBz,
                        weight_new, ptot_pf_new, pb_pf_new, grid_new, downstream_new, inj_new,
                        xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new, tcut_new, x_PT_cm_new,
                        l_save, weight_sav, ptot_pf_sav, pb_pf_sav, grid_sav, downstream_sav, inj_sav,
                        xn_per_sav, prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav, x_PT_cm_sav,
                        use_custom_εB, x_grid_stop,
                        uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
                        γ_ef_grid, β_ef_grid, btot_grid, θ_grid,
                        pxx_flux, pxz_flux, energy_flux, pₓ_esc_feb, energy_esc_feb,
                        do_rad_losses, do_retro, do_tcuts, dont_DSA, dont_scatter, use_custom_frg,
                        inj_fracs, xn_per_fine, xn_per_coarse,
                        x_grid_cm, therm_grid, therm_ptot_sk, therm_pₓ_sk, therm_weight,
                        spectra_pf, spectra_sf, tcuts, pcuts_use, age_max,
                        ε_target, energy_transfer_pool, electron_weight_fac,
                        weight_coupled, spectra_coupled, esc_energy_eff, esc_num_eff, esc_flux,
                        esc_psd_feb_upstream, esc_psd_feb_downstream,
                       )

                    if !l_save[i_prt]
                        particle_finish!(pₓ_esc_feb, energy_esc_feb, esc_energy_eff, esc_num_eff,
                                         esc_flux, esc_psd_feb_downstream, esc_psd_feb_upstream,
                                         i_reason, i_iter, i_ion,
                                         num_psd_mom_bins, num_psd_θ_bins,
                                         psd_bins_per_dec_mom, psd_bins_per_dec_θ,
                                         psd_mom_min, psd_θ_min, psd_cos_fine, Δcos,
                                         aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad,
                                         uₓ_sk, uz_sk, utot,
                                         γᵤ_sf, b_cosθ, b_sinθ, weight, mc)
                    end

                    # Particle counting
                    if i_cut == 1 && (i_prt == 1 || i_prt % n_print_pt == 0)
                        @info("Particle = $i_prt")
                    end

                    i_fin += 1
                    #if i_fin % 16 == 0
                    #    print_progress_bar(i_fin, n_pts_use)
                    #end

                end # loop_pt
                #$omp end parallel do
                #--------------------------------------------------------------------
                # Conclusion of particle loop

                break_pcut, n_pts_target, n_saved = pcut_finalize(
                     i_iter, i_ion, i_cut, p_pcut_hi, n_pts_pcut, n_pts_pcut_hi, n_pts_use,
                     weight_running, l_save, t_start, pcuts_in, pcuts_use, outfile)
                if break_pcut
                    break
                end
                (
                 grid_new, tcut_new, downstream_new, inj_new, weight_new, ptot_pf_new, pb_pf_new,
                 x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new,
                 n_pts_use, weight_running
                ) = new_pcut(
                    n_pts_target, n_saved, l_save, grid_sav, downstream_sav, inj_sav,
                    weight_sav, ptot_pf_sav, pb_pf_sav, x_PT_cm_sav, xn_per_sav,
                    prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav,
                    n_pts_use, weight_running)


            end # loop_pcut
            #----------------------------------------------------------------------
            # Conclusion of pcuts loop

            (;
                dNdp_therm, dNdp_therm_pvals, dNdp_cr, zone_pop,
                pressure_psd_par, pressure_psd_perp, energy_density_psd,
            ) = ion_finalize(
                outfile, nc_unit,
                aa, esc_psd_feb_upstream, esc_psd_feb_downstream,
                jet_rad_pc, jet_sph_frac, m_ion, aa_ion, zz_ion, T₀_ion, n₀_ion,
                (u₀, β₀, γ₀), u₂, n_ions,
                do_multi_dNdps, do_photons,
                mc,
                n_grid, x_grid_cm, uₓ_sk_grid,
                i_iter,
                i_ion,
                γ_sf_grid,
                therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight, num_crossings, n_cr_count,
                num_psd_mom_bins, psd_mom_bounds,
                psd, psd_lin_cos_bins, num_psd_θ_bins, psd_θ_bounds,
                psd_bins_per_dec_mom, psd_mom_min, psd_bins_per_dec_θ, psd_cos_fine, Δcos, psd_θ_min,
                n_shell_endpoints,
                flux_px_upstream, flux_energy_upstream, btot_grid,
                zone_vol, therm_energy_density, energy_density,
            )

        end # loop_ion
        #------------------------------------------------------------------------
        # Conclusion of species loop

        #DEBUGLINE
        for i in 1:n_grid
            @debug("",
                   i_iter, i,
                   #x_grid_rg[i],
                   therm_energy_density[i,1:n_ions]/zone_vol[i],
                   therm_energy_density[i,n_ions]/max(sum(therm_energy_density[i,1:n_ions]),1e-99),
                   sum(therm_energy_density[i,1:n_ions]/zone_vol[i]),
                   energy_density[i,1:n_ions]/zone_vol[i],
                   energy_density[i,n_ions]/max(sum(energy_density[i,1:n_ions]),1e-99),
                   sum(energy_density[i,1:n_ions]/zone_vol[i]))
        end
        #print_plot_vals(888)

        iter_finalize(
            i_iter, i_shock, n_grid, outfile, Γ₂_RH, x_grid_cm, x_grid_rg,
            Γ_grid, uz_sk_grid, θ_grid,
            pxx_flux, energy_flux, pₓ_esc_flux_upstream, pₓ_esc_upstream, flux_px_upstream,
            energy_esc_flux_upstream, energy_esc_upstream, energy_density_psd,
            flux_energy_upstream,
            pressure_psd_par, pressure_psd_perp,
            Γ_downstream, ∑P_downstream, ∑KEdensity_downstream,
            q_esc_cal_pₓ, q_esc_cal_energy,
            uₓ_sk_grid, γ_sf_grid, btot_grid, utot_grid,
            γ_ef_grid, β_ef_grid, εB_grid,
            r_comp, r_RH, u₀, β₀, γ₀, species, γ₂, β₂, u₂,
        )
        # If tcut tracking was enabled, print out the results here
        if do_tcuts
            tcut_print(weights_file, spectra_file, n_tcuts, tcuts, n_ions,
                       num_psd_mom_bins, psd_mom_bounds,
                       i_iter, weight_coupled, spectra_coupled)
        end

    end # loop_itr
    #--------------------------------------------------------------------------
    # Conclusion of iteration loop

end
