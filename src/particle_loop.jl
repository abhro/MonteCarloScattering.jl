function particle_loop(
        i_iter::Int, i_ion::Int, i_cut::Int, i_prt::Int, i_grid_feb::Int, i_shock::Int,
        n_ions::Int, n_pcuts::Int, n_pts_max::Int, n_xspec::Int, n_grid::Int,
        n_print_pt::Int, num_psd_mom_bins::Int, num_psd_θ_bins::Int,
        psd_cos_fine, Δcos, psd_θ_min,
        species,
        γ₀::Float64, β₀::Float64, u₀, u₂, bmag₂,
        pₑ_crit, γₑ_crit, η_mfp,
        psd_mom_min, psd_bins_per_dec_mom, psd_bins_per_dec_θ,
        energy_transfer_frac, energy_recv_pool,
        psd, num_crossings, x_spec, feb_upstream, feb_downstream,
        energy_esc_upstream, pₓ_esc_upstream, pcut_prev, i_fin, ∑P_downstream,
        ∑KEdensity_downstream,
        aa, zz, m, mc,
        nc_unit, n_cr_count, pmax_cutoff,
        B_CMBz,
        weight_new, ptot_pf_new, pb_pf_new, grid_new, downstream_new, inj_new,
        xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new, tcut_new, x_PT_cm_new,
        l_save, weight_sav, ptot_pf_sav, pb_pf_sav, grid_sav, downstream_sav, inj_sav,
        xn_per_sav, prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav, x_PT_cm_sav,
        use_custom_εB, x_grid_stop,
        uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid, γ_ef_grid, β_ef_grid, btot_grid, θ_grid,
        pxx_flux, pxz_flux, energy_flux, pₓ_esc_feb, energy_esc_feb,
        do_rad_losses, do_retro, do_tcuts, dont_DSA, dont_scatter, use_custom_frg,
        inj_fracs, xn_per_fine, xn_per_coarse,
        x_grid_cm, therm_grid, therm_ptot_sk, therm_pₓ_sk, therm_weight,
        spectra_pf, spectra_sf, tcuts, pcuts_use, age_max,
        ε_target, energy_transfer_pool, electron_weight_fac,
        weight_coupled, spectra_coupled, esc_energy_eff, esc_num_eff, esc_flux,
        esc_psd_feb_upstream, esc_psd_feb_downstream,
    )
    # To maintain identical results between OpenMP and serial versions,
    # set RNG seed based on current iteration/ion/pcut/particle number
    iseed_mod =  (  (i_iter - 1)*n_ions*n_pcuts*n_pts_max
                  + (i_ion - 1)        *n_pcuts*n_pts_max
                  + (i_cut - 1)                *n_pts_max
                  +  i_prt)
    Random.seed!(iseed_mod)

    # Reset the counter for number of times through the main loop
    helix_count = 0


    # Get the properties of the particle we're about to treat
    weight       = weight_new[i_prt]
    ptot_pf      = ptot_pf_new[i_prt]
    pb_pf        = pb_pf_new[i_prt]
    i_grid       = grid_new[i_prt]
    i_grid_old   = i_grid  # Needed for energy transfer
    l_downstream = downstream_new[i_prt]
    inj          = inj_new[i_prt]
    xn_per       = xn_per_new[i_prt]
    prp_x_cm     = prp_x_cm_new[i_prt]
    acctime_sec  = acctime_sec_new[i_prt]
    φ_rad        = φ_rad_new[i_prt]
    tcut_curr    = tcut_new[i_prt]

    r_PT_cm = SVector{3,LengthCGS}(x_PT_cm_new[i_prt], # x
                                   # Not currently tracked, but could be added in at later date
                                   0.0cm, # y
                                   0.0cm) # z

    γₚ_pf = hypot(1, ptot_pf/mc)

    # Constant that will be used repeatedly during loop
    # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
    gyro_denom = 1 / (zz * btot_grid[i_grid]) # zz already has units of charge
    if use_custom_εB && r_PT_cm.x > x_grid_stop
        gyro_denom *= √(r_PT_cm.x / x_grid_stop)
    end

    # Gyroradius assuming all motion is perpendicular to B field;
    # pitch-angle-correct gyroradius is gyro_rad_cm
    gyro_rad_tot_cm = ptot_pf * c * gyro_denom |> cm

    # Gyroperiod in seconds
    gyro_period_sec = 2π * γₚ_pf * aa*mp * c * gyro_denom


    # Get the properties of the grid zone the particle's in
    uₓ_sk = uₓ_sk_grid[i_grid]
    uz_sk = uz_sk_grid[i_grid]
    utot  =  utot_grid[i_grid]
    γᵤ_sf =  γ_sf_grid[i_grid]
    γᵤ_ef =  γ_ef_grid[i_grid]
    βᵤ_ef =  β_ef_grid[i_grid]
    bmag  =  btot_grid[i_grid]
    bθ    =     θ_grid[i_grid]

    b_sinθ = sin(bθ)
    b_cosθ = cos(bθ)


    #------------------------------------------------------------------
    # Helix loop: The following loop is a large-scale restructuring of the original code's
    # loops 5702 and 2001.
    # The original code was slightly clearer in intent, but the goto statements made it
    # impossible to parallelize for CPU computing; this code has sacrificed some of the
    # clarity in exchange for potential speed gains down the road.
    #
    #    ORIGINAL                       THIS CODE
    #  5702 continue                 first_iter = true
    #
    #   Code Block 1                 while keep_looping
    #                                    if first_iter || condition 1
    #  2001 continue                         first_iter = false
    #                                        Code Block 1
    #   Code Block 2                     else
    #                                        Code Block 3
    #  if (condition 1) goto 5702        end
    #
    #   Code Block 3                     Code Block 2
    #                                end
    #  goto 2001
    #
    # Code Block 1: start of helix loop; minor setup
    # Code Block 2: particle movement; check on downstream PRP return
    # Code Block 3: grid zone changes; non-downstream escape (exception: if no
    #   scattering); radiative losses; momentum splitting; scattering
    #
    # loop_helix in this code is the old loop 2001. The first part of loop 5702 has been
    # subsumed into the if statement, which determines whether Code Block 3 would be run or
    # whether the code would cycle back to run Code Block 1 again. Following every
    # combination of logical choices should result in the same path through the two loops.
    #------------------------------------------------------------------
    keep_looping = true
    first_iter   = true
    # Control variable for "condition 1":
    #  -1: default
    #   0: particle escapes
    #   1: particle returns from downstream via convection
    #   2: particle didn't enter return calculations
    i_return = -1
    i_reason = 0
    lose_pt  = false

    local t_step, p_perp_b_pf, gyro_rad_cm
    r_PT_old = SVector{3,LengthCGS}(0cm, 0cm, 0cm)
    while keep_looping # loop_helix

        # Track number of times through the main loop. This will only be
        # needed for electrons at high energies when using radiative losses
        helix_count += 1

        if first_iter || i_return == 1
            first_iter = false

            # Code Block 1: start of helix loop; minor setup
            #--------------------------------------------------------------
            # Calculate momentum perpendicular to magnetic field...
            p_perp_b_pf = perpendicular_momentum(ptot_pf, pb_pf)

            # ...and turn that into a gyroradius
            gyro_rad_cm = p_perp_b_pf * c * gyro_denom |> cm
            #--------------------------------------------------------------
            # End of Code Block 1

        else

            # Code Block 3: grid zone changes; non-downstream escape (exception:
            # if no scattering); radiative losses; momentum splitting; scattering
            #--------------------------------------------------------------
            # Store values from the previous run through loop_helix; only
            # values needed by subroutine transform_p_PSP are kept here
            uₓ_sk_old = uₓ_sk
            uz_sk_old = uz_sk
            utot_old  = utot
            γᵤ_sf_old = γᵤ_sf
            b_sin_old = b_sinθ
            b_cos_old = b_cosθ


            # Get new values for this run through loop_helix
            uₓ_sk  = uₓ_sk_grid[i_grid]
            uz_sk  = uz_sk_grid[i_grid]
            utot   =  utot_grid[i_grid]
            γᵤ_sf  =  γ_sf_grid[i_grid]
            γᵤ_ef  =  γ_ef_grid[i_grid]
            βᵤ_ef  =  β_ef_grid[i_grid]
            bmag   =  btot_grid[i_grid]
            bθ     =     θ_grid[i_grid]
            b_sinθ = sin(bθ)
            b_cosθ = cos(bθ)

            # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
            if use_custom_εB && (r_PT_cm.x > x_grid_stop)
                bmag = btot_grid[n_grid] * √(x_grid_stop / r_PT_cm.x)
            end

            gyro_denom = 1 / (zz * bmag)

            # If particle crossed a velocity gradient, find its new shock frame properties
            if uₓ_sk != uₓ_sk_old
                (
                 ptot_pf, ptot_sk, p_sk, pb_sk, p_perp_b_sk, γₚ_sk,
                 pb_pf, p_perp_b_pf, γₚ_pf, φ_rad
                ) = transform_p_PSP(
                                    aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad,
                                    uₓ_sk_old, uz_sk_old, utot_old, γᵤ_sf_old,
                                    b_cos_old, b_sin_old, uₓ_sk, uz_sk, utot, γᵤ_sf,
                                    b_cosθ, b_sinθ, mc)

                gyro_rad_cm = p_perp_b_pf * c * gyro_denom |> cm
                gyro_rad_tot_cm = ptot_pf * c * gyro_denom |> cm
                pₓ_sk, py_sk, pz_sk = p_sk
            end


            # Energy transfer from ions to electrons in the upstream region; only bother with this if
            #  1. Electrons are a separate species and energy transfer is enabled
            #  2. Particles are still on their first trip to the downstream region of the shock structure
            #  3. Particles have entered a new grid zone, and energy should be either subtracted or added
            # TODO factor this out
            if energy_transfer_frac > 0 && !inj && r_PT_old.x ≤ 0cm && i_grid_old != i_grid
                i_start = i_grid_old
                i_stop  = min(i_grid, i_shock)

                # Subtract energy from the ions and add it to the pool of energy for this grid zone.
                #@debug "" i_start i_stop i_grid i_shock
                if aa ≥ 1 && maximum(ε_target[i_start+1:i_stop]) > 0

                    # Subtract energy based on the difference between current
                    # and previous values of ε_electron
                    γ_pf_i = hypot(1, ptot_pf/mc)
                    γ_pf_f = 1 + (γ_pf_i - 1) * (1-ε_target[i_stop])/(1-ε_target[i_start])

                    # Split the energy equally among all grid zones crossed during
                    # this scattering step; the donated energy is weighted by weight.
                    n_split = count(ε_target[i_start+1:i_stop] .> 0)
                    #$omp critical
                    energy_increment = (γ_pf_i - γ_pf_f) * aa * E₀ₚ * weight / n_split
                    for i in i_start+1:i_stop
                        if ε_target[i] > 0
                            energy_transfer_pool[i] += energy_increment
                        end
                    end
                    #$omp end critical

                    # Calculate the new momentum based on the new energy,
                    # and rescale components accordingly
                    ptot_pf_f = aa*mp*c * √(γ_pf_f^2 - 1) |> (g*cm/s)
                    scale_fac = ptot_pf_f / ptot_pf

                    pb_pf *= scale_fac
                    p_perp_b_pf *= scale_fac
                    ptot_pf = ptot_pf_f
                    γₚ_pf = γ_pf_f

                elseif maximum(energy_recv_pool[i_start+1:i_stop], init=0erg) > 0erg

                    # For electrons, add pooled energy. Include energy from all cells electron
                    # crossed in this scattering step. Also modify the amount of energy to reflect
                    # the number of electrons each MC particle represents
                    energy_to_transfer = sum(@view(energy_recv_pool[i_start+1:i_stop])) * electron_weight_fac

                    γ_pf_i = hypot(1, ptot_pf/mc)
                    γ_pf_f = γ_pf_i + energy_to_transfer/(aa*E₀ₚ)

                    # Calculate the new momentum based on the new energy, and
                    # rescale components accordingly
                    ptot_pf_f = aa*mp*c * √(γ_pf_f^2 - 1) |> (g*cm/s)
                    scale_fac = ptot_pf_f / ptot_pf

                    pb_pf *= scale_fac
                    p_perp_b_pf *= scale_fac
                    ptot_pf = ptot_pf_f
                    γₚ_pf = γ_pf_f

                end  # check on particle species

                # Since the plasma-frame momenta have changed,
                # recalculate the shock-frame momenta
                ptot_sk, p_sk, γₚ_sk = transform_p_PS(aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad,
                                                      uₓ_sk, uz_sk, utot, γᵤ_sf,
                                                      b_cosθ, b_sinθ, mc)
                pₓ_sk = p_sk.x
                pz_sk = p_sk.z
            end
            # check on grid zone
            # end check on whether particles are thermal and upstream
            # end check on energy_transfer_frac

            # Particle escape: downstream with scattering disabled
            if dont_scatter && r_PT_cm.x > 10gyro_rad_cm
                i_return = 0
                i_reason = 1

                keep_looping = false
                continue
            end

            # Particle escape: pmax (note that effects are same as for escape at upstream FEB)
            if ptot_pf > pmax_cutoff
                # Transform plasma frame momentum into shock frame to test there also
                ptot_sk, p_sk, γₚ_sk = transform_p_PS(aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad,
                                                      uₓ_sk, uz_sk, utot, γᵤ_sf, b_cosθ, b_sinθ, mc)

                if ptot_sk > pmax_cutoff
                    i_reason = 2

                    keep_looping = false
                    continue
                end  # check on ptot_sk
            end  # check on ptot_pf

            # Particle escape: upstream FEB (note that effects are same as for escape by pmax)
            if inj && r_PT_cm.x < feb_upstream
                i_reason = 2

                keep_looping = false
                continue
            end

            # Particle escape: age_max
            if age_max > 0s && acctime_sec > age_max
                i_reason = 3

                keep_looping = false
                continue
            end

            # Particle escape: transverse motion
            #TODO: implement new keyword related to transverse motion
            #if yz_pos_max > 0 && hypot(y_PT, z_PT) > yz_pos_max
            #    i_reason = 2
            #    keep_looping = false
            #    continue
            #end

            # Synchrotron/IC/CMB losses for electrons
            if do_rad_losses && aa < 1
                # Store previous ptot_pf value
                ptot_pf_old = ptot_pf

                ptot_pf = electron_radiation_loss(B_CMBz, γᵤ_ef, bmag, ptot_pf, t_step)

                # Catch electrons that have somehow lost all their energy in a single time step
                if ptot_pf ≤ 0g*cm/s
                    ptot_pf = 1e-99g*cm/s
                    pb_pf = 1e-99g*cm/s
                    p_perp_b_pf = 1e-99g*cm/s
                    γₚ_pf = 1
                    i_reason = 4
                    keep_looping = false
                    continue
                end

                # Recalculate γₚ_pf since that's needed elsewhere
                γₚ_pf = hypot(ptot_pf/mc, 1)

                # Modify components of ptot_pf due to losses
                pb_pf       *= ptot_pf/ptot_pf_old
                p_perp_b_pf *= ptot_pf/ptot_pf_old

                # Also recalculate gyroradii
                gyro_rad_tot_cm =     ptot_pf * c * gyro_denom |> cm
                gyro_rad_cm     = p_perp_b_pf * c * gyro_denom |> cm
            end


            #---------------------------------------------
            # -- SCATTERING -- SCATTERING -- SCATTERING --
            #---------------------------------------------
            if !dont_scatter
                gyro_period_sec, pb_pf, p_perp_b_pf, φ_rad = scattering(
                    aa, gyro_denom, ptot_pf, γₚ_pf, xn_per, pb_pf, p_perp_b_pf, φ_rad,
                    use_custom_frg, pₑ_crit, γₑ_crit, η_mfp,
                )
            end

            # Update acceleration time in explosion frame, so convert t_step from plasma frame
            # Only start the clock once a particle has crossed the shock for the first time
            if l_downstream
                acctime_sec += t_step*γᵤ_ef

                if do_tcuts && acctime_sec ≥ tcuts[tcut_curr]
                    tcut_track!(weight_coupled, spectra_coupled, tcut_curr,
                                weight, ptot_pf, i_ion, num_psd_mom_bins,
                                psd_mom_min, psd_bins_per_dec_mom)
                    tcut_curr += 1
                end

                # Remove ions at splitting momentum
                if ptot_pf > pcuts_use[i_cut]
                    @info("Removing ions at splitting momentum",
                          pcuts_use[i_cut], i_cut, ptot_pf)
                    l_save[i_prt] = true

                    weight_sav[i_prt]       = weight
                    ptot_pf_sav[i_prt]      = ptot_pf
                    pb_pf_sav[i_prt]        = pb_pf
                    x_PT_cm_sav[i_prt]      = r_PT_cm.x
                    grid_sav[i_prt]         = i_grid
                    downstream_sav[i_prt]   = l_downstream
                    inj_sav[i_prt]          = inj
                    xn_per_sav[i_prt]       = xn_per
                    prp_x_cm_sav[i_prt]     = r_PT_cm.x < prp_x_cm ? prp_x_cm : r_PT_cm.x * 1.1cm # Ensure particle is within PRP
                    acctime_sec_sav[i_prt]  = acctime_sec
                    φ_rad_sav[i_prt]        = φ_rad
                    tcut_sav[i_prt]         = tcut_curr

                    keep_looping = false
                    continue
                end
            end

            # Shift between coarse/fine xn_per; right now the code uses only one value of
            # xn_per everywhere in the ptot/x_PT phase space, so this branch does nothing.
            xn_per = r_PT_cm.x > gyro_rad_tot_cm ? xn_per_coarse : xn_per_fine
            #--------------------------------------------------------------
            # End of Code Block 3

        end  # check on i_return


        # Code Block 2: particle movement; fluxes; downstream escape/return
        #---------------------------------------------------------------
        # Store old x/y/z/φ positions
        r_PT_old  = r_PT_cm
        φ_rad_old = φ_rad


        # Find time step
        t_step = gyro_period_sec / xn_per |> s

        # Odd little loop that only matters if DSA has been disabled per the input file
        #------------------------------------------------------------------------------
        (; pb_pf, φ_rad, x_move_bpar, r_PT_cm) = no_DSA_loop(
            φ_rad, xn_per, pb_pf, t_step, γₚ_pf, aa, mp, dont_DSA, inj_fracs,
            r_PT_old, r_PT_cm, b_cosθ, b_sinθ, γᵤ_sf, φ_rad_old, uₓ_sk, inj, i_ion, gyro_rad_cm)
        #----------------------------------------------------------------
        # DSA injection prevented if specified


        # If particle crosses shock going upstream -> downstream
        if r_PT_old.x < 0cm && r_PT_cm.x ≥ 0cm
            l_downstream = true

            # Ensure the downstream region is sufficiently long to allow particle to isotropize.
            # This simple equation is decidedly non-trivial, and comes from two assumptions:
            # 1. The particle's diffusion coefficient D may be described by D = ⅓⋅η_mfp⋅r_g⋅v_pt
            #    (by default; a different f(r_g) may be specified as desired in place of η⋅r_g)
            # 2. The relation between the diffusion coefficient D and the diffusion length L is
            #    L = D/⟨u⟩, where ⟨u⟩ is the average speed of diffusion. Assuming isotropic
            #    particles in the downstream frame, ⟨u⟩ = u₂ since the average thermal *velocity* of
            #    the population is 0.
            #TODO: include f(r_g) in place of η*r_g to allow for arbitrary diffusion
            L_diff = η_mfp/3 * gyro_rad_tot_cm * ptot_pf/(aa*mp*γₚ_pf * u₂) |> cm

            @info("Particle crossing shock going upstream → downstream",
                  prp_x_cm, r_PT_old.x, r_PT_cm.x, L_diff)
            prp_x_cm = max(prp_x_cm, L_diff)
        end


        # Test for injection into the acceleration process
        if l_downstream && r_PT_cm.x < 0cm
            inj = true
        end


        # Calculate fluxes due to this motion; also locates new grid zone
        (i_grid, i_grid_old, n_cr_count, pₓ_esc_upstream, energy_esc_upstream) = all_flux!(
            i_prt, aa, pb_pf, p_perp_b_pf, ptot_pf, γₚ_pf, φ_rad,
            weight, i_grid, uₓ_sk, uz_sk, utot, γᵤ_sf, b_cosθ, b_sinθ,
            r_PT_cm.x, r_PT_old.x, inj, nc_unit,
            i_grid_feb, pxx_flux, pxz_flux,
            energy_flux, energy_esc_upstream, pₓ_esc_upstream, spectra_sf, spectra_pf,
            n_cr_count, num_crossings, psd,
            n_xspec, x_spec, feb_upstream, γ₀, u₀, mc,
            n_grid, x_grid_cm,
            therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight,
            psd_bins_per_dec_mom, psd_bins_per_dec_θ, psd_mom_min,
            num_psd_mom_bins, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min,
        )


        # Downstream escape/return test; assume we'll call subroutine
        # prob_return unless told otherwise by particle location/info
        do_prob_ret, i_return = downstream_test(
            i_return, feb_downstream, gyro_denom, gyro_rad_tot_cm, prp_x_cm, r_PT_cm, aa,
            ptot_pf, pₑ_crit, γₚ_pf, γₑ_crit, u₂, η_mfp)

        if do_prob_ret
            (
             i_return, lose_pt, tcut_curr, _x_PT_cm, prp_x_cm, ptot_pf, γₚ_pf,
             gyro_denom, pb_pf, p_perp_b_pf, acctime_sec, φ_rad
            ) = prob_return(
                            i_ion, num_psd_mom_bins, B_CMBz, r_PT_old.x, aa, zz, gyro_denom,
                            r_PT_cm.x, prp_x_cm, ptot_pf, γₚ_pf, pb_pf, p_perp_b_pf,
                            acctime_sec, φ_rad, helix_count, pcut_prev, weight, tcut_curr,
                            x_grid_stop, u₂, use_custom_εB, η_mfp, do_retro, bmag₂, mc,
                            do_rad_losses, do_tcuts, tcuts,
                            n_grid, uₓ_sk_grid, γ_sf_grid, γ_ef_grid, θ_grid, btot_grid,
                            psd_mom_min, psd_bins_per_dec_mom, weight_coupled, spectra_coupled,
                           )
            r_PT_cm = SVector(_x_PT_cm, r_PT_cm.y, r_PT_cm.z)
        end

        # If particle escaped downstream, handle final calculations here
        if i_return == 0

            vel = ptot_pf / (aa*mp) |> cm/s
            if (γₚ_pf - 1) ≥ E_rel_pt
                vel /= γₚ_pf
            end

            ∑P_downstream += ptot_pf/3 * vel * weight * density(species[i_ion]) # downstream pressure
            ∑KEdensity_downstream += (γₚ_pf - 1) * aa*E₀ₚ * weight * density(species[i_ion])

            i_reason = 1
            if lose_pt
                i_reason = 4
            end

            keep_looping = false
            continue
        end
        #----------------------------------------------------------------
        # End of Code Block 2

    end # loop_helix
    #------------------------------------------------------------------
    # End of loop moving/tracking particles on/off grid
    return (i_fin, i_reason, pb_pf, p_perp_b_pf, γₚ_pf, γᵤ_sf, b_cosθ, b_sinθ, weight,
            uₓ_sk, uz_sk, utot, φ_rad,
            ∑P_downstream, ∑KEdensity_downstream)
end

function no_DSA_loop(
        φ_rad, xn_per, pb_pf, t_step, γₚ_pf, aa, mp, dont_DSA,
        inj_fracs, r_PT_old, r_PT_cm, b_cosθ, b_sinθ, γᵤ_sf, φ_rad_old, uₓ_sk, inj, i_ion,
        gyro_rad_cm)
    local x_move_bpar
    while true # loop_no_DSA

      # Update position/phase angle.
      # Phase angle is easy: add appropriate amount and make sure result is in [0,2π)
      # For position, use inverse Lorentz transformations where primed frame is the
      # plasma frame, moving to the right along the x axis, and the unprimed frame is the
      # (stationary) shock frame. Take frames to be coincident
      # at t = t' = 0. Then
      #    x  =  γᵤ_sf * (x' + uₓ_sk*t'),
      # where t' is t_step and x' is distance moved in plasma
      # frame (i.e. x_move_bpar below).
      # Remember to take gyration about magnetic field into account
      φ_rad = mod2pi(φ_rad + 2π/xn_per)

      x_move_bpar = pb_pf * t_step / (γₚ_pf * aa*mp) |> cm

      Δx_PT_cm = γᵤ_sf * (x_move_bpar * b_cosθ
                          - gyro_rad_cm * b_sinθ * (cos(φ_rad)-cos(φ_rad_old))
                          + uₓ_sk * t_step)
      r_PT_cm = SVector{3,LengthCGS}(r_PT_old.x + Δx_PT_cm, r_PT_cm.y, r_PT_cm.z)

      # Don't care about transverse motion unless it's being tracked
      #TODO: add new keyword for tracking transverse motion
      #if yz_pos_max > 0
      #    y_PT_cm = y_PT_old + gyro_rad_cm * (sin(φ_rad) - sin(φ_rad_old))
      #    z_PT_cm = z_PT_old + γᵤ_sf * (x_move_bpar*b_sinθ - gyro_rad_cm*b_cosθ*(cos(φ_rad)-cos(φ_rad_old)) + uz_sk*t_step)
      #end


      #TODO: set flags for energy transfer here?
      #Search original code for "out of ions"

      # Reflect particles under two conditions:
      #  1. They've been downstream and are crossing back upstream
      #  2. Either DSA is disabled or they fail an injection check
      # continue again if particles get reflected
      if r_PT_cm.x ≤ 0cm && r_PT_old.x > 0cm && !inj && (dont_DSA || inj_fracs[i_ion] < 1)
          if dont_DSA || (Random.rand() > inj_fracs[i_ion])
              # Reflect particles with negative pb_pf; randomize phase of the rest
              if pb_pf < 0g*cm/s
                  pb_pf = -pb_pf
              else
                  φ_rad = Random.rand() * 2π
              end
          else
              break
          end
      else
          break
      end
    end # loop_no_DSA
    return (; pb_pf, φ_rad, x_move_bpar, r_PT_cm)
end

"""
    electron_radiation_loss(B_CMBz, γᵤ_ef, bmag, p, Δt)

TODO
"""
function electron_radiation_loss(B_CMBz, γᵤ_ef, bmag, p, Δt)
    # Compute effective magnetic field for radiative losses. When squared to get
    # energy density, one factor of γᵤ_ef in B_CMB_loc represents increased number
    # density, while the other is increased energy per photon.
    B_CMB_loc = B_CMBz * γᵤ_ef

    # Note that here Δln_p_synch = Δp/p. If this value is too large we will
    # directly integrate from p_i to get p_f, since the discrete approach would
    # result in too high a loss in a single time step
    Δln_p_synch = rad_loss_fac * (bmag^2 + B_CMB_loc^2) * p * Δt |> NoUnits

    # Correction to make sure electrons don't lose too much energy in a single time step
    if Δlnp_synch > 1e-2
        p /= 1 + Δln_p_synch
    else
        p *= 1 - Δln_p_synch # Put second factor of ptot_pf back into dp_synch
    end
    return p
end


function downstream_test(
        i_return_orig, feb_downstream, gyro_denom, gyro_rad_tot_cm,
        prp_x_cm, r_PT_cm, aa, ptot_pf, pₑ_crit, γₚ_pf, γₑ_crit, u₂, η_mfp)
    do_prob_ret = true

    local i_return = i_return_orig

    if feb_downstream > 0cm && r_PT_cm.x > feb_downstream

        # Particle flagged as escaping downstream
        i_return    = 0
        do_prob_ret = false

    elseif r_PT_cm.x > 1.1*prp_x_cm

        # The following simple equation is decidedly non-trivial, and comes from two
        # assumptions:
        # 1. The particle's diffusion coefficient D may be described by D = ⅓⋅η_mfp⋅r_g⋅v_pt
        #    (by default; a different f(r_g) may be specified as desired in place of η⋅r_g)
        # 2. The relation between the diffusion coefficient D and the diffusion length L is
        #    L = D/⟨u⟩, where ⟨u⟩ is the average speed of diffusion. Assuming isotropic
        #    particles in the downstream frame, ⟨u⟩ = u₂ since the average thermal *velocity* of
        #    the population is 0.
        #TODO: include f(r_g) in place of η*r_g to allow for arbitrary diffusion
        if aa < 1 && ptot_pf < pₑ_crit
            gyro_fac = pₑ_crit * c * gyro_denom
            v_fac = gyro_fac        * pₑ_crit / (aa*mp*γₑ_crit * u₂)
        else
            v_fac = gyro_rad_tot_cm * ptot_pf / (aa*mp*γₚ_pf   * u₂)
        end
        L_diff = η_mfp/3 * v_fac

        if r_PT_cm.x > 6.91 * L_diff
            i_return    = 0
            do_prob_ret = false
        end

    end
    return do_prob_ret, i_return
end

function perpendicular_momentum(ptot_pf, pb_pf)
    if ptot_pf < abs(pb_pf)
        p_perp_b_pf = 1e-6 * ptot_pf
        pb_pf = copysign(√(ptot_pf^2 - p_perp_b_pf^2), pb_pf)
        # XXX CHECKTHIS: does this *ever* happen?!
        @warn("ptot_pf < pb_pf at top of loop_helix", ptot_pf, pb_pf)
    else
        p_perp_b_pf = √(ptot_pf^2 - pb_pf^2)
    end

    return p_perp_b_pf
end
