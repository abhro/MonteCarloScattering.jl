function helix_loop!(
        helix_count, keep_looping,
        first_iter, n_cr_count, i_return, i_prt, i_grid, i_grid_old,
        pt_pf, pb_pf, γ_pt_pf,
        γ_usf,
        gyro_denom, gyro_period_sec,
        r_PT_cm, φ_rad, xn_per,
        ux_sk, uz_sk, utot,
        b_cosθ, b_sinθ,
        aa, zz, o_o_mc,
        prp_x_cm, acctime_sec, pcut_prev, tcut_curr, pmax_cutoff,
        l_DwS, inj,
        xwt, energy_esc_UpS, px_esc_UpS,
        nc_unit,
    )

    # Track number of times through the main loop. This will only be
    # needed for electrons at high energies when using radiative losses
    helix_count += 1

    if first_iter || i_return == 1
        first_iter = false

        # Code Block 1: start of helix loop; minor setup
        #--------------------------------------------------------------
        # Calculate momentum perpendicular to magnetic field...
        if pt_pf < abs(pb_pf)
            global p_perp_b_pf = 1e-6 * pt_pf
            pb_pf = copysign( √( pt_pf^2 - p_perp_b_pf^2 ), pb_pf )
            #CHECKTHIS: does this *ever* happen?!
            @warn("pt_pf < pb_pf at top of loop_helix")
        else
            global p_perp_b_pf = √( pt_pf^2  -  pb_pf^2 )
        end

        # ...and turn that into a gyroradius
        global gyro_rad_cm = p_perp_b_pf * c_cgs * gyro_denom
        #--------------------------------------------------------------
        # End of Code Block 1

    else

        # Code Block 3: grid zone changes; non-DwS escape (exception:
        # if no scattering); radiative losses; momentum splitting; scattering
        #--------------------------------------------------------------
        # Store values from the previous run through loop_helix; only
        # values needed by subroutine transform_p_PSP are kept here
        ux_sk_old  = ux_sk
        uz_sk_old  = uz_sk
        utot_old  = utot
        γ_usf_old = γ_usf
        b_sin_old = b_sinθ
        b_cos_old = b_cosθ


        # Get new values for this run through loop_helix
        ux_sk   = ux_sk_grid[i_grid]
        uz_sk   = uz_sk_grid[i_grid]
        utot   = utot_grid[i_grid]
        γ_usf  = γ_sf_grid[i_grid]
        γ_uef  = γ_ef_grid[i_grid]
        β_uef  = β_ef_grid[i_grid]
        bmag   = btot_grid[i_grid]
        bθ     = θ_grid[i_grid]
        b_sinθ = sin(bθ)
        b_cosθ = cos(bθ)

        # Exponent of 0.5 corresponds to Blandford-McKee solution,
        # where e ∝ 1/χ ∝ 1/r
        if use_custom_εB && (r_PT_cm.x > x_grid_stop)
            bmag = btot_grid[n_grid] * √(x_grid_stop / r_PT_cm.x)
        end

        gyro_denom = 1 / (zz*qp_cgs * bmag)

        # If particle crossed a velocity gradient, find its new shock frame properties
        if ux_sk != ux_sk_old
            pt_pf, pt_sk, px_sk, py_sk, pz_sk, pb_sk, p_perp_b_sk, γ_pt_sk, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad = transform_p_PSP(
                aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad,
                ux_sk_old, uz_sk_old, utot_old, γ_usf_old,
                b_cos_old, b_sin_old, ux_sk, uz_sk, utot, γ_usf,
                b_cosθ, b_sinθ,
                o_o_mc)

            global gyro_rad_cm = p_perp_b_pf * c_cgs * gyro_denom
            global gyro_rad_tot_cm = pt_pf * c_cgs * gyro_denom
        end


        # Energy transfer from ions to electrons in the UpS region; only bother with this if
        #  1. Electrons are a separate species and energy transfer is enabled
        #  2. Particles are still on their first trip to the DwS region of the shock structure
        #  3. Particles have entered a new grid zone, and energy should be either subtracted or added
        if sc_electron && energy_transfer_frac > 0 && !inj && r_PT_old.x ≤ 0 && i_grid_old != i_grid

            i_start = i_grid_old
            i_stop  = min(i_grid, i_shock)

            # Subtract energy from the ions and add it to the pool of energy for this grid zone.
            @info "" i_start i_stop i_grid i_shock
            if aa ≥ 1 && maximum(ε_target[i_start+1:i_stop]) > 0

                # Subtract energy based on the difference between current
                # and previous values of ε_electron
                γ_pf_i = √( 1  +  (pt_pf*o_o_mc)^2 )
                γ_pf_f = 1 + (γ_pf_i - 1) * (1-ε_target[i_stop])/(1-ε_target[i_start])

                # Split the energy equally among all grid zones crossed during
                # this scattering step; the donated energy is weighted by xwt.
                n_split = count( ε_target[i_start+1:i_stop] .> 0 )
                #$omp critical
                for i in i_start+1:i_stop
                    if ε_target[i] > 0
                        energy_transfer_pool[i] += (γ_pf_i - γ_pf_f) * aa * E₀_proton * xwt / n_split
                    end
                end
                #$omp end critical

                # Calculate the new momentum based on the new energy,
                # and rescale components accordingly
                pt_pf_f   = aa*mp_cgs * c_cgs * √( γ_pf_f^2 - 1 )
                scale_fac = pt_pf_f / pt_pf

                pb_pf *= scale_fac
                p_perp_b_pf *= scale_fac
                pt_pf = pt_pf_f
                γ_pt_pf = γ_pf_f

            elseif maximum(energy_recv_pool[i_start+1:i_stop], init=0) > 0

                # For electrons, add pooled energy. Include energy from
                # all cells electron crossed in this scattering step.
                # Also modify the amount of energy to reflect the number
                # of electrons each MC particle represents
                energy_to_transfer = sum( energy_recv_pool[i_start+1:i_stop] )
                energy_to_transfer = energy_to_transfer * electron_wt_fac

                γ_pf_i = √( 1  +  (pt_pf*o_o_mc)^2 )
                γ_pf_f = γ_pf_i  +  energy_to_transfer / (aa*E₀_proton)

                # Calculate the new momentum based on the new energy, and
                # rescale components accordingly
                pt_pf_f   = aa*mp_cgs * c_cgs * √( γ_pf_f^2 - 1 )
                scale_fac = pt_pf_f / pt_pf

                pb_pf *= scale_fac
                p_perp_b_pf *= scale_fac
                pt_pf = pt_pf_f
                γ_pt_pf = γ_pf_f

            end  # check on particle species

            # Since the plasma-frame momenta have changed,
            # recalculate the shock-frame momenta
            pt_sk, px_sk, pz_sk, γ_pt_sk = transform_p_PS(
                aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad,
                ux_sk, uz_sk, utot, γ_usf,
                b_cosθ, b_sinθ, oblique, o_o_mc)
        end
        # check on grid zone
        # end check on whether particles are thermal and UpS
        # end check on sc_electron and energy_transfer_frac


        # Particle escape: DwS with scattering disabled
        if dont_scatter && r_PT_cm.x > 10gyro_rad_cm
            i_return = 0
            i_reason = 1

            keep_looping = false
            return keep_looping, helix_count, first_iter, i_grid_old
        end


        # Particle escape: pmax (note that effects are same as for escape at UpS FEB)
        if pt_pf > pmax_cutoff
            # Transform plasma frame momentum into shock frame to test there also
            pt_sk, px_sk, pz_sk, γ_pt_sk = transform_p_PS(
                aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad,
                ux_sk, uz_sk, utot, γ_usf, b_cosθ, b_sinθ,
                oblique, o_o_mc)

            if pt_sk > pmax_cutoff
                i_reason = 2

                keep_looping = false
                return keep_looping, helix_count, first_iter, i_grid_old
            end  # check on pt_sk
        end  # check on pt_pf


        # Particle escape: UpS FEB (note that effects are same as for escape by pmax)
        if inj && r_PT_cm.x < feb_UpS
            i_reason = 2

            keep_looping = false
            return keep_looping, helix_count, first_iter, i_grid_old
        end


        # Particle escape: age_max
        if age_max > 0 && acctime_sec > age_max
            i_reason = 3

            keep_looping = false
            return keep_looping, helix_count, first_iter, i_grid_old
        end


        # Particle escape: transverse motion
        #TODO: implement new keyword related to transverse motion
        #if yz_pos_max > 0 && hypot(y_PT, z_PT) > yz_pos_max
        #    i_reason = 2
        #    keep_looping = false
        #    return keep_looping, helix_count, first_iter, i_grid_old
        #end


        # Synchrotron/ICCMB losses for electrons
        if do_rad_losses && aa < 1

            # Store previous pt_pf value
            pt_pf_old = pt_pf

            # Compute effective magnetic field for radiative losses. When squared to get
            # energy density, one factor of γ_uef in B_CMB_loc represents increased number
            # density, while the other is increased energy per photon.
            B_CMB_loc = B_CMBz * γ_uef

            B_tot_sq = bmag^2 + B_CMB_loc^2

            # Note that here dp_synch is actually dp/p. If this value is too large we will
            # directly integrate from p_i to get p_f, since the discrete approach would
            # result in too high a loss in a single time step
            dp_synch = rad_loss_fac * B_tot_sq * pt_pf * t_step

            # Correction to make sure electrons don't lose too much
            # energy in a single time step
            if dp_synch > 1e-2
                pt_pf = pt_pf / (1 + dp_synch)
            else
                pt_pf = pt_pf * (1 - dp_synch) # Put second factor of pt_pf back into dp_synch
            end

            # Catch electrons that have somehow lost all their energy in a single time step
            if pt_pf ≤ 0
                pt_pf = 1e-99
                pb_pf = 1e-99
                global p_perp_b_pf = 1e-99
                γ_pt_pf = 1
                i_reason = 4
                keep_looping = false
                return keep_looping, helix_count, first_iter, i_grid_old
            end

            # Recalculate γ_pt_pf since that's needed elsewhere
            γ_pt_pf = hypot(pt_pf*o_o_mc, 1)

            # Modify components of pt_pf due to losses
            pb_pf       *= pt_pf/pt_pf_old
            p_perp_b_pf *= pt_pf/pt_pf_old

            # Also recalculate gyroradii
            global gyro_rad_tot_cm =   pt_pf * c_cgs * gyro_denom
            global gyro_rad_cm = p_perp_b_pf * c_cgs * gyro_denom
        end


        #---------------------------------------------
        # -- SCATTERING -- SCATTERING -- SCATTERING --
        #---------------------------------------------
        if !dont_scatter
            gyro_period_sec, pb_pf, p_perp_b_pf, φ_rad = scattering(
                aa, gyro_denom, pt_pf, γ_pt_pf, xn_per, pb_pf, p_perp_b_pf, φ_rad,
                use_custom_frg, p_electron_crit, γ_electron_crit, η_mfp,
            )
        end


        # Update acceleration time in explosion frame, so convert t_step from plasma frame
        # Only start the clock once a particle has crossed the shock for the first time
        if l_DwS
            acctime_sec += t_step*γ_uef

            if do_tcuts && acctime_sec ≥ tcuts[tcut_curr]
                tcut_track!(tcut_curr, xwt, pt_pf, wt_coupled, spectra_coupled, i_ion)
                tcut_curr += 1
            end

            # Remove ions at splitting momentum
            if pt_pf > pcuts_use[i_cut]
                @info "Removing ions at splitting momentum"
                l_save[i_prt] = true

                xwt_sav[i_prt]          = xwt
                pt_pf_sav[i_prt]        = pt_pf
                pb_pf_sav[i_prt]        = pb_pf
                x_PT_cm_sav[i_prt]      = r_PT_cm.x
                grid_sav[i_prt]         = i_grid
                DwS_sav[i_prt]          = l_DwS
                inj_sav[i_prt]          = inj
                xn_per_sav[i_prt]       = xn_per
                prp_x_cm_sav[i_prt]     = r_PT_cm.x < prp_x_cm ? prp_x_cm : r_PT_cm.x * 1.1 # Ensure particle is within PRP
                acctime_sec_sav[i_prt]  = acctime_sec
                φ_rad_sav[i_prt]        = φ_rad
                tcut_sav[i_prt]         = tcut_curr

                keep_looping = false
                return keep_looping, helix_count, first_iter, i_grid_old
            end
        end


        # Shift between coarse/fine xn_per; right now the code uses only one value of
        # xn_per everywhere in the ptot/x_PT phase space, so this branch does nothing.
        xn_per = r_PT_cm.x > gyro_rad_tot_cm ? xn_per_coarse : xn_per_fine
        #--------------------------------------------------------------
        # End of Code Block 3


    end  # check on i_return


    # Code Block 2: particle movement; fluxes; DwS escape/return
    #---------------------------------------------------------------
    # Store old x/y/z/φ positions
    global r_PT_old  = r_PT_cm
    global φ_rad_old = φ_rad


    # Find time step
    global t_step = gyro_period_sec / xn_per

    # Odd little loop that only matters if DSA has been disabled per the input file
    #------------------------------------------------------------------------------
    loop_again = true
    while loop_again # loop_no_DSA

        # Update position/phase angle.
        # Phase angle is easy: add appropriate amount and make sure result is in [0,2π)
        # For position, use inverse Lorentz transformations where primed frame is the
        # plasma frame, moving to the right along the x axis, and the unprimed frame is the
        # (stationary) shock frame. Take frames to be coincident
        # at t = t' = 0. Then
        #
        #    x  =  γ_usf * (x' + ux_sk*t'),
        #
        # where t' is t_step and x' is distance moved in plasma
        # frame (i.e. x_move_bpar below).
        # Remember to take gyration about magnetic field into account
        φ_rad = mod2pi(φ_rad  +  2π/xn_per)

        x_move_bpar = pb_pf * t_step / (γ_pt_pf * aa*mp_cgs)

        r_PT_cm = SVector(r_PT_old.x +  γ_usf * ( x_move_bpar * b_cosθ
                                       - gyro_rad_cm * b_sinθ * (cos(φ_rad)-cos(φ_rad_old))
                                       + ux_sk * t_step ),
                          r_PT_cm.y,
                          r_PT_cm.z)

        # Don't care about transverse motion unless it's being tracked
        #TODO: add new keyword for tracking transverse motion
        #if yz_pos_max > 0
        #    y_PT_cm = y_PT_old + gyro_rad_cm * ( sin(φ_rad) - sin(φ_rad_old) )
        #    z_PT_cm = z_PT_old + γ_usf * (x_move_bpar*b_sinθ - gyro_rad_cm*b_cosθ*(cos(φ_rad)-cos(φ_rad_old)) + uz_sk*t_step)
        #end


        #TODO: set flags for energy transfer here?
        #Search original code for "out of ions"


        # Reflect particles under two conditions:
        #  1. they've been DwS and are crossing back UpS
        #  2. either DSA is disabled or they fail an injection check
        # continue again if particles get reflected
        if r_PT_cm.x ≤ 0 && r_PT_old.x > 0 && !inj && (dont_DSA || inj_fracs[i_ion] < 1)
            rand = Random.rand()
            if dont_DSA || (rand > inj_fracs[i_ion])
                # Reflect particles with negative pb_pf; randomize phase of the rest
                if pb_pf < 0
                    pb_pf = -pb_pf
                else
                    rand = Random.rand()
                    φ_rad = rand * 2π
                end
            else
                loop_again = false
            end
        else
            loop_again = false
        end
    end # loop_no_DSA
    #----------------------------------------------------------------
    # DSA injection prevented if specified


    # If particle crosses shock going UpS -> DwS
    if r_PT_old.x < 0 && r_PT_cm.x ≥ 0

        l_DwS = true

        # Ensure the downstream region is sufficiently long to allow particle to isotropize.
        # This simple equation is decidedly non-trivial, and comes from two assumptions:
        # 1. the particle's diffusion coefficient D may be described by D = ⅓⋅η_mfp⋅r_g⋅v_pt
        #    (by default; a different f(r_g) may be specified as desired in place of η⋅r_g)
        # 2. the relation between the diffusion coefficient D and the diffusion length L is
        #    L = D/<u>, where <u> is the average speed of diffusion. Assuming isotropic
        #    particles in the DwS frame, <u> = u_2 since the average thermal *velocity* of
        #    the population is 0.
        #TODO: include f(r_g) in place of η*r_g to allow for arbitrary diffusion
        L_diff = η_mfp/3 * gyro_rad_tot_cm * pt_pf/(aa*mp_cgs*γ_pt_pf * u_2)

        @show prp_x_cm
        if prp_x_cm < L_diff
            prp_x_cm = L_diff
        end

    end


    # Test for injection into the acceleration process
    if l_DwS && r_PT_cm.x < 0
        inj = true
    end


    # Calculate fluxes due to this motion; also locates new grid zone
    (i_grid, i_grid_old, n_cr_count, px_esc_UpS, energy_esc_UpS) = all_flux!(
        i_prt, aa, pb_pf, p_perp_b_pf, pt_pf, γ_pt_pf, φ_rad,
        xwt, i_grid, ux_sk, uz_sk, utot, γ_usf, b_cosθ, b_sinθ,
        r_PT_cm.x, r_PT_old.x, inj, nc_unit,
        i_grid_feb, pxx_flux, pxz_flux,
        energy_flux, energy_esc_UpS, px_esc_UpS, spectra_sf, spectra_pf,
        n_cr_count, num_crossings, psd,
        n_xspec, x_spec, feb_UpS, oblique, γ_Z, u_Z, o_o_mc,
        n_grid, x_grid_cm,
        therm_grid, therm_px_sk, therm_pt_sk, therm_xwt,
    )


    # Downstream escape/return test; assume we'll call subroutine
    # prob_return unless told otherwise by particle location/info
    do_prob_ret = true

    if feb_DwS > 0 && r_PT_cm.x > feb_DwS

        # Particle flagged as escaping downstream
        i_return    = 0
        do_prob_ret = false

    elseif r_PT_cm.x > 1.1*prp_x_cm

        # The following simple equation is decidedly non-trivial, and comes from two
        # assumptions:
        # 1. the particle's diffusion coefficient D may be described by D = ⅓⋅η_mfp⋅r_g⋅v_pt
        #    (by default; a different f(r_g) may be specified as desired in place of η⋅r_g)
        # 2. the relation between the diffusion coefficient D and the diffusion length L is
        #    L = D/<u>, where <u> is the average speed of diffusion. Assuming isotropic
        #    particles in the DwS frame, <u> = u_2 since the average thermal *velocity* of
        #    the population is 0.
        #TODO: include f(r_g) in place of η*r_g to allow for arbitrary diffusion
        if aa < 1 && pt_pf < p_electron_crit
            gyro_fac = p_electron_crit * c_cgs * gyro_denom
            v_fac = gyro_fac * p_electron_crit / (aa*mp_cgs*γ_electron_crit * u_2)
        else
            v_fac = gyro_rad_tot_cm * pt_pf / (aa*mp_cgs*γ_pt_pf * u_2)
        end
        L_diff = η_mfp/3 * v_fac

        if r_PT_cm.x > 6.91 * L_diff
            i_return    = 0
            do_prob_ret = false
        end

    end

    if do_prob_ret
        (i_return, lose_pt, tcut_curr, _x_PT_cm, prp_x_cm, pt_pf, γ_pt_pf,
         gyro_denom, pb_pf, p_perp_b_pf, acctime_sec, φ_rad) = prob_return(
            rad_loss_fac, B_CMBz, r_PT_old.x, aa, zz, gyro_denom,
            r_PT_cm.x, prp_x_cm, pt_pf, γ_pt_pf, pb_pf, p_perp_b_pf,
            acctime_sec, φ_rad, helix_count, pcut_prev, xwt, tcut_curr,
            x_grid_stop, u_2, use_custom_εB, η_mfp, do_retro, bmag_2,
        )
        r_PT_cm = SVector(_x_PT_cm, r_PT_cm.y, r_PT_cm.z)
    end

    # If particle escaped DwS, handle final calculations here
    if i_return == 0

        vel = pt_pf / (aa*mp_cgs)
        if (γ_pt_pf - 1) ≥ energy_rel_pt
            vel /= γ_pt_pf
        end

        sum_pressure_DwS += pt_pf/3 * vel * xwt
        sum_KEdensity_DwS += (γ_pt_pf - 1) * aa*E₀_proton * xwt

        i_reason = 1
        if lose_pt
            i_reason = 4
        end

        keep_looping = false
        return keep_looping, helix_count, first_iter, i_grid_old
    end
    #----------------------------------------------------------------
    # End of Code Block 2
    return helix_count, keep_looping, first_iter, i_grid_old
end
