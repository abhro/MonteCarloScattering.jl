"""
    prob_return(...)

If the particle ends its movement downstream of the shock, perform a series
of tests to determine whether it will be culled from the simulation.

### Arguments

- `rad_loss_fac`: constant related to radiative losses; only used to pass to `retro_time` when needed
- `B_CMBz`: effective magnetic field due to CMB at redshift of source;
  only used to pass to function retro_time when needed
- `x_PT_old`: particle position before most recent move
- `aa`: atomic mass number of ion species
- `gyro_denom`: denominator of gyroradius fraction, q⋅B
- `helix_count`: counter for number of times through main propagation loop for current particle
- `i_cut`: current pcut; needed when electrons are undergoing radiative losses
- `pcut_prev`: momentum of previous pcut; needed when electrons are undergoing radiative losses
- `weight`: current particle weight; passed to retro_time if called

### Returns

- `i_return`: flag for fate of particle; see top of `loop_helix`
- `lose_pt`: true if particle hit zero energy due to radiative losses while in downstream region

### Modifies

- `x_PT_cm`: position of particle after recent motion & downstream adjustment
- `prp_x_cm`: location of probability-of-return plane
- `ptot_pf`: total plasma frame momentum of particle
- `γ`ₚ_pf: Lorentz factor associated with ptot_pf
- `gyro_denom`: denominator of gyroradius fraction, q⋅B
- `pb_pf`/`p_perp_b_pf`: components of `ptot_pf` parallel/perpendicular to B field
- `acctime_sec`: total accumulated acceleration time
- `φ_rad`: phase angle of particle's gyration
- `tcut_curr`: current tcut for particle tracking; passed to `retro_time` if called
"""
function prob_return(
        i_ion, num_psd_mom_bins, rad_loss_fac, B_CMBz, x_PT_old, aa, zz, gyro_denom,
        x_PT_cm, prp_x_cm, ptot_pf, γₚ_pf, pb_pf, p_perp_b_pf,
        acctime_sec, φ_rad, helix_count::Integer, pcut_prev, weight, tcut_curr,
        x_grid_stop, u₂, use_custom_εB, η_mfp, do_retro, B₂, mc,
        do_rad_losses, do_tcuts, tcuts,
        n_grid, uₓ_sk_grid, γ_sf_grid, γ_ef_grid, θ_grid, btot_grid,
        psd_mom_min, psd_bins_per_dec_mom,
        weight_coupled, spectra_coupled,
    )

    # Presume particle didn't enter probability of return calculation; change later as needed
    i_return = 2

    lose_pt = false

    # Test whether particle is still upstream from PRP position;
    # don't bother with return calculations if so
    if x_PT_cm < x_grid_stop

        # Do nothing

    # Particle has just crossed end of shock region as initially defined in input file
    elseif (x_PT_old < x_grid_stop) && (x_PT_cm ≥ x_grid_stop)

        # The following simple equation is decidedly non-trivial, and comes from
        # two assumptions:
        # 1. The particle's diffusion coefficient D may be described by D = ⅓⋅η_mfp⋅r_g⋅v_pt
        #    (by default; a different f(r_g) may be specified as desired in place of η⋅r_g)
        # 2. The relation between the diffusion coefficient D and the diffusion length L is:
        #    L = D/⟨u⟩, where ⟨u⟩ is the average speed of diffusion. Assuming
        #    isotropic particles in the downstream frame, ⟨u⟩ = u₂ since the
        #    average thermal *velocity* of the population is 0.
        # The calculation of gyro_tmp ensures that even particles that started upstream
        # of shock still use the downstream magnetic field for their diffusion length
        #TODO: include f(r_g) in place of η*r_g to allow for arbitrary diffusion
        # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
        if use_custom_εB && x_PT_cm > x_grid_stop
            gyro_tmp = √(x_grid_stop / x_PT_cm)
        else
            gyro_tmp = 1.0
        end
        gyro_rad_tot_cm = ptot_pf * c * gyro_tmp / (qcgs * B₂)
        L_diff          = η_mfp/3 * gyro_rad_tot_cm * ptot_pf/(aa*mp*γₚ_pf * u₂)

        # Make absolutely sure particles will have enough distance to isotropize before
        # encountering PRP; allow for three diffusion lengths beyond *current position*,
        # not just beyond end of grid
        prp_x_cm = x_PT_cm + 3L_diff


    # Particle has crossed PRP, and we need more complex calculations to determine if it returns
    elseif x_PT_old < prp_x_cm && x_PT_cm ≥ prp_x_cm
        vt_pf    = ptot_pf / (γₚ_pf * aa*mp)
        prob_ret = ((vt_pf - u₂) / (vt_pf + u₂))^2

        # If the particle's plasma frame velocity is less than u₂, or if the probability
        # of return calculation (see Jones & Ellison 1991 [1991SSRv...58..259J]) fails,
        # the particle will not return from the downstream region.
        if vt_pf < u₂ || Random.rand() > prob_ret

            i_return = 0

        # Particle will return from downstream region. Either analytically determine its properties
        # upon return or use retro-time calculation
        else

            i_return = 1

            # Track particle histories "explicitly" (see note in subroutine)
            if do_retro
                (lose_pt, φ_rad, tcut_curr, ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf,
                 gyro_denom, acctime_sec) = retro_time(
                    i_ion, num_psd_mom_bins, rad_loss_fac, B_CMBz, aa, zz, gyro_denom, prp_x_cm,
                    ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf, acctime_sec, weight,
                    tcut_curr,
                    use_custom_εB, x_grid_stop, do_rad_losses, do_tcuts, tcuts,
                    n_grid, uₓ_sk_grid, γ_sf_grid, γ_ef_grid, θ_grid, btot_grid,
                    mc, weight_coupled, spectra_coupled,
                    psd_mom_min, psd_bins_per_dec_mom,
                )

                # If electrons somehow lost all their energy due to radiative losses,
                # flag for removal
                if lose_pt
                    i_return = 0
                end

                # Particles return from retro_time at the location of the PRP
                x_PT_cm = prp_x_cm


            else # Analytically place particles back at PRP

                #TODO: check relations in Appendix A3 of Ellison, Baring & Jones [1996ApJ...473.1029E]
                # to verify that they are correct in relativistic limit.
                error("Code not set up for analytical PRP calculations. ",
                      "Must verify that EBJ1996 relations are correct in relativistic case")
            end

        end # check on prob_ret

    else

        # Particle is downstream from grid end, but didn't cross it this time step.
        # Also, it is upstream from the PRP. However, electrons experiencing
        # radiative losses may have a smaller L_diff during their propagation,
        # so a shorter PRP can be used. Test for that here. Two methods for doing that:
        # 1. If L_diff has dropped sufficiently, move the PRP far upstream from the particle's
        #    current position. The particle will be culled at the next time step
        # 2. Otherwise, calculate a new PRP location based on ratio of current momentum to
        #    minimum momentum for this pcut. The strong dependence on momentum (p⁵) is so
        #    that these electrons have time to isotropize downstream, even though
        #    the bulk of their motion occurred at a much higher energy and
        #    therefore mean free path
        if aa < 1 && ptot_pf < pcut_prev && helix_count % 1000 == 0
            gyro_rad_tot_cm = ptot_pf * c * gyro_denom
            L_diff          = η_mfp/3 * gyro_rad_tot_cm * ptot_pf/(aa*mp*γₚ_pf * u₂)

            if x_PT_cm > 2e3*L_diff
                prp_x_cm = 0.8 * x_PT_cm
            else
                prp_x_cm = min(prp_x_cm, x_grid_stop + L_diff*(pcut_prev/ptot_pf)^5)
            end
        end


    end # check on position vs x_grid_stop and prp_x_cm

    return (i_return, lose_pt, tcut_curr, x_PT_cm, prp_x_cm, ptot_pf, γₚ_pf,
            gyro_denom, pb_pf, p_perp_b_pf, acctime_sec, φ_rad)
end

# Floating point error can cause sin_Δφ to fall outside [-1,1]; place
# a limit on allowed values equal to the largest value of sin that can be
# distinguished from 1.0
##const sin_upper_limit = 1 - epsilon(1.0)

"""
    retro_time(...)

Explicitly tracks particles downstream of the probability of return plane.
Does so by the "retrodictive approach" explained in Jones (1978) [1978ApJ...222.1097J]
and also by Ellison, Jones, & Reynolds (1990) [1990ApJ...360..702E].

Since we know the particle will return to the PRP, we run time "backwards" by
placing the particle in a bulk flow opposite to its motion. Given a presumably
infinite downstream region of propagation, the particle is guaranteed to return to the
PRP over long enough time scales, and we assume that the "true" history of the
particle (its PRP1→downstream→PRP2 path) resembles on average its backwards history
(PRP1←downstream←PRP2).

### Arguments

- `rad_loss_fac`: constant related to radiative losses
- `B_CMBz`: effective magnetic field due to CMB at redshift of source
- `aa`: atomic mass number of ion species
- `zz`: charge of ion species
- `prp_x_cm`: starting (and ending) location for retro_time motion
- `weight`: current particle weight

### Returns

- `lose_pt`: whether particle is lost due to energy losses despite "returning" from a
  probabilistic standpoint
- `φ_rad`: phase angle of particle upon return

### Modifies

- `ptot_pf`: total plasma frame momentum of particle
- `pb_pf`/`p_perp_b_pf`: components of ptot_pf parallel/perpendicular to B field
- `γₚ_pf`: Lorentz factor associated with ptot_pf
- `gyro_denom`: denominator of gyroradius fraction, zz*B
- `acctime_sec`: total accumulated acceleration time
- `tcut_curr`: current tcut for particle tracking
"""
function retro_time(
        i_ion::Int, num_psd_mom_bins::Int, rad_loss_fac, B_CMBz, aa, zz, gyro_denom, prp_x_cm,
        ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf, acctime_sec, weight,
        tcut_curr,
        use_custom_εB, x_grid_stop, do_rad_losses, do_tcuts, tcuts,
        n_grid::Int, uₓ_sk_grid, γ_sf_grid, γ_ef_grid, θ_grid, btot_grid,
        mc, weight_coupled, spectra_coupled,
        psd_mom_min, psd_bins_per_dec_mom,
        # η_mfp,
    )

    # Set constants that will be used during the loop
    xn_per     = 10.0
    φ_step     = 2π / xn_per
    t_step_fac = 2π * aa*mp * c * gyro_denom / xn_per |> s  # t_step/γₚ_pf

    uₓ_sk = -uₓ_sk_grid[n_grid]
    γᵤ_sf =   γ_sf_grid[n_grid]
    γᵤ_ef =   γ_ef_grid[n_grid]
    B     =   btot_grid[n_grid]
    # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
    if use_custom_εB
        B *= √(x_grid_stop / prp_x_cm)
    end
    b_cosθ = cos(θ_grid[n_grid])
    b_sinθ = sin(θ_grid[n_grid])

    B_CMB_loc = B_CMBz * γᵤ_ef
    B²_tot    = B^2 + B_CMB_loc^2

    lose_pt = false

    # Initialize position, and phase angle
    x_PT = prp_x_cm

    φ_rad = Random.rand() * 2π


    # Main loop; note similarity to main loop of code, as it's doing most of
    # the same things with smaller scope
    while true

        # Store old values in preparation for the loop
        x_PT_old      = x_PT
        φ_rad_old     = φ_rad
        ptot_pf_old   = ptot_pf
        cos_old_pitch = pb_pf / ptot_pf
        sin_old_pitch = p_perp_b_pf / ptot_pf


        # Calculate true gyroradius and total gyroradius for particle;
        # if use_custom_εB is true, then also need to adjust magnetic field
        # for gyroradius and radiative cooling
        # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
        if use_custom_εB
            B          = btot_grid[n_grid] * √(x_grid_stop / x_PT)
            B²_tot     = B^2 + B_CMB_loc^2
            gyro_denom = 1 / (zz * B)
        end
        gyro_rad     = p_perp_b_pf * c * gyro_denom |> cm
        gyro_rad_tot =     ptot_pf * c * gyro_denom |> cm

        # Update φ_rad
        φ_rad = mod2pi(φ_rad_old + φ_step)

        # Calculate time step and movement distance; note that x_move_bpar =
        # pb_pf*t_step*m_pt/γₚ_pf, but t_step = t_step_fac*γₚ_pf, so the
        # factors of γₚ_pf divide out
        t_step = t_step_fac * γₚ_pf
        x_move_bpar = pb_pf * t_step_fac / (aa*mp) |> cm    # FIXME confirm formula


        # Move particle and update the acceleration time; note that we don't
        # care about y or z motion here, and that uₓ_sk is negative per the
        # definition above the loop
        x_PT = x_PT_old + γᵤ_sf * (x_move_bpar*b_cosθ - gyro_rad*b_sinθ*(cos(φ_rad)-cos(φ_rad_old)) + uₓ_sk*t_step)
        acctime_sec += t_step * γᵤ_ef

        # If tcut tracking is enabled, it should continue even during retro_time
        if do_tcuts && acctime_sec ≥ tcuts[tcut_curr]
            tcut_track!(weight_coupled, spectra_coupled, tcut_curr, weight, ptot_pf, i_ion, num_psd_mom_bins, psd_mom_min, psd_bins_per_dec_mom)
            tcut_curr += 1
        end

        # Large-angle scattering. Comment out pitch-angle diffusion section below if using LAS.
        φ_rad = 2π * Random.rand()   # Completely randomize phase angle

        pb_pf = (2Random.rand() - 1) * ptot_pf  # Completely randomize pb_pf
        p_perp_b_pf = √(ptot_pf^2 - pb_pf^2)

        # Pitch-angle diffusion. Comment out large-angle scattering section above if using PAD.
        #pitch_angle_diffusion()


        if do_rad_losses && aa < 1 # Radiative losses
            # Note that here dp_synch is actually dp/p. If this value is too large we will
            # directly integrate from p_i to get p_f, since the discrete approach would
            # result in too high a loss in a single time step
            dp_synch = rad_loss_fac * B²_tot * ptot_pf * t_step

            # Correction to make sure electrons don't lose too much energy in a single time step
            if dp_synch > 1e-2
                ptot_pf /= 1 + dp_synch
            else
                ptot_pf *= 1 - dp_synch # Put second factor of ptot_pf back into dp_synch
            end

        end  # check on radiative losses

        # Catch electrons that have somehow lost all their energy in a single time step,
        # and update the pitch angle of particles that remain
        if ptot_pf ≤ 0g*cm/s
            ptot_pf = 1e-99g*cm/s
            γₚ_pf = 1.0

            lose_pt = true
            break
        else
            pb_pf       = ptot_pf * cos_old_pitch
            p_perp_b_pf = ptot_pf * sin_old_pitch
            γₚ_pf = hypot(1, ptot_pf/mc)
        end


        # Check for return to PRP
        x_PT < prp_x_cm && break

    end  # retro time loop

    return lose_pt, φ_rad, tcut_curr, ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf, gyro_denom, acctime_sec
end

function pitch_angle_diffusion()

    # Compute maximum allowed pitch angle cosine
    vp_tg = 2π * gyro_rad_tot
    use_custom_frg && error("Use of custom f(r_g) not yet supported. Add functionality or use standard.")
    λ_mfp = η_mfp * gyro_rad_tot
    cos_max = cos(√(6vp_tg / (xn_per*λ_mfp)))

    # Compute change to pitch angle and roll
    Δcos = 1 - Random.rand()*(1 - cos_max)
    Δθ_scat = acos(Δcos)
    Δsin = sin(Δθ_scat)
    φ_scat = Random.rand()*2π - π

    # New pitch angle
    cos_new_pitch = cos_old_pitch * Δcos + sin_old_pitch * Δsin * cos(φ_scat)
    sin_new_pitch = √(1 - cos_new_pitch^2)

    # Adjust components of ptot_pf and phase angle
    pb_pf       = ptot_pf * cos_new_pitch
    p_perp_b_pf = ptot_pf * sin_new_pitch

    φ_p_old = φ_rad + π/2
    if sin_new_pitch != 0
        sin_Δφ = sin(φ_scat) * Δsin / sin_new_pitch
        if abs(sin_Δφ) > sin_upper_limit
            sin_Δφ = copysign(sin_upper_limit, sin_Δφ)
        end
        φ_p_new = φ_p_old + asin(sin_Δφ)
    else
        φ_p_new = φ_p_old
    end
    φ_rad = φ_p_new - π/2
end
