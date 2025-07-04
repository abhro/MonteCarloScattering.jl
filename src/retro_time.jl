using .parameters: na_c

# Floating point error can cause sin_Δφ to fall outside [-1,1]; place
# a limit on allowed values equal to the largest value of sin that can be
# distinguished from 1.0
#const sin_upper_limit = 1 - epsilon(1.0)

"""
Explicitly tracks particles downstream of the probability of return plane.
Does so by the "retrodictive approach" explained in Jones (1978) [1978ApJ...222.1097J]
and also by Ellison, Jones, & Reynolds (1990) [1990ApJ...360..702E].

Since we know the particle will return to the PRP, we run time "backwards" by
placing the particle in a bulk flow opposite to its motion. Given a presumably
infinite DwS region of propagation, the particle is guaranteed to return to the
PRP over long enough time scales, and we assume that the "true" history of the
particle (its PRP1→DwS→PRP2 path) resembles on average its backwards history
(PRP1←DwS←PRP2).

### Arguments

- rad_loss_fac: constant related to radiative losses
- B_CMBz: effective magnetic field due to CMB at redshift of source
- aa: atomic mass number of ion species
- zz: charge number of ion species
- prp_x_cm: starting (and ending) location for retro_time motion
- weight: current particle weight

### Returns

- lose_pt: whether particle is lost due to energy losses despite "returning" from a
  probabilistic standpoint
- φ_rad: phase angle of particle upon return

### Modifies

- ptot_pf: total plasma frame momentum of particle
- pb_pf/p_perp_b_pf: components of ptot_pf parallel/perpendicular to B field
- γₚ_pf: Lorentz factor associated with ptot_pf
- gyro_denom: denominator of gyroradius fraction, zz*qcgs*bmag
- acctime_sec: total accumulated acceleration time
- tcut_curr: current tcut for particle tracking
"""
function retro_time(
        rad_loss_fac, B_CMBz, aa, zz, gyro_denom, prp_x_cm,
        ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf, acctime_sec, weight,
        tcut_curr,
        use_custom_εB, x_grid_stop, do_rad_losses, do_tcuts, tcuts,
        n_grid, uₓ_sk_grid, γ_sf_grid, γ_ef_grid, θ_grid, btot_grid,
        mc,
        # η_mfp,
    )

    # Set constants that will be used during the loop
    xn_per     = 10.0
    φ_step     = 2π / xn_per
    t_step_fac = 2π * aa*mp * c * gyro_denom / xn_per  # t_step/γₚ_pf

    uₓ_sk = -uₓ_sk_grid[n_grid]
    γᵤ_sf =   γ_sf_grid[n_grid]
    γᵤ_ef =   γ_ef_grid[n_grid]
    bmag  =   btot_grid[n_grid]
    # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
    if use_custom_εB
        bmag *= √(x_grid_stop / prp_x_cm)
    end
    b_cosθ = cos(θ_grid[n_grid])
    b_sinθ = sin(θ_grid[n_grid])

    B_CMB_loc = B_CMBz * γᵤ_ef
    B_tot_sq  = bmag^2 + B_CMB_loc^2

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
        ptot_pf_old     = ptot_pf
        cos_old_pitch = pb_pf / ptot_pf
        sin_old_pitch = p_perp_b_pf / ptot_pf


        # Calculate true gyroradius and total gyroradius for particle;
        # if use_custom_εB is true, then also need to adjust magnetic field
        # for gyroradius and radiative cooling
        # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
        if use_custom_εB
            bmag       = btot_grid[n_grid] * √(x_grid_stop / x_PT)
            B_tot_sq   = bmag^2 + B_CMB_loc^2
            gyro_denom = 1 / (zz * bmag)
        end
        gyro_rad_cm     = p_perp_b_pf * c * gyro_denom |> cm
        gyro_rad_tot_cm =     ptot_pf * c * gyro_denom |> cm

        # Update φ_rad
        φ_rad = mod2pi(φ_rad_old + φ_step)

        # Calculate time step and movement distance; note that x_move_bpar =
        # pb_pf*t_step*m_pt/γₚ_pf, but t_step = t_step_fac*γₚ_pf, so the
        # factors of γₚ_pf divide out
        t_step      =         t_step_fac * γₚ_pf
        @debug "" pb_pf t_step_fac aa mp gyro_denom xn_per
        x_move_bpar = pb_pf * t_step_fac * aa*mp |> cm


        # Move particle and update the acceleration time; note that we don't
        # care about y or z motion here, and that uₓ_sk is negative per the
        # definition above the loop
        x_PT = x_PT_old + γᵤ_sf * (x_move_bpar*b_cosθ - gyro_rad_cm*b_sinθ*(cos(φ_rad)-cos(φ_rad_old)) + uₓ_sk*t_step)
        acctime_sec += t_step * γᵤ_ef

        # If tcut tracking is enabled, it should continue even during retro_time
        if do_tcuts && acctime_sec ≥ tcuts[tcut_curr]
            tcut_track(tcut_curr, weight, ptot_pf)
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
            dp_synch = rad_loss_fac * B_tot_sq * ptot_pf * t_step

            # Correction to make sure electrons don't lose too much energy in a single time step
            if dp_synch > 1e-2
                ptot_pf /= 1 + dp_synch
            else
                ptot_pf *= 1 - dp_synch # Put second factor of ptot_pf back into dp_synch
            end

        end  # check on radiative losses

        # Catch electrons that have somehow lost all their energy in a single time step,
        # and update the pitch angle of particles that remain
        if ptot_pf ≤ 0
            ptot_pf   = 1e-99
            γₚ_pf = 1.0

            lose_pt = true
            break
        else
            pb_pf       = ptot_pf * cos_old_pitch
            p_perp_b_pf = ptot_pf * sin_old_pitch
            γₚ_pf     = hypot(1, ptot_pf/mc)
        end


        # Check for return to PRP
        x_PT < prp_x_cm && break

    end  # retro time loop

    return lose_pt, φ_rad, tcut_curr, ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf, gyro_denom, acctime_sec
end

function pitch_angle_diffusion()

    # Compute maximum allowed pitch angle cosine
    vp_tg = 2π * gyro_rad_tot_cm
    use_custom_frg && error("Use of custom f(r_g) not yet supported. Add functionality or use standard.")
    λ_mfp = η_mfp * gyro_rad_tot_cm
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
