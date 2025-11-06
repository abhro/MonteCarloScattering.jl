# Floating point error can cause sin_Δφ to fall outside [-1,1]; place a limit on
# allowed values equal to the largest value of sin that can be distinguished from 1.0
const sin_upper_limit = prevfloat(1.0)

"""
    scattering(...)

This is a combination of two subroutines from the old code: `prob_scat` and `scattering`.
Randomly moves the particle's momentum vector along the surface of the unit sphere.

### Arguments

- `aa`: particle atomic mass
- `gyro_denom`: q*B, the denominator of gyroradius formula
- `ptot_pf`: total plasma frame momentum
- `γₚ_pf`: Lorentz factor associated with `ptot_pf`
- `xn_per`: number of time steps a gyroperiod is divided into

### Returns

- `gyro_period_sec`

### Modifies

- `pb_pf`: component of ptot_pf parallel to magnetic field
- `p_perp_b_pf`: component of ptot_pf perpendicular to magnetic field
- `φ_rad`: phase angle of gyration
"""
function scattering(
        rng, aa, gyro_denom, ptot_pf, γₚ_pf, xn_per,
        pb_pf, p_perp_b_pf, φ_rad,
        use_custom_frg, pₑ_crit, γₑ_crit, η_mfp)

    # If particle is an electron and p < pₑ_crit, use a constant MFP for scattering.
    # Note that instead addition to calculating the gyro period in seconds,
    # the code keeps vt_pf times the gyro period to find Δθ_max.
    if aa < 1 && ptot_pf < pₑ_crit
        gyro_rad_tot_cm =      pₑ_crit *         c * gyro_denom
        gyro_period_sec = 2π * γₑ_crit * aa*mp * c * gyro_denom
    else
        gyro_rad_tot_cm =    ptot_pf *         c * gyro_denom
        gyro_period_sec = 2π * γₚ_pf * aa*mp * c * gyro_denom
    end
    vp_tg = 2π * gyro_rad_tot_cm

    # To calculate collision time, need to know how MFP depends on gyro radius. Can either
    # use η_mfp⋅r_g (the default) or a user-specified custom f(r_g). Note that instead of
    # calculating the collision time in seconds, the code determines vt_pf times the
    # collision time, i.e. the mean free path
    if use_custom_frg
        error("Use of custom f(r_g) not yet supported. Add functionality or use standard.")
    else
        λ_mfp = η_mfp * gyro_rad_tot_cm
    end

    # Calculate the maximum allowed change in pitch angle; this formula is
    # slightly different from that used in previous version of the code
    cos_max = cos(√(6vp_tg / (xn_per*λ_mfp)))

    # Compute the actual change in pitch angle, as well as its modulation due to a
    # randomly-selected phase angle adjustment. See Ellison+ (1990) [1990ApJ...360..702E]
    # for more information.
    cos_old_pitch = pb_pf / ptot_pf
    sin_old_pitch = p_perp_b_pf / ptot_pf

    cos_Δθ = 1 - Random.rand(rng)*(1 - cos_max)  # Change of cos between [0, cos_max]
    sin_Δθ = √(1 - cos_Δθ^2)

    φ_scat = Random.rand(rng)*2π - π

    # Spherical law of cosines
    cos_new_pitch = cos_old_pitch*cos_Δθ + sin_old_pitch*sin_Δθ*cos(φ_scat)
    sin_new_pitch = √(1 - cos_new_pitch^2)

    # Adjust the components of ptot_pf and the phase angle
    pb_pf       = ptot_pf * cos_new_pitch
    p_perp_b_pf = ptot_pf * sin_new_pitch

    φ_p_old = φ_rad + π/2

    if sin_new_pitch != 0

        sin_Δφ = sin(φ_scat) * sin_Δθ / sin_new_pitch

        if abs(sin_Δφ) > sin_upper_limit
            sin_Δφ = copysign(sin_upper_limit, sin_Δφ)
        end

        φ_p_new = φ_p_old + asin(sin_Δφ)

    else
        # sin_new_pitch = 0, so no change at all
        φ_p_new = φ_p_old
    end

    φ_rad = φ_p_new - π/2

    return gyro_period_sec, pb_pf, p_perp_b_pf, φ_rad
end
