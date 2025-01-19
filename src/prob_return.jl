using .constants: mₚ_cgs, qₚ_cgs, c_cgs
include("retro_time.jl")

"""
If the particle ends its movement downstream of the shock, perform a series
of tests to determine whether it will be culled from the simulation.

### Arguments

- rad_loss_fac: constant related to radiative losses; only used to pass to retro_time when needed
- B_CMBz: effective magnetic field due to CMB at redshift of source; only used to pass to
  function retro_time when needed
- x_PT_old: particle position before most recent move
- aa: atomic mass number of ion species
- gyro_denom: denominator of gyroradius fraction, zz*qₚ_cgs*bmag
- helix_count: counter for number of times through main propagation loop for current particle
- i_cut: current pcut; needed when electrons are undergoing radiative losses
- pcut_prev: momentum of previous pcut; needed when electrons are undergoing radiative losses
- weight: current particle weight; passed to retro_time if called

### Returns

- i_return: flag for fate of particle; see top of loop_helix
- lose_pt: true if particle hit zero energy due to radiative losses while in DwS region

### Modifies

- x_PT_cm: position of particle after recent motion & DwS adjustment
- prp_x_cm: location of probability-of-return plane
- ptot_pf: total plasma frame momentum of particle
- γₚ_pf: Lorentz factor associated with ptot_pf
- gyro_denom: denominator of gyroradius fraction, zz*qₚ_cgs*bmag
- pb_pf/p_perp_b_pf: components of ptot_pf parallel/perpendicular to B field
- acctime_sec: total accumulated acceleration time
- φ_rad: phase angle of particle's gyration
- tcut_curr: current tcut for particle tracking; passed to retro_time if called
"""
function prob_return(
        rad_loss_fac, B_CMBz, x_PT_old, aa, zz, gyro_denom,
        x_PT_cm, prp_x_cm, ptot_pf, γₚ_pf, pb_pf, p_perp_b_pf,
        acctime_sec, φ_rad, helix_count::Integer, pcut_prev, weight, tcut_curr,
        x_grid_stop, u₂, use_custom_εB, η_mfp, do_retro, bmag₂, mc)

    # Presume particle didn't enter probability of return calculation; change later as needed
    i_return = 2

    lose_pt = false

    # Test whether particle is still UpS from PRP position; don't bother with
    # return calculations if so
    if x_PT_cm < x_grid_stop

        # Do nothing

    # Particle has just crossed end of shock region as initially defined in input file
    elseif (x_PT_old < x_grid_stop) && (x_PT_cm ≥ x_grid_stop)

        # The following simple equation is decidedly non-trivial, and comes from
        # two assumptions:
        # 1. The particle's diffusion coefficient D may be described by D = ⅓⋅η_mfp⋅r_g⋅v_pt
        #    (by default; a different f(r_g) may be specified as desired in place of η⋅r_g)
        # 2. The relation between the diffusion coefficient D and the diffusion length L is:
        #    L = D/<u>, where <u> is the average speed of diffusion. Assuming isotropic
        #    particles in the DwS frame, <u> = u₂ since the average thermal *velocity* of
        #    the population is 0.
        # The calculation of gyro_tmp ensures that even particles that started UpS of shock
        # still use the DwS magnetic field for their diffusion length
        #TODO: include f(r_g) in place of η*r_g to allow for arbitrary diffusion
        # Square root corresponds to Blandford-McKee solution, where e ∝ 1/χ ∝ 1/r
        if use_custom_εB && x_PT_cm > x_grid_stop
            gyro_tmp = √(x_grid_stop / x_PT_cm)
        else
            gyro_tmp = 1.0
        end
        gyro_rad_tot_cm = ptot_pf * c_cgs * gyro_tmp / (qₚ_cgs * bmag₂)
        L_diff          = η_mfp/3 * gyro_rad_tot_cm * ptot_pf/(aa*mₚ_cgs*γₚ_pf * u₂)

        # Make absolutely sure particles will have enough distance to isotropize before
        # encountering PRP; allow for three diffusion lengths beyond *current position*,
        # not just beyond end of grid
        prp_x_cm = x_PT_cm + 3L_diff


    # Particle has crossed PRP, and we need more complex calculations to determine if it returns
    elseif x_PT_old < prp_x_cm && x_PT_cm ≥ prp_x_cm
        vt_pf    = ptot_pf / (γₚ_pf * aa*mₚ_cgs)
        prob_ret = ((vt_pf - u₂) / (vt_pf + u₂))^2

        # If the particle's plasma frame velocity is less than u₂, or if the probability
        # of return calculation (see Jones & Ellison 1991 [1991SSRv...58..259J]) fails,
        # the particle will not return from the DwS region.
        if vt_pf < u₂ || Random.rand() > prob_ret

            i_return = 0

        # Particle will return from DwS region. Either analytically determine its properties
        # upon return or use retro-time calculation
        else

            i_return = 1

            # Track particle histories "explicitly" (see note in subroutine)
            if do_retro
                (lose_pt, φ_rad, tcut_curr, ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf,
                 gyro_denom, acctime_sec) = retro_time(
                    rad_loss_fac, B_CMBz, aa, zz, gyro_denom, prp_x_cm,
                    ptot_pf, pb_pf, p_perp_b_pf, γₚ_pf, acctime_sec, weight,
                    tcut_curr,
                    use_custom_εB, x_grid_stop, do_rad_losses, do_tcuts, tcuts,
                    n_grid, uₓ_sk_grid, γ_sf_grid, γ_ef_grid, θ_grid, btot_grid,
                    mc,
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

        # Particle is DwS from grid end, but didn't cross it this time step. Also, it is UpS
        # from the PRP. However, electrons experiencing radiative losses may have a smaller
        # L_diff during their propagation, so a shorter PRP can be used. Test for that here.
        # Two methods for doing that:
        # 1. If L_diff has dropped sufficiently, move the PRP far UpS from the particle's
        #    current position. The particle will be culled at the next time step
        # 2. Otherwise, calculate a new PRP location based on ratio of current momentum to
        #    minimum mometum for this pcut. The strong dependence on momentum (p⁵) is so
        #    that these electrons have time to isotropize DwS, even though the bulk of their
        #    motion occurred at a much higher energy and therefore mean free path
        if aa < 1 && ptot_pf < pcut_prev && helix_count % 1000 == 0
            gyro_rad_tot_cm = ptot_pf * c_cgs * gyro_denom
            L_diff          = η_mfp/3 * gyro_rad_tot_cm * ptot_pf/(aa*mₚ_cgs*γₚ_pf * u₂)

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
