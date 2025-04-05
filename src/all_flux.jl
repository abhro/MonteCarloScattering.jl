using .constants: E₀_proton
using .parameters: energy_rel_pt, na_cr

all_flux_spike_away = 1000.0 # Max value for 1/cosine

"""
Tracks particle flux due to motion on the grid. Also finds number of new grid zone

### Arguments
TODO
- i_prt
- nc_unit
- i_grid_feb
- inj
- aa
- pb_pf
- p_perp_b_pf
- ptot_pf
- γₚ_pf
- φ_rad
- weight
- uₓ_sk
- uz_sk
- utot
- γᵤ_sf
- b_cosθ
- b_sinθ
- x_PT_cm
- x_PT_old

### Modifies
TODO
- i_grid
- i_grid_old
- n_cr_count
- num_crossings (array modified in-place)
- pₓ_esc_UpS
- energy_esc_UpS
- pₓₓ_flux (array modified in-place)
- pxz_flux (array modified in-place)
- energy_flux (array modified in-place)
- spectra_sf (array modified in-place)
- spectra_pf (array modified in-place)
- psd (array modified in-place)
"""
function all_flux!(
        i_prt, aa, pb_pf, p_perp_b_pf, ptot_pf, γₚ_pf, φ_rad,
        weight, i_grid, uₓ_sk, uz_sk, utot, γᵤ_sf, b_cosθ, b_sinθ,
        x_PT_cm, x_PT_old, inj, nc_unit,
        i_grid_feb, pₓₓ_flux, pxz_flux, energy_flux, energy_esc_UpS, pₓ_esc_UpS,
        spectra_sf, spectra_pf, n_cr_count, num_crossings, psd,
        # from controls
        n_xspec, x_spec, feb_UpS, γ₀, u₀, mc,
        # from grid_vars
        n_grid, x_grid_cm,

        # from species_vars
        # The therm_*** arrays can be pulled from the module because they are inherently
        # thread-safe. Only n_cr_count and num_crossings need protection from race
        # conditions, and so need to be included explicitly in the arguments of all_flux.
        therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight,
    )

    # Very early check to see if particle crossed a grid zone boundary
    i_grid_old = i_grid

    # Find the next boundary the particle hasn't crossed
    if x_PT_cm > x_PT_old
        i_grid = findnext(>(x_PT_cm), x_grid_cm, i_grid+1) - 1
    else
        i_grid = findprev(≤(x_PT_cm), x_grid_cm, i_grid)
    end
    if isnothing(i_grid)
        error("Failed to find crossing in ", x_grid_cm, " with x_PT_cm = ", x_PT_cm)
    end

    # Don't bother with *any* of the other computations if
    # 1. Grid zone hasn't changed, and
    # 2. There's no possibility of crossing an intra-grid detector
    ##TODO: remove extra return statement, folding rest of computation into
    # "else" block of if-then
    if i_grid == i_grid_old && i_grid > i_grid_feb && iszero(n_xspec)
        return (i_grid, i_grid_old, n_cr_count, pₓ_esc_UpS, energy_esc_UpS)
    end

    # Convert plasma frame momentum to shock frame; determine a few values
    # that will be reused during the call to all_flux
    ptot_sk, p_sk, γₚ_sk = transform_p_PS(
        aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad, uₓ_sk, uz_sk, utot, γᵤ_sf,
        b_cosθ, b_sinθ, mc)

    if ptot_sk > abs(p_sk.x*all_flux_spike_away)
        pt_o_pₓ_sk = all_flux_spike_away
        # Minimum shock frame velocity is a small fraction of local bulk flow speed
        abs_inv_vx_sk = abs(all_flux_spike_away/uₓ_sk)
    else
        pt_o_pₓ_sk = ptot_sk / p_sk.x
        abs_inv_vx_sk = abs(γₚ_sk * aa*mp / p_sk.x)
    end

    pt_o_pₓ_pf = min(abs(ptot_pf/pb_pf), all_flux_spike_away)

    # Kinetic energy only; rest mass energy NOT included
    if (γₚ_sk - 1) > energy_rel_pt
        energy_flux_add = (γₚ_sk - 1) * aa*E₀_proton * weight
    else
        energy_flux_add = ptot_sk^2 / (2 * aa*mp) * weight
    end


    # Calculate spectrum at x_spec locations if needed
    if n_xspec > 0
        i_pt    = get_psd_bin_momentum(ptot_sk, psd_bins_per_dec_mom,
                                       # from some module
                                       psd_mom_min, num_psd_mom_bins)
        i_pt_pf = get_psd_bin_momentum(ptot_pf, psd_bins_per_dec_mom,
                                       # from some module
                                       psd_mom_min, num_psd_mom_bins)

        for i in 1:n_xspec
            if ((x_PT_old < x_spec[i] && x_PT_cm  ≥ x_spec[i]) ||
                (x_PT_cm  ≤ x_spec[i] && x_PT_old > x_spec[i]))

                # Spectrum in shock frame
                spectra_sf[i_pt, i] += weight * pt_o_pₓ_sk


                # Spectrum in plasma frame; flux_weight_fac corresponds to vx_pf/vx_sk and
                # measures relative likelihood of crossing in the plasma frame
                # given a known crossing in the shock frame
                flux_weight_fac = abs(pb_pf/p_sk.x) * (γₚ_sk/γₚ_pf)

                spectra_pf[i_pt_pf, i] += weight * pt_o_pₓ_pf * flux_weight_fac
            end
        end
    end



    # Main loop for all_flux: depending on motion of particle, travel UpS or DwS and update
    # fluxes across all zone boundaries the particle crossed. WARNING: the particle quantity
    # "weight" is the fraction of the far UpS density each particle represents. However, the
    # actual flux is
    #     γ₀ * u₀ * n₀,
    # which means that the flux contribution of each particle must be increased by a factor
    # of γ₀⋅u₀.
    #--------------------------------------------------------------------------
    # Downstream first
    if x_PT_cm > x_PT_old
        i_range = i_grid_old+1:i_grid
        inj_check = false
        sign_fac = 1
    else # Particle has moved upstream
        i_range = i_grid_old:-1:i_grid+1
        inj_check = true
        sign_fac = -1
    end  # check on direction of motion
    n_cr_count = flux_stream!(
        i_range, inj_check, sign_fac,
        pₓₓ_flux, pxz_flux, energy_flux,
        p_sk, ptot_sk, weight, energy_flux_add, inj, n_cr_count, abs_inv_vx_sk)
    #--------------------------------------------------------------------------
    # Finished updating fluxes for crossed grid boundaries


    # One final task: update tracker for escaping flux at FEB if needed;
    # don't forget that flux must be rescaled
    #$omp critical
    if inj && x_PT_cm < feb_UpS && x_PT_old ≥ feb_UpS
        energy_esc_UpS += energy_flux_add * γ₀*u₀
        pₓ_esc_UpS     -= p_sk.x * weight * γ₀*u₀
    end
    #$omp end critical

    return (i_grid, i_grid_old, n_cr_count, pₓ_esc_UpS, energy_esc_UpS)
end

function flux_stream!(
        i_range, inj_check, sign_fac,
        pₓₓ_flux, pxz_flux, energy_flux,
        p_sk, ptot_sk, weight, energy_flux_add, inj, n_cr_count, abs_inv_vx_sk,
    )

    # Check if particle has already been injected into acceleration process;
    # if it has get the PSD bins it will be placed into
    if inj
        i_pt = get_psd_bin_momentum(ptot_sk, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
        jθ = get_psd_bin_angle(p_sk.x, ptot_sk, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)
    end

    # Loop over grid zones boundaries that the particle has crossed
    for i in i_range

        # Is particle UpS of free escape boundary after crossing DwS?
        # Then don't count flux contributions to grid zones UpS from FEB
        inj_check && inj && i ≤ i_grid_feb && continue

        # Update fluxes; note the minus signs in sign_fac force pₓₓ_flux to increase
        # (b/c particle is moving UpS and p_sk.x < 0) and force energy_flux to decrease
        pₓₓ_flux[i]    += sign_fac *  p_sk.x *weight * γ₀*u₀
        pxz_flux[i]    +=         abs(p_sk.z)*weight * γ₀*u₀
        energy_flux[i] += sign_fac * energy_flux_add * γ₀*u₀

        # Update PSD, or add to list of tracked thermal particles
        if inj
            psd[i_pt, jθ, i] += weight * abs_inv_vx_sk
        else
            # Particle is a thermal particle. Increment crossing counter and attempt to add
            # entry to the various arrays. If arrays already full, write out crossing data to
            # scratch file for tracking them. Then update array holding number of crossings.
            #$omp critical
            if n_cr_count < na_cr
                n_cr_count += 1
                therm_grid[n_cr_count] = i
                therm_pₓ_sk[n_cr_count] = p_sk.x
                therm_ptot_sk[n_cr_count] = ptot_sk
                therm_weight[n_cr_count]  = weight * abs_inv_vx_sk
            else
                # Need to write to scratch file, formatted or otherwise
                println(nc_unit, i, i_prt, p_sk.x, ptot_sk, weight * abs_inv_vx_sk)
            end  # check on n_cr_count
            #$omp end critical

            num_crossings[i] += 1

        end  # check on inj
    end  # loop over grid zones
    return n_cr_count
end
