using .constants: mp_cgs, E₀_proton
using .parameters: na_grid, psd_max, energy_rel_pt, na_cr

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
- pt_pf
- γ_pt_pf
- φ_rad
- xwt
- ux_sk
- uz_sk
- utot
- γ_usf
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
- px_esc_UpS
- energy_esc_UpS
- pxx_flux (array modified in-place)
- pxz_flux (array modified in-place)
- energy_flux (array modified in-place)
- spectra_sf (array modified in-place)
- spectra_pf (array modified in-place)
- psd (array modified in-place)
"""
function all_flux!(
        i_prt, aa, pb_pf, p_perp_b_pf, pt_pf, γ_pt_pf, φ_rad,
        xwt, i_grid, ux_sk, uz_sk, utot, γ_usf, b_cosθ, b_sinθ,
        x_PT_cm, x_PT_old, inj, nc_unit,
        i_grid_feb, pxx_flux, pxz_flux, energy_flux, energy_esc_UpS, px_esc_UpS,
        spectra_sf, spectra_pf, n_cr_count, num_crossings, psd,
        # from controls
        n_xspec, x_spec, feb_UpS, oblique, γ_Z, u_Z, o_o_mc,
        # from grid_vars
        n_grid, x_grid_cm,

        # from species_vars
        # The therm_*** arrays can be pulled from the module because they are inherently
        # thread-safe. Only n_cr_count and num_crossings need protection from race
        # conditions, and so need to be included explicitly in the arguments of all_flux.
        therm_grid, therm_px_sk, therm_pt_sk, therm_xwt,
    )

    # Very early check to see if particle crossed a grid zone boundary
    i_grid_old = i_grid

    # Find the next boundary the particle hasn't crossed
    if x_PT_cm > x_PT_old
        i_grid = findfirst(>(x_PT_cm), x_grid_cm[i_grid+1:n_grid+1]) - 1
    else
        i_grid = findlast( ≤(x_PT_cm), x_grid_cm[i_grid:-1:1])
    end

    # Don't bother with *any* of the other computations if
    # 1. Grid zone hasn't changed, and
    # 2. There's no possibility of crossing an intra-grid detector
    ##TODO: remove extra return statement, folding rest of computation into
    # "else" block of if-then
    if i_grid == i_grid_old && i_grid > i_grid_feb && iszero(n_xspec)
        return (i_grid, i_grid_old, n_cr_count, px_esc_UpS, energy_esc_UpS)
    end

    # Convert plasma frame momentum to shock frame; determine a few values
    # that will be reused during the call to all_flux
    pt_sk, px_sk, pz_sk, γ_pt_sk = transform_p_PS(
        aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad, ux_sk, uz_sk, utot, γ_usf,
        b_cosθ, b_sinθ, oblique, o_o_mc)

    if pt_sk > abs(px_sk*all_flux_spike_away)
        pt_o_px_sk = all_flux_spike_away
        # Minimum shock frame velocity is a small fraction of local bulk flow speed
        abs_inv_vx_sk = abs(all_flux_spike_away/ux_sk)
    else
        pt_o_px_sk = pt_sk / px_sk
        abs_inv_vx_sk = abs(γ_pt_sk * aa*mp_cgs / px_sk)
    end

    if pt_pf > abs(pb_pf*all_flux_spike_away)
        pt_o_px_pf = all_flux_spike_away
    else
        pt_o_px_pf = abs(pt_pf / pb_pf)
    end

    # Kinetic energy only; rest mass energy NOT included
    if (γ_pt_sk - 1) > energy_rel_pt
        energy_flux_add = (γ_pt_sk - 1) * aa*E₀_proton * xwt
    else
        energy_flux_add = pt_sk^2 / (2 * aa*mp_cgs) * xwt
    end


    # Calculate spectrum at x_spec locations if needed
    if n_xspec > 0
        i_pt    = get_psd_bin_momentum(pt_sk, psd_bins_per_dec_mom,
                                       # from some module
                                       psd_mom_min, num_psd_mom_bins)
        i_pt_pf = get_psd_bin_momentum(pt_pf, psd_bins_per_dec_mom,
                                       # from some module
                                       psd_mom_min, num_psd_mom_bins)

        for i in 1:n_xspec
            if (( x_PT_old < x_spec[i] && x_PT_cm  ≥ x_spec[i] ) ||
                ( x_PT_cm  ≤ x_spec[i] && x_PT_old > x_spec[i] ))

                # Spectrum in shock frame
                spectra_sf[i_pt, i] += xwt * pt_o_px_sk


                #TODO: if the shock is oblique, need to replace pb_pf with px_pf below,
                # as well as changing the outputs from transform_p_PS to include px_pf at all
                oblique && error("Cannot calculate flux_wt_fac when pb_pf != px_pf")
                # Spectrum in plasma frame; flux_wt_fac corresponds to vx_pf/vx_sk and
                # measures relative likelihood of crossing in the plasma frame
                # given a known crossing in the shock frame
                flux_wt_fac = abs(pb_pf/px_sk) * (γ_pt_sk/γ_pt_pf)

                spectra_pf[i_pt_pf, i] += xwt * pt_o_px_pf * flux_wt_fac
            end
        end
    end


    # Check if particle has already been injected into acceleration process;
    # if it has get the PSD bins it will be placed into
    if inj
        i_pt = get_psd_bin_momentum(pt_sk, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
        local jθ = get_psd_bin_angle(px_sk, pt_sk, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)
    end


    # Main loop for all_flux: depending on motion of particle, travel UpS or
    # DwS and update fluxes across all zone boundaries the particle crossed
    # WARNING: the particle quantity "xwt" is the fraction of the far UpS
    # density each particle represents. However, the actual flux is
    #     γ_Z * u_Z * density_Z,
    # which means that the flux contribution of each particle must be
    # increased by a factor of γ_Z*u_Z.
    #--------------------------------------------------------------------------
    # Downstream first
    if x_PT_cm > x_PT_old
        flux_downstream!(pxx_flux, pxz_flux, energy_flux, i_grid_old, i_grid)
    else # Particle has moved upstream
        flux_upstream!(pxx_flux, pxz_flux, energy_flux, i_grid_old, i_grid)
    end  # check on direction of motion
    #--------------------------------------------------------------------------
    # Finished updating fluxes for crossed grid boundaries


    # One final task: update tracker for escaping flux at FEB if needed;
    # don't forget that flux must be rescaled
    #$omp critical
    if inj && x_PT_cm < feb_UpS && x_PT_old ≥ feb_UpS
        energy_esc_UpS += energy_flux_add * γ_Z*u_Z
        px_esc_UpS     -= px_sk * xwt      * γ_Z*u_Z
    end
    #$omp end critical

    return (i_grid, i_grid_old, n_cr_count, px_esc_UpS, energy_esc_UpS)
end

function flux_downstream!(pxx_flux, pxz_flux, energy_flux, i_grid_old, i_grid)
    # Loop over grid zones boundaries that the particle has crossed
    for i in i_grid_old+1:i_grid

        # Update fluxes
        pxx_flux[i]    +=      px_sk *xwt * γ_Z*u_Z
        pxz_flux[i]    +=  abs(pz_sk)*xwt * γ_Z*u_Z
        energy_flux[i] += energy_flux_add * γ_Z*u_Z

        # Update PSD, or add to list of tracked thermal particles
        if inj
            psd[i_pt, jθ, i] += xwt * abs_inv_vx_sk
        else
            # Particle is a thermal particle. Increment crossing counter and attempt to add
            # entry to the various arrays. If arrays already full, write out crossing data to
            # scratch file for tracking them. Then update array holding number of crossings.
            #$omp critical
            if n_cr_count < na_cr
                n_cr_count += 1
                therm_grid[n_cr_count] = i
                therm_px_sk[n_cr_count] = px_sk
                therm_pt_sk[n_cr_count] = pt_sk
                therm_xwt[n_cr_count]  = xwt * abs_inv_vx_sk
            else
                # Need to write to scratch file, formatted or otherwise
                println(nc_unit, i, i_prt, px_sk, pt_sk, xwt * abs_inv_vx_sk)
            end  # check on n_cr_count
            #$omp end critical

            num_crossings[i] += 1

        end  # check on inj

    end  # loop over grid zones
end

function flux_upstream!(pxx_flux, pxz_flux, energy_flux, i_grid_old, i_grid)

    # Loop over grid zones boundaries that the particle has crossed
    for i in i_grid_old:-1:i_grid+1

        # Is particle UpS of free escape boundary after crossing DwS?
        # Then don't count flux contributions to grid zones UpS from FEB
        inj && i ≤ i_grid_feb && continue

        # Update fluxes; note the minus signs force pxx_flux to increase (b/c particle
        # is moving UpS and px_sk < 0) and force energy_flux to decrease
        pxx_flux[i]    -=      px_sk *xwt * γ_Z*u_Z
        pxz_flux[i]    +=  abs(pz_sk)*xwt * γ_Z*u_Z
        energy_flux[i] -= energy_flux_add * γ_Z*u_Z

        # Update PSD, or add to list of tracked thermal particles
        if inj
            psd[i_pt, jθ, i] += xwt * abs_inv_vx_sk
        else

            # Particle is a thermal particle. Increment crossing counter and attempt to add
            # entry to the various arrays. If arrays already full, write out crossing data to
            # scratch file for tracking them. Then update array holding number of crossings.
            #$omp critical
            if n_cr_count < na_cr
                n_cr_count += 1
                therm_grid[n_cr_count] = i
                therm_px_sk[n_cr_count] = px_sk
                therm_pt_sk[n_cr_count] = pt_sk
                therm_xwt[n_cr_count]  = xwt * abs_inv_vx_sk
            else
                # Need to write to scratch file, formatted or otherwise
                println(nc_unit, i, i_prt, px_sk, pt_sk, xwt * abs_inv_vx_sk)
            end  # check on n_cr_count
            #$omp end critical

            num_crossings[i] += 1

        end  # check on inj
    end  # loop over grid zones
end
