"""
Reads in the PSD -- as well as the thermal arrays & scratch file -- to
calculate the pressure (which may be anisotropic) everywhere on the grid.

### Arguments

- num_crossings: array containing number of thermal particle crossings recorded at each grid zone
- n_cr_count: total number (up to na_cr) of thermal particle crossings; signals whether
  scratch file was used
- therm_grid: array of grid zone number for thermal crossings
- therm_px_sk: array of x-momentum for thermal crossings
- therm_pt_sk: array of total momentum for thermal crossings
- therm_xwt: array of flux-adjusted particle weights for thermal particle crossings
- nc_unit: unit number for the scratch file holding crossing data
- psd: 3-D phase space distribution holding shock-frame ptot and cos(θ) for all recorded CR
  grid crossings
- zone_pop: (Lorentz-invariant) number of particles in each grid zone

### Returns

- pressure_psd_par: component of pressure parallel to magnetic field
- pressure_psd_perp: component of pressure perpendicular to magnetic field
- energy_density_psd: local kinetic energy density of fluid in every grid zone
"""
function thermo_calcs(
        num_crossings, n_cr_count, therm_grid, therm_px_sk,
        therm_pt_sk, therm_xwt, nc_unit, psd, zone_pop,
        #pressure_psd_par, pressure_psd_perp, energy_density_psd,
        aa_ion, zz_ion, tZ_ion, denZ_ion, psd_lin_cos_bins,
        γ_Z, β_Z, tZ_electron, sc_electron
    )

    # Output arguments
    pressure_psd_par = zeros(na_grid)
    pressure_psd_perp = zeros(na_grid)
    energy_density_psd = zeros(na_grid)

    d2N_pf = fill(1e-99, (0:psd_max, 0:psd_max, na_grid))

    # Initialize d2N_pf to "zero"
    # Dimensions of arrays (chosen for maximum speed during loops):
    #     1 - cos(θ)
    #     2 - ptot
    #     3 - grid position
    # The omission of "dpdcos" is not a typo; the array will only ever hold a
    # total count of particles, not a true dN/dp
    # Set constants to be used repeatedly, including calculating the center
    # points of all bins to save time later
    rest_mass_energy = aa_ion[i_ion] * E₀_proton  # Rest mass-energy of the current particle species

    cos_center = zeros(0:psd_max)
    for jθ in 0:num_psd_θ_bins
        # Determine current cosines, remembering that psd_θ_bounds has both a
        # linearly-spaced region in cosine and logarithmically-spaced region in θ.
        if jθ > (num_psd_θ_bins - psd_lin_cos_bins)
            cos_hi = psd_θ_bounds[jθ]
            cos_lo = psd_θ_bounds[jθ + 1]
        elseif jθ == (num_psd_θ_bins - psd_lin_cos_bins)
            cos_hi = cos(psd_θ_bounds[jθ])
            cos_lo = psd_θ_bounds[jθ + 1]
        else
            cos_hi = cos(psd_θ_bounds[jθ])
            cos_lo = cos(psd_θ_bounds[jθ + 1])
        end

        # Minus sign needed because finest gradations must point UpS
        cos_center[jθ] = - ( cos_lo + cos_hi ) / 2
    end

    pt_cgs_center = zeros(0:psd_max)
    for k in 0:num_psd_mom_bins
        # Convert from log to linear space
        pt_cgs_center[k] = exp10((psd_mom_bounds[k] + psd_mom_bounds[k+1]) / 2)
    end


    # Read the scratch files associated with thermal particles. Bin them into
    # the combined d2N/(dp-dcos) for the shock frame
    #----------------------------------------------------------------------------
    # Set up the arrays that will hold the crossing data
    max_cross = maximum(num_crossings,1)
    therm_px = zeros(max_cross, 0:n_grid+1)
    therm_pt = zeros(max_cross, 0:n_grid+1)
    therm_wt = zeros(max_cross, 0:n_grid+1)


    # Fill d2N_pf with data from thermal particles
    #----------------------------------------------------------------------------
    ntot_crossings = sum(num_crossings,1)
    # n_cross_fill needs to be an array because we will be skipping around the grid
    # as we move through the data rather than filling a single zone at a time
    n_cross_fill = zeros(Int, na_grid)

    rewind(nc_unit)


    # Handle crossings stored within the crossing arrays
    for i in 1:n_cr_count
        i_grid = therm_grid[i]

        n_cross_fill[i_grid] += 1

        # Note: ordering of coordinates chosen to make memory accesses in next loop faster
        therm_px[n_cross_fill[i_grid], i_grid] = therm_px_sk[i]
        therm_pt[n_cross_fill[i_grid], i_grid] = therm_pt_sk[i]
        therm_wt[n_cross_fill[i_grid], i_grid] = therm_xwt[i]
    end

    if ntot_crossings > na_cr

        # Need to go into the scratch file for the remainder of the crossings
        for i in 1:ntot_crossings-na_cr
            read(nc_unit, i_grid, idum, px_sk, pt_sk, cell_wt)

            n_cross_fill[i_grid] += 1

            # Note: ordering of coordinates chosen to make memory accesses in next loop faster
            therm_px[n_cross_fill[i_grid], i_grid] = px_sk
            therm_pt[n_cross_fill[i_grid], i_grid] = pt_sk
            therm_wt[n_cross_fill[i_grid], i_grid] = cell_wt
        end

    end
    #-------------------------------------------------------------------------
    # Arrays filled


    # Loop over grid locations, filling in the array d2N_pf
    #----------------------------------------------------------------------------
    for i in 1:na_grid
        num_crossings[i] == 0 && continue # Ignore zones that no thermal particles crossed

        # Get Lorentz factor and speed relating the shock and plasma frames
        γ_u  = γ_sf_grid[i]
        β_u = ux_sk_grid[i] / c_cgs


        # Loop over all crossings for this grid zone, binning the thermal particles.
        for n in 1:num_crossings[i]

            # Transform the center of the zone into the new frame
            pt_sk    = therm_pt[n,i]
            px_sk    = therm_px[n,i]
            etot_sk = hypot( pt_sk*c_cgs, rest_mass_energy )

            px_Xf    = γ_u * (px_sk - β_u*etot_sk/c_cgs)
            pt_Xf    = √( (pt_sk^2 - px_sk^2) + px_Xf^2 )
            # In rare cases, floating point roundoff can cause pt_Xf to be smaller than abs(px_Xf)
            if abs(px_Xf) > pt_Xf
                px_Xf = copysign(pt_Xf, px_Xf)
            end

            # Get location of center in transformed d2N_pf
            kpt_Xf = get_psd_bin_momentum(pt_Xf, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
            jθ_Xf = get_psd_bin_angle(px_Xf, pt_Xf, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)


            # And add particle count to the appropriate bin in d2N_pf
            d2N_pf[jθ_Xf, kpt_Xf, i] += therm_wt[n,i]
        end

    end
    #-------------------------------------------------------------------------
    # Plasma frame d2N filled from thermal arrays


    # Deallocate the crossing data arrays
    #deallocate(therm_px)
    #deallocate(therm_pt)
    #deallocate(therm_wt)
    #----------------------------------------------------------------------------
    # Thermal particles binned


    # Now, get the CRs into d2N_pf. Do this by transforming the bin *centers* of
    # PSD into each frame, using just the center location to rebin.
    # TODO: consider doing an exact calculation of bin overlaps rather than
    # just the center-point rebinning that currently happens
    #----------------------------------------------------------------------------
    d2N_pop = zeros(0:psd_max)
    # Loop over grid locations, rebinning all cells with particles
    for i in 1:n_grid

        # Get Lorentz factor and speed relating the shock and current frames
        γ_u = γ_sf_grid[i]
        β_u = ux_sk_grid[i] / c_cgs


        # Loop over cells in d2N_sf
        for jθ in 0:num_psd_mom_bins
            for k in 0:num_psd_θ_bins

                psd[k,jθ,i] ≤ 1e-66 && continue # Skip empty zones

                cell_wt = psd[k,jθ,i]

                # Transform the center of the zone into the new frame
                cos_θ_sk = cos_center[jθ]
                pt_sk     = pt_cgs_center[k]
                px_sk     = pt_sk * cos_θ_sk
                etot_sk  = hypot( pt_sk*c_cgs,  rest_mass_energy )

                px_Xf    = γ_u * (px_sk - β_u*etot_sk/c_cgs)
                pt_Xf    = √( pt_sk^2 - px_sk^2 + px_Xf^2 )

                # Get location of center in transformed d2N_pf
                kpt_Xf = get_psd_bin_momentum(pt_Xf, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
                jθ_Xf = get_psd_bin_angle(px_Xf, pt_Xf, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)

                # And add particle count to the appropriate bin in d2N_pf
                d2N_pf[jθ_Xf, kpt_Xf, i] += cell_wt

            end # loop over angles
        end # loop over ptot


        norm_fac = sum( d2N_pf[:,:,i], MASK = (d2N_pf[:,:,i]>1e-66) )
        if iszero(num_crossings[i]) && norm_fac > 0
            norm_fac += denZ_ion[i_ion] / ux_sk_grid[i]
        end
        if norm_fac > 0
            norm_fac = zone_pop[i] / norm_fac
        end


        for k in 0:psd_max-1, jθ in 0:psd_max-1
            if d2N_pf[jθ,k,i] > 1e-66
                d2N_pf[jθ,k,i] *= norm_fac
            end
        end

        d2N_pop[i] = sum( d2N_pf[:,:,i], MASK = (d2N_pf[:,:,i]>1e-66) )
    end # loop over grid locations
    #----------------------------------------------------------------------------
    # Transformed d2N_pf calculated


    # Now sum over d2N_pf to get the pressure at all grid locations.
    # With angular information we can get both parallel and perpendicular components
    #----------------------------------------------------------------------------
    # Find the velocity associated with each momentum bin
    vel_ptot = @. pt_cgs_center * c_cgs * o_o_mc /  √( 1  +  (pt_cgs_center*o_o_mc)^2 )


    # Loop over grid zones and find pressure components everywhere
    #----------------------------------------------------------------------------
    for i in 1:n_grid

        density_loc = γ_Z*β_Z*denZ_ion[i_ion] / √( γ_sf_grid[i]^2 - 1 )

        # Determine the normalization factors to use during the pressure calculation
        # Three cases: (1) no detected particles, (2) only CRs detected, and
        # (3) thermal particles detected
        #-----------------------------------------------------------------------
        if maximum( d2N_pf[:,:,i] ) < 1e-66 && iszero(num_crossings[i])

            # Case (1): No particles of any kind detected at this grid location.
            # Thermal particles must have passed through, so find their density
            # and analytically determine components of pressure
            # #assumecold: using γ_adiab = 5/3 in the pressure calc
            pressure_loc = density_loc^(5//3) * kB_cgs*tZ_ion[i_ion]

            if aa_ion[i_ion] ≥ 1
                density_electron = density_loc * zz_ion[i_ion]
            end

            # Add pressure due to electrons if they aren't a separate species
            if ! sc_electron
                pressure_loc += density_electron * kB_cgs*tZ_electron
            end

            # Update components of pressure, assuming isotropy; twice as many
            # perpendicular components as parallel (y & z vs x), thus the factor
            # of 2 in the computation
            pressure_psd_par[i]  += 1/3 * pressure_loc
            pressure_psd_perp[i] += 2/3 * pressure_loc

            # With only thermal particles at this grid location, the energy density is easy to determine
            # #assumecold
            energy_density_psd[i] += 1.5 * pressure_loc

            # Don't bother running through the d2N calculation
            continue


        elseif num_crossings[i] == 0

            # Case (2): no thermal particles, but some CRs detected.
            # Find pressure due to un-tracked thermal particles
            # #assumecold: using γ_adiab = 5/3 in the pressure calc
            pressure_loc = density_loc^(5//3) * kB_cgs*tZ_ion[i_ion]

            if aa_ion[i_ion] > 1
                density_electron *= zz_ion[i_ion]
            end

            # Add pressure due to electrons if they aren't a separate species
            if ! sc_electron
                pressure_loc += density_electron * kB_cgs*tZ_electron
            end

            # Scale the contributions from the thermal particles by the fraction of the
            # zone's population they represent, and add that to the running total
            pressure_loc *= 1 - d2N_pop[i]/zone_pop[i]
            pressure_psd_par[i]  += 1/3 * pressure_loc
            pressure_psd_perp[i] += 2/3 * pressure_loc

            # Set the normalization factor for the CRs that were tracked during
            # the iteration; summation over CRs will total d2N_pop[i], which
            # makes up for the amount subtracted from pressure_loc above
            norm_fac = density_loc / zone_pop[i]

            # Finally, calculate contribution of thermal particles to energy density.
            # This will be updated during the following loop to include CR contribution.
            # #assumecold
            energy_density_psd[i] += 3//2 * pressure_loc

        else

            # Case (3): thermal particles detected. If no CRs are present that
            # just means they didn't propagate UpS to this zone. Whether CRs
            # are present or not, d2N_pf represents a true count of the total
            # population of this grid zone. Set normalization factor accordingly
            norm_fac = density_loc / zone_pop[i]

        end
        #-----------------------------------------------------------------------
        # Appropriate normalization constant(s) found


        # Find the pressure using the density calculated above
        for k in 0:num_psd_mom_bins

            # These factors are the same in all cells of this ptot column
            pressure_fac  = 1/3 * pt_cgs_center[k] * vel_ptot[k] * norm_fac
            γ_tmp         = √( 1  +  (pt_cgs_center[k]*o_o_mc)^2 )
            energy_density_fac = (γ_tmp - 1) * rest_mass_energy

            for jθ in 0:num_psd_θ_bins

                # Don't bother with empty cells
                d2N_pf[jθ,k,i] < 1e-66 && continue

                pressure_par_tmp  = d2N_pf[jθ,k,i] * pressure_fac * cos_center[jθ]^2
                # The factor of two below is because there are two perpendicular components of pressure
                pressure_perp_tmp = d2N_pf[jθ,k,i] * pressure_fac * (1 - cos_center[jθ]^2)

                pressure_psd_par[i]  += pressure_par_tmp
                pressure_psd_perp[i] += pressure_perp_tmp


                # Update the energy density also
                energy_density_psd[i] += energy_density_fac * d2N_pf[jθ,k,i] * norm_fac
            end
        end

    end  # loop over i
    #-------------------------------------------------------------------------
    # Loop over grid zones finished

    #deallocate(d2N_pf)
    return pressure_psd_par, pressure_psd_perp, energy_density_psd
end
