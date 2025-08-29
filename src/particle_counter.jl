module particle_counter
using Unitful, UnitfulAstro
using Unitful: cm, c, mp
using UnitfulAstro: kpc, pc
using OffsetArrays: OffsetVector

using ..constants: E₀ₚ
using ..parameters: psd_max, na_cr, num_therm_bins
using ..transformers: get_transform_dN, transform_psd_corners
#using ..debug: zone_vol, therm_energy_density, energy_density
import ..print_plot_vals

"""
    get_dNdp_cr(...)

Calculate dN/dp (NOT normalized) in plasma frame for particles that *have* been
injected into acceleration process.

### Arguments
FIXME

None, but does use PSD information from module psd_vars

### Returns
- `dNdp_cr`: 3-D array, containing 1-D array for each grid zone, of dN/dp for
  particles injected into acceleration process, for three different reference frames
"""
function get_dNdp_cr(
        m_ion, psd_lin_cos_bins, γ₀, num_psd_θ_bins, psd_θ_bounds,
        psd, psd_mom_bounds, num_psd_mom_bins, n_grid, γ_sf_grid, i_ion)

    # Zero out the output array before summing in loops to follow
    # Note on dNdp_cr's third dimension:
    #     1  -  Shock frame
    #     2  -  Plasma frame
    #     3  -  ISM frame
    dNdp_cr = zeros(0:psd_max, n_grid, 3)

    ct_bounds = fill(-2.0, 0:num_psd_θ_bins+1)

    # Rest mass of particle species, used here in binning pitch angle
    rest_mass = m_ion[i_ion]


    # Set values needed to bin pitch angle
    for l in 0:num_psd_θ_bins+1
        # Determine current cosine, remembering that psd_θ_bounds has both a linearly-spaced
        # region in cosine and logarithmically-spaced region in θ. Also need to remember
        # that the most finely spaced bins should occur in the upstream-pointing direction,
        # so need to negate psd_θ_bounds to get true cosine value.
        if l > (num_psd_θ_bins - psd_lin_cos_bins)
            ct_bounds[l] = -psd_θ_bounds[l]
        else
            ct_bounds[l] = -cos(psd_θ_bounds[l])
        end
    end


    # Degree of approximation to use in distributing cell_weight if NOT in shock
    # frame (where uniform distribution may be assumed):
    #   i_approx = 0:  assume uniform distribution
    #   i_approx = 1:  assume isosceles triangular distribution of cell_weight
    #   i_approx = 2:  assume scalene triangular distribution of cell_weight, with
    #                  peak of triangle centered on mean of ct_hi_pt and ct_lo_pt
    #   i_approx = 3:  exact calculation of fractional area, i.e. use no approximations
    i_approx = 2

    # Loop over grid zones and cos(θ)/ptot space to build dN_cr
    for k in 1:n_grid

        # Handle shock frame separately, as no special treatment is needed to
        # convert from PSD to dN(p). Conversion to dN/dp will happen at the
        # end of this subroutine.
        #----------------------------------------------------------------------
        for j in 0:psd_max, i in 0:psd_max
            if psd[i,j,k] > 0
                dNdp_cr[i,k,1] += psd[i,j,k]
            end
        end
        #-------------------------------------------------------------------------
        # Shock frame dN(p) found


        # Now handle plasma and ISM frames. Do this sequentially, since the
        # process is exactly the same but requires a single input that differs
        # between the two frames.
        for m in 2:3

            # Transform corners of PSD into correct frame for use in finding
            # p_cell_lo and p_cell_hi in next block of code
            #----------------------------------------------------------------------
            if m == 2       #  Plasma frame
                γᵤ = γ_sf_grid[k]
            elseif m == 3   #  ISM frame
                γᵤ = γ₀
            else
                error("ERROR with frame selection in get_dNdp_cr")
            end

            transform_corner_pt, transform_corner_ct = transform_psd_corners(γᵤ)
            #----------------------------------------------------------------------
            # corners transformed


            # Determine minimum and maximum plasma-frame momenta particles can have.
            # Allocate the cos(θ) arrays based on those values, since particles
            # will be binned according to both cos(θ) and (coarse-grained) momentum.
            # Then zero them out.
            i_ct_ptot_sk_min = floor(Int, psd_mom_bounds[1])
            i_ct_ptot_sk_max =  ceil(Int, psd_mom_bounds[num_psd_mom_bins+1])
            cθ_sk_xw = zeros(0:psd_max, i_ct_ptot_sk_min:i_ct_ptot_sk_max)

            i_ct_ptot_pf_min = floor(Int, minimum(transform_corner_pt[begin+1:end, :]))
            i_ct_ptot_pf_max =  ceil(Int, maximum(transform_corner_pt[begin+1:end, :]))
            cθ_pf_xw = zeros(0:psd_max, i_ct_ptot_pf_min:i_ct_ptot_pf_max)
            cθ_ef_xw = zeros(0:psd_max, i_ct_ptot_pf_min:i_ct_ptot_pf_max)

            # Note that dNdp_cr as calculated here is dN(p), *not* dN/dp.
            # That conversion happens at the end of this subroutine
            dNdp_cr[:,k,m] = get_transform_dN(psd[:,:,k], m, transform_corner_pt, transform_corner_ct, γᵤ, i_approx)

        end # loop over frames


        # If tracking pitch angles, and any bins were filled while doing so, plot histogram
        # of their values. Because pitch angle was broken down by momentum, loop over all
        # columns that aren't empty in ct_**_xw.
        # WARNING: this must be re-checked and re-tested since it was brought
        # into new version of code
        #------------------------------------------------------------------------
        ## First, histograms of individual momentum decades in the shock frame
        #for i in i_ct_ptot_sk_min:i_ct_ptot_sk_max
        #    for l in 0:num_psd_θ_bins
        #        if cθ_sk_xw[l,i] > 0
        #            cθ_sk_xw[l,i] /= (ct_bounds[l] - ct_bounds[l+1])
        #        end
        #    end
        #    ∑cθ_sk_xw = sum(cθ_sk_xw[:,i])
        #    iszero(∑cθ_sk_xw) && continue   # skip empty rows
        #
        #    k_unit = k + 8000 + 100*(numion-1)
        #    j_plot = 0
        #    for l in 0:num_psd_θ_bins
        #        j_plot += 1
        #
        #        θ = ct_bounds[l] == 1 ? 1e-99 : acos(ct_bounds[l])
        #        θ_next = acos(ct_bounds[l+1])
        #
        #        write(k_unit, k, j_plot,               # fort.80**
        #                      ct_bounds[l],            # 1
        #                      θ,                       # 2
        #                log10(θ),                      # 3
        #                      i,                       # 4
        #                      cθ_sk_xw[l,i]/∑cθ_sk_xw) # 5
        #
        #        j_plot = j_plot + 1
        #        if l < num_psd_θ_bins
        #            write(k_unit, k, j_plot,               # fort.80**
        #                          ct_bounds[l+1],          # 1
        #                          θ_next,                  # 2
        #                    log10(θ_next),                 # 3
        #                          i,                       # 4
        #                          cθ_sk_xw[l,i]/∑cθ_sk_xw) # 5
        #        end
        #    end
        #
        #    print_plot_vals(k_unit)
        #end

        ## Next, histogram for *all* momenta in the shock frame
        #∑cθ_sk_xw = sum(cθ_sk_xw)
        #if ∑cθ_sk_xw > 0  # skip grid locations w/no CRs
        #    k_unit = k + 8000 + 100*(numion-1)
        #    j_plot = 0
        #    for l in 0:num_psd_θ_bins
        #        j_plot += 1
        #
        #        θ = ct_bounds[l] == 1 ? 1e-99 : acos(ct_bounds[l])
        #        θ_next = acos(ct_bounds[l+1])
        #
        #        write(k_unit,"(2i4,1p48e15.7)") k, j_plot,  # fort.80**
        #                      ct_bounds[l],                 # 1
        #                      θ,                            # 2
        #                log10(θ),                           # 3
        #                  sum(cθ_sk_xw[l,:])/∑cθ_sk_xw   # 4
        #
        #        j_plot = j_plot + 1
        #        if l < num_psd_θ_bins
        #            write(k_unit,"(2i4,1p48e15.7)") k, j_plot,  # fort.80**
        #                          ct_bounds[l+1],               # 1
        #                          θ_next,                       # 2
        #                    log10(θ_next),                      # 3
        #                      sum(cθ_sk_xw[l,:])/∑cθ_sk_xw   # 4
        #        end
        #    end
        #    print_plot_vals(k_unit)
        #end


        ## Now histograms of individual momentum decades, in plasma frame
        #for i in i_ct_ptot_pf_min:i_ct_ptot_pf_max
        #
        #    for l in 0:num_psd_θ_bins
        #        if cθ_pf_xw[l,i] > 0
        #            cθ_pf_xw[l,i] /= (ct_bounds[l] - ct_bounds[l+1])
        #        end
        #    end
        #    ∑cθ_pf_xw = sum(cθ_pf_xw[:,i])
        #    ∑cθ_pf_xw == 0 && continue   # skip empty rows
        #
        #    k_unit = k + 9000 + 100*(numion-1); j_plot = 0
        #    for l in 0:num_psd_θ_bins
        #        j_plot += 1
        #
        #        θ = ct_bounds[l] == 1 ? 1e-99 : acos(ct_bounds[l])
        #        θ_next = acos(ct_bounds[l+1])
        #
        #        write(k_unit,      # fort.90**
        #              k, j_plot,
        #              ct_bounds[l],                  # 1
        #              θ,                             # 2
        #              log10(θ),                      # 3
        #              i,                             # 4
        #              cθ_pf_xw[l,i]/∑cθ_pf_xw)       # 5
        #
        #        j_plot += 1
        #        if l < num_psd_θ_bins
        #            write(k_unit,    # fort.90**
        #                  k, j_plot,
        #                  ct_bounds[l+1],            # 1
        #                  θ_next,                    # 2
        #                  log10(θ_next),             # 3
        #                  i,                         # 4
        #                  cθ_pf_xw[l,i]/∑cθ_pf_xw)   # 5
        #        end
        #
        #    end
        #
        #    print_plot_vals(k_unit)
        #end


        ## Finally, histogram for *all* momenta in the plasma frame
        #∑cθ_pf_xw = sum(cθ_pf_xw)
        #if ∑cθ_pf_xw > 0  # skip grid locations w/no CRs
        #    k_unit = k + 9000 + 100*(numion-1)
        #    j_plot = 0
        #    for l in 0:num_psd_θ_bins
        #        j_plot += 1
        #
        #        θ = ct_bounds[l] == 1 ? 1e-99 : acos(ct_bounds[l])
        #        θ_next = acos(ct_bounds[l+1])
        #
        #        write(k_unit,      # fort.90**
        #              k, j_plot,
        #              ct_bounds[l],                  # 1
        #              θ,                             # 2
        #              log10(θ),                      # 3
        #              sum(cθ_pf_xw[l,:])/∑cθ_pf_xw)  # 4
        #
        #        j_plot += 1
        #        if l < num_psd_θ_bins
        #            write(k_unit,
        #                  k, j_plot,    # fort.90**
        #                  ct_bounds[l+1],                # 1
        #                  θ_next,                        # 2
        #                  log10(θ_next),                 # 3
        #                  sum(cθ_pf_xw[l,:])/∑cθ_pf_xw)  # 4
        #        end
        #
        #    end
        #
        #    print_plot_vals(k_unit)
        #end
        #------------------------------------------------------------------------
        # pitch angle histograms plotted

    end # loop over grid locations


    # Finally, convert dNdp_cr into true dN/dp by dividing by dp for each cell
    for m in 1:3, k in 1:n_grid, l in 0:psd_max-1

        # Skip empty cells in PSD
        if dNdp_cr[l,k,m] < 1e-66
            dNdp_cr[l,k,m] = 1e-99
            continue
        end

        dNdp_cr[l, k, m] /= exp10(psd_mom_bounds[l+1]) - exp10(psd_mom_bounds[l])
    end

    return dNdp_cr
end

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

"""
    get_dNdp_2D(...)

Generates d²N/(dp dcos) (or something like it) for electrons in the explosion frame.

First, use the scratch file and the phase space distribution to generate d²N/(dp dcos),
the 2-D version of dN/dp (with extra information about the angular distribution). Combine
the results from thermal and CR particles into a single d²N/(dp dcos) to save memory.
After d²N/(dp dcos) is found in the shock frame, transform it to the plasma and ISM frames
by re-binning the *center* of each shock frame bin. Leave the computationally difficult
problem of overlap for a later date. Finally, condense each d²N/(dp dcos) into dN/dp for
comparison against the results of the original dN/dp subroutines.

!!! warning

    d²N/(dp dcos) isn't technically accurate. To get number, just multiply by dp
    and sum over the cos(θ) column.

!!! note

    Right now, this is only used to calculate the ISM-frame distribution of electrons;
    all other frame/particle pairs can be done with the 1-D dN/dp subroutines.

### Arguments

- `nc_unit`: unit number for the scratch file holding crossing data
- `zone_pop`: (Lorentz-invariant) number of particles in each grid zone

### Returns

- `d²N_dpdcos_ef`: explosion-frame array holding dN/dp spread out across angular dimension
"""
function get_dNdp_2D(
        nc_unit, zone_pop, m_ion, n_ions, n₀_ion,
        psd_lin_cos_bins, γ₀, β₀, psd, num_psd_θ_bins,
        psd_θ_bounds, num_psd_mom_bins, psd_mom_bounds, n_grid,
        γ_sf_grid, i_ion, num_crossings, n_cr_count, therm_grid,
        therm_pₓ_sk, therm_ptot_sk, therm_weight,
        psd_bins_per_dec_mom, psd_mom_min, psd_bins_per_dec_θ, psd_cos_fine,
        Δcos, psd_θ_min,
    )

    # Initialize d²N_dpdcos_sf to "zero"
    # Dimensions of d²N_dpdcos arrays (chosen for maximum speed during loops):
    #     1 - cos(θ)
    #     2 - ptot
    #     3 - grid position
    # The omission of "dpdcos" in the plasma frame array is not a typo; the
    # array will only ever hold a total count of particles, not a true dN/dp
    d²N_dpdcos_sf = fill(1e-99, (0:psd_max, 0:psd_max, n_grid))
    d²N_pf        = fill(1e-99, (0:psd_max, 0:psd_max, n_grid))
    d²N_dpdcos_ef = fill(1e-99, (0:psd_max, 0:psd_max, n_grid))


    # Set constants to be used repeatedly
    E₀ = m_ion[i_ion] * c^2 # Rest mass-energy of the current particle species
    Δp = zeros(0:num_psd_mom_bins)
    for k in eachindex(Δp)
        Δp[k] = exp10(psd_mom_bounds[k+1]) - exp10(psd_mom_bounds[k])
    end


    # Read the scratch files associated with thermal particles. Bin them into
    # the combined d²N/(dp-dcos) for the shock frame
    #----------------------------------------------------------------------------
    # Set up the arrays that will hold the crossing data
    max_cross = maximum(num_crossings)
    therm_pₓ = zeros(max_cross, 0:n_grid+1)
    therm_pt = zeros(max_cross, 0:n_grid+1)
    therm_weight = zeros(max_cross, 0:n_grid+1)


    # Fill the arrays with data from the scratch file
    #----------------------------------------------------------------------------
    ntot_crossings = sum(num_crossings)
    # n_cross_fill needs to be an array because we will be skipping around the
    # grid as we move through the data rather than filling a single zone at a time
    n_cross_fill = zeros(Int, n_grid)

    ##rewind(nc_unit) # XXX Fortran holdover

    # Handle crossings stored within the crossing arrays
    for i in 1:n_cr_count
        i_grid = therm_grid[i]
        n_cross_fill[i_grid] += 1

        # Note: ordering of coordinates chosen to make memory accesses in next loop faster
        therm_pₓ[n_cross_fill[i_grid], i_grid] = therm_pₓ_sk[i]
        therm_pt[n_cross_fill[i_grid], i_grid] = therm_ptot_sk[i]
        therm_weight[n_cross_fill[i_grid], i_grid] = therm_weight[i]
    end

    if ntot_crossings > na_cr

        # Need to go into the scratch file for the remainder of the crossings
        for i in 1:ntot_crossings-na_cr
            i_grid, _, pₓ_sk, ptot_sk, cell_weight = read(nc_unit)

            n_cross_fill[i_grid] += 1

            # Note: ordering of coordinates chosen to make memory accesses in next loop faster
            therm_pₓ[n_cross_fill[i_grid], i_grid] = pₓ_sk
            therm_pt[n_cross_fill[i_grid], i_grid] = ptot_sk
            therm_weight[n_cross_fill[i_grid], i_grid] = cell_weight
        end
    end

    # With the read-in complete, the scratch file can be closed
    close(nc_unit)
    #-------------------------------------------------------------------------
    # Arrays filled


    # Loop over grid locations, filling in the array d²N_dpdcos_sf
    #-------------------------------------------------------------------------
    for i in 1:n_grid
        if num_crossings[i] == 0
            d²N_dpdcos_sf[:,:,i] .= 1e-99
            continue     # Ignore zones that no thermal particles crossed
        end

        # Loop over all crossings for this grid zone, binning the thermal particles.
        for n in 1:num_crossings[i]
            k = get_psd_bin_momentum(therm_pt[n,i],
                                     psd_bins_per_dec_mom,
                                     psd_mom_min, num_psd_mom_bins)
            jθ = get_psd_bin_angle(therm_pₓ[n,i], therm_pt[n,i],
                                   psd_bins_per_dec_θ,
                                   num_psd_mom_bins, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)
            d²N_dpdcos_sf[jθ,k,i] += therm_weight[n,i]
        end
    end
    #-------------------------------------------------------------------------
    # Shock frame d²N_dpdcos filled from thermal arrays


    # Deallocate the crossing data arrays
    #deallocate(therm_pₓ)
    #deallocate(therm_pt)
    #deallocate(therm_weight)
    #----------------------------------------------------------------------------
    # Thermal particles binned


    # With the thermal particles taken care of in the shock frame, move on to the cosmic ray
    # particles. Note that the slices of psd need to be transposed from (ptot,θ) to (θ,ptot)
    # order. Once that is taken care of, convert dN into dN/dp by dividing by dp
    #-------------------------------------------------------------------------
    for i in 1:n_grid+1, k in 0:psd_max-1, jθ in 0:psd_max-1
        if psd[k,jθ,i] > 1e-66
            d²N_dpdcos_sf[jθ,k,i] = psd[k,jθ,i]
        end
    end
    #-------------------------------------------------------------------------
    # Cosmic rays finished


    # Calculate the density of particles represented by d²N_dpdcos_sf
    density_tot = zeros(n_grid)
    for i in 1:n_grid
        mask = (d²N_dpdcos_sf[:,0:psd_max-1,i] .> 1e-66)
        density_tot[i] = sum(d²N_dpdcos_sf[mask...,i])
    end


    # Can convert from dN to dN/dp now
    for k in 0:psd_max-1, jθ in 0:psd_max
        if d²N_dpdcos_sf[jθ,k,i] > 1e-66
            d²N_dpdcos_sf[jθ,k,i] /= Δp[k]
        end
    end


    # Re-scale d²N_dpdcos_sf according to the normalization factor. When dealing
    # with fast push, need to count thermal particle population even for zones
    # upstream of fast push location. Fortunately the number of particles in
    # this case is approximately n₀_ion[i_ion]. Note that this assumes minimal
    # shock modification upstream of the fast push zone -- if (e.g. for
    # non-relativistic shocks) the shock is extensively modified this will not
    # be correct.
    for i in 1:n_grid

        # No thermal particles measured, so add them in proper amount
        if num_crossings[i] == 0 && density_tot[i] > 0
            density_tot[i] += n₀_ion[i_ion]
        end

        # Determine rescaling factor
        if density_tot[i] > 0
            norm_factor = zone_pop[i] / density_tot[i]
        else
            norm_factor = 0.0
        end

        # Re-normalize all nonempty cells of d²N_dpdcos_sf
        for k in 0:psd_max, jθ in 0:psd_max
            if d²N_dpdcos_sf[jθ,k,i] > 1e-99 && norm_factor > 0
                d²N_dpdcos_sf[jθ,k,i] *= norm_factor
            else
                d²N_dpdcos_sf[jθ,k,i] = 1e-99
            end
        end
    end


    # Now, generate d²N_dpdcos_ef and d²n_dpdcos_pf. Do this by transforming the bin
    # *centers* of d²N_dpdcos_sf into each frame, using just the center location to re-bin.
    # TODO: consider doing an exact calculation of bin overlaps rather than just the
    # center-point re-binning that currently happens
    #----------------------------------------------------------------------------
    # Calculate the center points of all the bins to save time later
    cos_center = zeros(0:psd_max)
    for jθ in 0:num_psd_θ_bins
        # Determine current cosines, remembering that psd_θ_bounds has both a
        # linearly-spaced region in cosine and logarithmically-spaced region in θ.
        if jθ > (num_psd_θ_bins - psd_lin_cos_bins)
            cos_hi = psd_θ_bounds[jθ]
            cos_lo = psd_θ_bounds[jθ+1]
        elseif jθ == (num_psd_θ_bins - psd_lin_cos_bins)
            cos_hi = cos(psd_θ_bounds[jθ])
            cos_lo = psd_θ_bounds[jθ+1]
        else
            cos_hi = cos(psd_θ_bounds[jθ])
            cos_lo = cos(psd_θ_bounds[jθ+1])
        end

        # Minus sign needed because finest gradations actually point upstream
        cos_center[jθ] = -(cos_lo + cos_hi)/2
    end
    ptot_center = zeros(0:num_psd_mom_bins+1)
    for k in 0:num_psd_mom_bins
        # Convert from log to linear space
        ptot_center[k] = exp10((psd_mom_bounds[k] + psd_mom_bounds[k+1])/2)
    end

    # Loop extents reflect desired frames:
    #  (1,1)    only plasma frame
    #  (1,2)    both plasma and explosion frame
    #  (2,2)    only explosion frame
    m_max = (i_ion < n_ions) ? 1 : 2
    for m in 1:m_max,  # Loop over plasma and/or ISM frames
        i in 1:n_grid

        # Get Lorentz factor and speed relating the shock and current frames
        if m == 1       # Plasma frame
            γᵤ = γ_sf_grid[i]
            βᵤ = √(1 - 1/γᵤ^2)
        elseif m == 2   # ISM frame
            γᵤ = γ₀
            βᵤ = β₀
        end

        # Loop over cells in d²N_dpdcos_sf
        for k in 0:num_psd_mom_bins
            for jθ in 0:num_psd_θ_bins

                # Skip empty zones
                d²N_dpdcos_sf[jθ,k,i] ≤ 1e-66 && continue

                # Return dN/dp to dN by multiplying by dp
                cell_weight = d²N_dpdcos_sf[jθ,k,i] * Δp[k]

                # Transform the center of the zone into the new frame
                cos_θ_sf    = cos_center[jθ]
                ptot_sf   = ptot_center[k]
                pₓ_sf   = ptot_sf * cos_θ_sf
                etot_sf = hypot(ptot_sf*c, E₀)
                # Get location of center in transformed d²N_dpdcos
                pₓ_Xf = γᵤ * (pₓ_sf - βᵤ*etot_sf/c)
                ptot_Xf = √(ptot_sf^2 - pₓ_sf^2 + pₓ_Xf^2)
                k_Xf = get_psd_bin_momentum(ptot_Xf, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
                jθ_Xf = get_psd_bin_angle(pₓ_Xf, ptot_Xf, psd_bins_per_dec_θ,
                                          num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)

                # And add shock frame number to the appropriate bin in target frame;
                # note that d²N_pf is NOT being converted back to d²N/dp
                if m == 1       # Plasma frame
                    d²N_pf[jθ_Xf, k_Xf, i] += cell_weight
                elseif m == 2   # ISM frame
                    d²N_dpdcos_ef[jθ_Xf, k_Xf, i] += cell_weight / Δp[k_Xf]
                end
            end # loop over angles
        end # loop over ptot
    end # loop over grid locations and reference frames
    #----------------------------------------------------------------------------
    # Transformed d²N_dpdcos calculated


    # Create NetCDF copies of all three arrays, noting that the file name
    # must include the particle species to avoid overwriting.
    #write(tmp1, i_ion)
    ## Shock frame
    #ncfilename = "dNdp_i" // tmp1 // "_sf"
    #output_netcdf(d²N_dpdcos_sf[0:num_psd_θ_bins, 0:num_psd_mom_bins, 1:n_grid],
    #              num_psd_mom_bins, num_psd_θ_bins, n_grid,
    #              ncfilename)
    ## Plasma frame
    #ncfilename = "dNdp_i" // tmp1 // "_pf"
    #output_netcdf(d²N_dpdcos_pf[0:num_psd_θ_bins, 0:num_psd_mom_bins, 1:n_grid],
    #              num_psd_mom_bins, num_psd_θ_bins, n_grid,
    #              ncfilename)
    ## ISM frame
    #ncfilename = "dNdp_i" // tmp1 // "_ef"
    #output_netcdf(d²N_dpdcos_ef[0:num_psd_θ_bins, 0:num_psd_mom_bins, 1:n_grid],
    #              num_psd_mom_bins, num_psd_θ_bins, n_grid,
    #              ncfilename)

    #deallocate(d²N_dpdcos_sf)
    #deallocate(d²N_pf)
    return d²N_dpdcos_ef
end # get_dNdp_2D

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

"""
    get_normalized_dNdp(...)

Computes the actual number of particles in each bin of dN/dp (which is divided by dp,
remember). Computes total area under curve, normalizes it against number of particles
upstream using plasma-frame density and volume, and then uses fractional area of each bin
to determine number of particles in it. Handles non-injected (i.e. thermal) particles
differently than injected ones, due to very small p.f. spread in momenta.

The array that would be `dNdp_cr_pvals` is already set, as `psd_mom_bound`s

### Arguments

- `nc_unit`: unit number for the scratch file holding crossing data

### Returns

- `dNdp_therm`: 3-D array, containing 1-D array for each grid zone of dN/dp for the thermal particles
- `dNdp_therm_pvals`: array of momentum bin boundaries for each row of `dNdp_therm`;
  each row handled separately to maximize resolution of what may be an extremely
  narrow peak at radically different energy from the upstream population
- `dNdp_cr`: 3-D array, containing 1-D array for each grid zone of dN/dp for the population
  of accelerated particles
- `zone_pop`: (Lorentz-invariant) number of particles in each grid zone
"""
function get_normalized_dNdp(
        nc_unit,
        jet_rad_pc, jet_sph_frac, m_ion, n₀_ion, β₀, γ₀, n_ions, do_multi_dNdps,
        num_psd_mom_bins, psd_mom_bounds,
        n_grid, x_grid_cm, uₓ_sk_grid,
        i_iter,
        i_ion,
        γ_sf_grid,
        therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight, num_crossings, n_cr_count,
        psd, psd_lin_cos_bins, num_psd_θ_bins, psd_θ_bounds,
        zone_vol, therm_energy_density, energy_density,
    )

    therm_temp = zeros(0:psd_max, 3)
    therm_pvals_temp = zeros(0:psd_max, 3)
    dNdp_therm_rebin = zeros(0:psd_max, 3)

    # Administrative constants to be used during main computation loops
    #---------------------------------------------------------------------------

    # Number of bins to use in histogram for ptot_pf & cθ_pf; needs to be divisor of
    # num_therm_bins to minimize binning artifacts in this subroutine
    num_hist_bins = num_therm_bins ÷ 2
    jet_rad_cm = ustrip(cm, jet_rad_pc*pc)

    # Get the non-normalized dN/dp's from the crossing arrays, scratch file (if necessary),
    # and the phase space distribution,
    dNdp_therm, dNdp_therm_pvals = get_dNdp_therm(
        num_hist_bins, nc_unit, m_ion, γ₀, β₀, n_grid, γ_sf_grid, i_ion,
        therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight, num_crossings, n_cr_count)

    dNdp_cr = get_dNdp_cr(
        m_ion, psd_lin_cos_bins, γ₀, num_psd_θ_bins, psd_θ_bounds,
        psd, psd_mom_bounds, num_psd_mom_bins, n_grid, γ_sf_grid, i_ion)

    # Now have non-normalized dN/dp for both thermal and CR populations.
    # Determine the total number of particles in each grid zone by using
    # shock frame flux, area, crossing time:
    #     #  =  flux * area * (distance/speed)
    local i_shock
    for i in 1:n_grid
        if iszero(x_grid_cm[i]) || (x_grid_cm[i]*x_grid_cm[i+1] < 0)
            # Either current grid zone is exactly at shock, or current and next grid zones
            # straddle shock
            i_shock = i
            break
        end
    end

    x_grid_cm_diff = diff(x_grid_cm)

    # Work upstream from shock and find volume of each grid zone
    rad_min = jet_rad_cm - x_grid_cm[i_shock] # Inner radius, including fact that upstream has x < 0
    surf_area = Vector{Float64}(undef, n_grid)
    for i in i_shock-1:-1:1
        # outer radius is inner radius plus zone width *in ISM frame*
        rad_max = rad_min + x_grid_cm_diff[i]/γ₀

        # Use rad_max and rad_min to get area of jet surface; will be same in
        # shock frame if shock motion is purely parallel to shock normal
        #surf_area[i] = 4π * ((rad_max+rad_min)/2)^2 * jet_sph_frac
        surf_area[i] = π * (rad_max+rad_min)^2 * jet_sph_frac

        # finally, set rad_min for next cycle through loop
        rad_min = rad_max
    end

    # Work downstream from shock and find volume of each grid zone
    rad_max = jet_rad_cm - x_grid_cm[i_shock]
    for i in i_shock:n_grid
        # inner radius is outer radius minus zone width *in ISM frame*
        rad_min = rad_max - x_grid_cm_diff[i]/γ₀

        # Use rad_max and rad_min to get area of jet surface
        #surf_area[i] = 4π * ((rad_max+rad_min)/2)^2 * jet_sph_frac
        surf_area[i] = π * (rad_max+rad_min)^2 * jet_sph_frac

        # finally, set rad_max for next cycle through loop
        rad_max = rad_min
    end


    # Now calculate the number of particles in each zone, i.e. the total area under dN/dp
    zone_pop = zeros(n_grid)
    for i in 1:n_grid

        # Crossing time in the shock frame, upstream particle flux in the shock frame
        # (conserved quantity throughout the shock structure), and finally number of
        # particles in this grid zone
        dwell_time = x_grid_cm_diff[i] / uₓ_sk_grid[i]

        flux_upstream = γ₀ * n₀_ion[i_ion] * β₀*c

        zone_pop[i] = flux_upstream * surf_area[i] * dwell_time

        #DEBUGLINE
        density_pf = γ₀ * uₓ_sk_grid[1] / (γ_sf_grid[i] * uₓ_sk_grid[i])
        zone_vol[i] = zone_pop[i] / density_pf
    end


    # For each grid zone, integrate the area under dNdp_therm and dNdp_cr to
    # find the non-normalized total area under the two curves. Then
    # normalize each bin of the dNdp curves to get the correct number of
    # particles in each grid zone.
    # IMPORTANT: the differential term is dp, not p^2 dp# Angular components
    # of momentum (the p^2 sin(θ) dφ dθ) have already been handled.
    # Outer loop manages each of the three possible frames of interest:
    #     1  -  Shock frame
    #     2  -  Plasma frame
    #     3  -  ISM frame
    #-------------------------------------------------------------------------
    for m in 1:3, i in 1:n_grid

        # Calculate the total area under the two curves
        area_tot_therm = 0.0
        for j in 0:num_hist_bins-1
            # Increment by dN/dp * dp
            if dNdp_therm[j,i,m] > 1e-99
                area_tot_therm += dNdp_therm[j,i,m] * (dNdp_therm_pvals[j+1,i,m] - dNdp_therm_pvals[j,i,m])
            end
        end
        area_tot_cr    = 0.0
        for j in 0:num_psd_mom_bins
            # Increment by dN/dp * dp
            if dNdp_cr[j,i,m] > 1e-99
                area_tot_cr += dNdp_cr[j,i,m] * (exp10(psd_mom_bounds[j+1]) - exp10(psd_mom_bounds[j]))
            end
        end

        # Re-scale each curve according to the normalization factor. DO NOT include dp, for
        # easier comparison against previous subroutines that performed similar tasks; the
        # dp can be added back in as necessary during integration of dN/dp. When dealing
        # with fast push, need to include thermal particle population even for zones upstream
        # of fast push location. Fortunately, area_tot_therm in this case is approximately
        # density_pf/uₓ, the compressed plasma-frame density divided by the local shock speed.
        #TODO: plot area_tot_therm against local density to see if this holds even once
        # particles start being heated
        if area_tot_therm == 0 && area_tot_cr > 0
            density_pf = n₀_ion[i_ion] * γ₀ * uₓ_sk_grid[1] / (γ_sf_grid[i] * uₓ_sk_grid[i])
            area_tot = density_pf/uₓ_sk_grid[i] + area_tot_cr
        else
            area_tot = area_tot_therm + area_tot_cr
        end

        norm_factor = area_tot > 0 ? zone_pop[i] / area_tot : 0.0

        for j in 0:num_hist_bins-1
            # dN/dp renormalized
            if dNdp_therm[j,i,m] > 1e-99
                dNdp_therm[j,i,m] *= norm_factor
            end
        end
        for j in 0:num_psd_mom_bins
            # dN/dp renormalized
            if dNdp_cr[j,i,m] > 1e-99
                dNdp_cr[j,i,m] *= norm_factor
            end
        end

    end # loop over grid zones and frames
    #-------------------------------------------------------------------------
    # Normalized dN/dp's found


    # Plot the dN/dps as a check-by-eye on their reasonability. Include every grid
    # zone and frame, which is probably more data than necessary in most cases
    #-------------------------------------------------------------------------
    if i_ion == 1
        if do_multi_dNdps
            therm_fileunit = open("mc_dNdp_grid_therm_$i_iter.dat")
            CR_fileunit = open("mc_dNdp_grid_CR_$i_iter.dat")
        else
            therm_fileunit = open("mc_dNdp_grid_therm.dat")
            CR_fileunit = open("mc_dNdp_grid_CR.dat")
        end
    end

    for i in 66:66 #1, n_grid #DEBUGLINE: we only care about this grid zone

        # Thermal particles, all three frames
        if maximum(dNdp_therm[:,i,:]) > 1e-66
            #DEBUGLINE
            for j in 0:num_hist_bins-1
                if dNdp_therm[j,i,2] > 1e-66
                    p_avg = (dNdp_therm_pvals[j,i,2] + dNdp_therm_pvals[j+1,i,2]) /2
                    energy_avg = hypot(m_ion[i_ion]*c^2, p_avg*c)

                    therm_energy_density[i,i_ion] += (dNdp_therm[j,i,2] *
                                                      (dNdp_therm_pvals[j+1,i,2] - dNdp_therm_pvals[j,i,2]) *
                                                      (energy_avg - m_ion[i_ion]*c^2))
                    energy_density[i,i_ion] += (dNdp_therm[j,i,2] *
                                                (dNdp_therm_pvals[j+1,i,2] - dNdp_therm_pvals[j,i,2]) *
                                                (energy_avg - m_ion[i_ion]*c^2))
                end
            end
          j_plot  = 0
          #DEBUGLINE: commented out to save disk space
          #for j in 0:num_hist_bins-1
          #
          #    j_plot += 1
          #    write(therm_fileunit,
          #      i, j_plot,
          #      i_ion,                                         # 1
          #      # Shock frame
          #      log10(dNdp_therm_pvals[j,i,1]),                # 2 (cgs units)
          #      log10(dNdp_therm_pvals[j,i,1] / (mp*c)),       # 3 (natural units)
          #      log10(dNdp_therm[j, i, 1]),                    # 4
          #      # Plasma frame
          #      log10(dNdp_therm_pvals[j,i,2]),                # 5 (cgs units)
          #      log10(dNdp_therm_pvals[j,i,2] / (mp*c)),       # 6 (natural units)
          #      log10(dNdp_therm[j, i, 2]),                    # 7
          #      # ISM frame
          #      log10(dNdp_therm_pvals[j,i,3]),                # 8 (cgs units)
          #      log10(dNdp_therm_pvals[j,i,3] / (mp*c)),       # 9 (natural units)
          #      log10(dNdp_therm[j,i,3]))                      # 10
          #
          #    j_plot += 1
          #    if j < num_hist_bins-1
          #        write(therm_fileunit,
          #          i, j_plot,
          #          i_ion,                                     # 1
          #          # Shock frame
          #          log10(dNdp_therm_pvals[j+1,i,1]),          # 2 (cgs)
          #          log10(dNdp_therm_pvals[j+1,i,1] / (mp*c)), # 3 (natural)
          #          log10(dNdp_therm[j, i, 1]),                # 4
          #          # Plasma frame
          #          log10(dNdp_therm_pvals[j+1,i,2]),          # 5 (cgs)
          #          log10(dNdp_therm_pvals[j+1,i,2] / (mp*c)), # 6 (natural)
          #          log10(dNdp_therm[j, i, 2]),                # 7
          #          # ISM frame
          #          log10(dNdp_therm_pvals[j+1,i,3]),          # 8 (cgs)
          #          log10(dNdp_therm_pvals[j+1,i,3] / (mp*c)), # 9 (natural)
          #          log10(dNdp_therm[j,i,3]))                  # 10
          #    end
          #end # loop over num_hist_bins
          #
          #print_plot_vals(therm_fileunit)

        end  # check on existence of thermal particles


        # Cosmic rays
        if maximum(dNdp_cr[:,i,:]) > 1e-66

            # Convert thermal particles, if any,  to bins used for cosmic rays
            if maximum(dNdp_therm[:,i,:]) > 1e-66
                #DEBUGLINE
                for j in 0:num_psd_mom_bins
                    if dNdp_cr[j,i,2] > 1e-66
                        p_avg = √(exp10(psd_mom_bounds[j+1]) - exp10(psd_mom_bounds[j]))
                        energy_avg = hypot(m_ion[i_ion]*c^2, p_avg*c)

                        energy_density[i,i_ion] += (dNdp_cr[j,i,2] *
                                                    (exp10(psd_mom_bounds[j+1]) - exp10(psd_mom_bounds[j])) *
                                                    (energy_avg - m_ion[i_ion]*c^2))
                    end
                end
                therm_temp       .= dNdp_therm[:,i,:]
                therm_pvals_temp .= dNdp_therm_pvals[:,i,:]

                dNdp_therm_rebin = rebin_dNdp_therm(
                    num_hist_bins, therm_temp, therm_pvals_temp,
                    num_psd_mom_bins, psd_mom_bounds)
            else
                fill!(dNdp_therm_rebin, 1e-99)
            end

            j_plot  = 0
            for j in 0:num_psd_mom_bins
                j_plot += 1
                write(CR_fileunit,
                      i, j_plot,
                      i_ion,                                            # 1
                      psd_mom_bounds[j],                                # 2 (cgs units)
                      psd_mom_bounds[j] - log10(mp*c),                  # 3 (natural units)
                      # Shock frame, Plasma frame, ISM frame
                      log10.(dNdp_cr[j, i, 1:3]),                       # 4-6
                      # Summed therm+CR dN/dps in all three frames
                      log10.(dNdp_cr[j,i,1:3]+dNdp_therm_rebin[j,1:3])) # 7-9

                j_plot += 1
                if j < num_psd_mom_bins
                    write(CR_fileunit,
                      i, j_plot,
                      i_ion,                                            # 1
                      psd_mom_bounds[j+1],                              # 2 (cgs)
                      psd_mom_bounds[j+1] - log10(mp*c),                # 3 (nat.)
                      # Shock frame, Plasma frame, ISM frame
                      log10.(dNdp_cr[j, i, 1:3]),                       # 4-6
                      # Summed therm+CR dN/dps in all three frames
                      log10.(dNdp_cr[j,i,1:3]+dNdp_therm_rebin[j,1:3])) # 7-9
                end
            end # loop over num_psd_mom_bins

            print_plot_vals(CR_fileunit)

        end # check on existence of CRs

    end # loop over grid zones

    i_ion == n_ions && close(therm_fileunit)
    i_ion == n_ions && close(CR_fileunit)
    #-------------------------------------------------------------------------
    # end of plotting section

    return dNdp_therm, dNdp_therm_pvals, dNdp_cr, zone_pop
end # get_normalized_dNdp


"""
    get_dNdp_therm(...)

Calculate dN/dp (NOT normalized) for particles that haven't been injected
into acceleration process.

First create arrays to hold crossing information from this ion species, and fill
them with data from the scratch file.
Next, loop over grid locations, performing three main tasks:
1. Transform shock frame values into plasma frame
2. Seek out maximum and minimum momenta in plasma frame
3. Bin crossings according to plasma frame total momentum
Also perform the above process for the ISM frame, and store shock frame
values for comparison against earlier results.

### Arguments

- `num_hist_bins`: number of bins to use in histograms
- `nc_unit`: unit number of scratch file for crossing data

### Returns

- `dNdp_therm`: 3-D array, containing 1-D array for each grid zone of dN/dp for the thermal particles
- `dNdp_therm_pvals`: array of momentum bin boundaries for each row of `dNdp_therm`;
  each row handled separately to maximize resolution of what may be an extremely
  narrow peak at radically different energy from the upstream population
"""
function get_dNdp_therm(
        num_hist_bins, nc_unit,
        m_ion, γ₀, β₀,
        n_grid, γ_sf_grid,
        i_ion, therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight, num_crossings, n_cr_count,
    )

    # Set a constant to be used repeatedly
    E₀ = m_ion[i_ion] * c^2  # Rest mass-energy of the current particle species

    # Also "zero" out the two output arrays to prevent issues later
    dNdp_therm       = fill(1e-99, (0:psd_max,n_grid,3))
    dNdp_therm_pvals = fill(1e-99, (0:psd_max,n_grid,3))

    # Allocate histogram arrays to be used during subroutine
    ptot_sk_bins = zeros(num_hist_bins)
    ptot_sk_vals = zeros(num_hist_bins)
    ptot_pf_bins = zeros(num_hist_bins)
    ptot_pf_vals = zeros(num_hist_bins)
    ptot_ef_bins = zeros(num_hist_bins)
    ptot_ef_vals = zeros(num_hist_bins)

    cθ_sk_bins = zeros(num_hist_bins+1)
    cθ_pf_bins = zeros(num_hist_bins+1)
    cθ_ef_bins = zeros(num_hist_bins+1)

    cθ_sk_weight = zeros(num_hist_bins)
    cθ_pf_weight = zeros(num_hist_bins)
    cθ_ef_weight = zeros(num_hist_bins)


    # Create the arrays to hold the crossing information, and fill them with
    # data from the scratch file
    #-----------------------------------------------------------------------
    max_cross = maximum(num_crossings)
    therm_pₓ = zeros(max_cross,n_grid)
    therm_pt = zeros(max_cross,n_grid)
    therm_weight = zeros(max_cross,n_grid)


    # Fill the arrays with data from the crossing arrays and (if needed) from the scratch file
    ntot_crossings = sum(num_crossings)

    ##rewind(nc_unit) # XXX Fortran holdover

    # n_cross_fill needs to be an array because we will be skipping around the grid
    # as we move through the data rather than filling a single zone at a time
    n_cross_fill = zeros(Int, n_grid)

    # Handle crossings stored within the crossing arrays
    for i in 1:n_cr_count
        i_grid = therm_grid[i]

        n_cross_fill[i_grid] += 1

        # Note: ordering of coordinates chosen to make memory accesses in next loop faster
        therm_pₓ[n_cross_fill[i_grid], i_grid] = therm_pₓ_sk[i]
        therm_pt[n_cross_fill[i_grid], i_grid] = therm_ptot_sk[i]
        therm_weight[n_cross_fill[i_grid], i_grid] = therm_weight[i]
    end

    if ntot_crossings > na_cr
        # Need to go into the scratch file for the remainder of the crossings
        for i in 1:ntot_crossings-na_cr
            i_grid, _, pₓ_sk, ptot_sk, cell_weight = read(nc_unit)

            n_cross_fill[i_grid] += 1

            # Note: ordering of coordinates chosen to make memory accesses in next loop faster
            therm_pₓ[n_cross_fill[i_grid], i_grid] = pₓ_sk
            therm_pt[n_cross_fill[i_grid], i_grid] = ptot_sk
            therm_weight[n_cross_fill[i_grid], i_grid] = cell_weight
        end
    end

    #-------------------------------------------------------------------------
    # Arrays created and read in


    # Main loop over grid locations
    #-----------------------------------------------------------------------
    for i in 1:n_grid
        if iszero(num_crossings[i])
            dNdp_therm[:,i,:] .= 1e-99
            continue      # Ignore zones that no thermal particles crossed
        end

        # 1. Transform crossing data from shock frame into plasma & ISM frames
        #------------------------------------------------------------------------
        # Get current flow speed
        γᵤ = γ_sf_grid[i]
        βᵤ = √(1 - 1/γᵤ^2)

        # Create arrays to hold plasma frame and ISM frame values, then
        # initialize them. Since the arrays are handled independently for each
        # grid zone, all positions in the arrays should be filled with correct
        # data. However, initializing them to non-physical values serves as an
        # additional check against error.
        ptot_pf = fill(-1.0, num_crossings[i])
        cθ_pf = fill(-2.0, num_crossings[i])
        ptot_ef = fill(-1.0, num_crossings[i])
        cθ_ef = fill(-2.0, num_crossings[i])

        # Set up min and max values for total momentum and cos(θ) in the shock frame
        cθ_sk_min =  2.0
        cθ_sk_max = -2.0
        ptot_sk_min =  1e99
        ptot_sk_max = -1e99

        # Convert information about shock frame values into plasma and ISM frames. Even
        # though only total momentum is needed for calculating dN/dp, histogram of pitch
        # angle values can provide an additional check on whether the distribution is
        # isotropic in the other frames. So convert ptot_sk & pₓ_sk into ptot_pf/cθ_pf and
        # ptot_ef/cθ_ef. Also, find minimum and maximum values of ptot_sk and cθ_sk.
        for j in 1:num_crossings[i]
            pₓ_sk = therm_pₓ[j, i]
            ptot_sk = therm_pt[j, i]
            ptot_sk_max = max(ptot_sk, ptot_sk_max)
            ptot_sk_min = min(ptot_sk, ptot_sk_min)

            cθ_sk = pₓ_sk / ptot_sk
            cθ_sk_max = max(cθ_sk, cθ_sk_max)
            cθ_sk_min = min(cθ_sk, cθ_sk_min)

            etot_sk = hypot(ptot_sk*c, E₀)

            pₓ_pf = γᵤ * (pₓ_sk - βᵤ*etot_sk/c)
            ptot_pf = √(ptot_sk^2 - pₓ_sk^2 + pₓ_pf^2)

            pₓ_ef = γ₀ * (pₓ_sk - β₀*etot_sk/c)
            ptot_ef = √(ptot_sk^2 - pₓ_sk^2 + pₓ_ef^2)

            ptot_pf[j] = ptot_pf
            cθ_pf[j] = pₓ_pf / ptot_pf

            ptot_ef[j] = ptot_ef
            cθ_ef[j] = pₓ_ef / ptot_ef
        end

        # Error checks
        count(ptot_pf .== -1) > 0 && error("Unfilled zones in ptot_pf")
        count(cθ_pf .== -2) > 0 && error("Unfilled zones in cθ_pf")
        count(ptot_ef .== -1) > 0 && error("Unfilled zones in ptot_ef")
        count(cθ_ef .== -2) > 0 && error("Unfilled zones in cθ_ef")
        #------------------------------------------------------------------------
        # Plasma frame, ISM frame values determined; Shock frame extrema found.


        # 2. Seek out maximum and minimum total momentum in plasma frame and ISM
        # frame. Create histogram bins in both momentum and angle as desired.
        #------------------------------------------------------------------------
        # Already have extrema for shock frame
        Δptot_sk = (ptot_sk_max - ptot_sk_min) / num_hist_bins

        ptot_pf_max = maximum(ptot_pf, dims=1)
        ptot_pf_min = minimum(ptot_pf, dims=1)
        Δptot_pf = (ptot_pf_max - ptot_pf_min) / num_hist_bins

        ptot_ef_max = maximum(ptot_ef, dims=1)
        ptot_ef_min = minimum(ptot_ef, dims=1)
        Δptot_ef = (ptot_ef_max - ptot_ef_min) / num_hist_bins

        Δcθ_sk = (cθ_sk_max - cθ_sk_min) / num_hist_bins

        for k in 1:num_hist_bins
            ptot_sk_bins[k] = ptot_sk_min + (k-1)*Δptot_sk
            ptot_pf_bins[k] = ptot_pf_min + (k-1)*Δptot_pf
            ptot_ef_bins[k] = ptot_ef_min + (k-1)*Δptot_ef

            cθ_sk_bins[k] = cθ_sk_min + (k-1) * Δcθ_sk
            cθ_pf_bins[k] = -1 + (k-1)*2/num_hist_bins
            cθ_ef_bins[k] = -1 + (k-1)*2/num_hist_bins
        end
        cθ_sk_bins[num_hist_bins+1] = cθ_sk_max
        cθ_pf_bins[num_hist_bins+1] = 1.0
        cθ_ef_bins[num_hist_bins+1] = 1.0

        ptot_sk_vals[1:num_hist_bins] = 0.0
        ptot_pf_vals[1:num_hist_bins] = 0.0
        ptot_ef_vals[1:num_hist_bins] = 0.0
        cθ_sk_weight[1:num_hist_bins] = 0.0
        cθ_pf_weight[1:num_hist_bins] = 0.0
        cθ_ef_weight[1:num_hist_bins] = 0.0
        #------------------------------------------------------------------------
        # Histogram bins set


        # 3. Bin particles according to total momentum and cos(θ) in all three frames
        #
        # For plasma and ISM frames, divide therm_cwt and by γ of flow speed as
        # part of Lorentz transformation of phase-space density (Rybicki & Lightman, p.146)
        #------------------------------------------------------------------------
        for j in 1:num_crossings[i]
            # Determine bin in shock frame momentum; add value to correct bin
            k = min(num_hist_bins, trunc(Int, (therm_pt[j,i]-ptot_sk_min)/Δptot_sk) + 1)
            ptot_sk_vals[k] += therm_weight[j,i]

            # Determine bin in plasma frame momentum; add value to correct bin
            k = min(num_hist_bins, trunc(Int, (ptot_pf[j]-ptot_pf_min)/Δptot_pf) + 1)
            ptot_pf_vals[k] += therm_weight[j,i]/γᵤ

            # Determine bin in ISM frame momentum; add value to correct bin
            k = min(num_hist_bins, trunc(Int, (ptot_ef[j]-ptot_ef_min)/Δptot_ef) + 1)
            ptot_ef_vals[k] += therm_weight[j,i]/γ₀

            # Determine shock frame bin of cos(θ); add value to correct bin
            cθ_sk = therm_pₓ[j,i] / therm_pt[j,i]
            k = min(num_hist_bins, trunc(Int, (cθ_sk - cθ_sk_min)/Δcθ_sk) + 1)
            k = max(k, 1) # odd floating point error makes some k < 1
            cθ_sk_weight[k] += therm_weight[j,i]

            # Determine plasma frame bin of cos(θ); add value to correct bin
            k = min(num_hist_bins, trunc(Int, (cθ_pf[j] + 1) * num_hist_bins/2) + 1)
            cθ_pf_weight[k] += therm_weight[j,i]/γᵤ

            # Determine ISM frame bin of cos(θ); add value to correct bin
            k = min(num_hist_bins, trunc(Int, (cθ_ef[j] + 1) * num_hist_bins/2) + 1)
            cθ_ef_weight[k] += therm_weight[j,i]/γ₀
        end

        # Finish off by converting pt**_vals into dN/dp, which requires dividing
        # by the momentum spread of the bin
        for k in 1:num_hist_bins
            ptot_sk_lo = ptot_sk_bins[k]
            ptot_pf_lo = ptot_pf_bins[k]
            ptot_ef_lo = ptot_ef_bins[k]
            if k == num_hist_bins
                ptot_sk_hi = ptot_sk_max
                ptot_pf_hi = ptot_pf_max
                ptot_ef_hi = ptot_ef_max
            else
                ptot_sk_hi = ptot_sk_bins[k+1]
                ptot_pf_hi = ptot_pf_bins[k+1]
                ptot_ef_hi = ptot_ef_bins[k+1]
            end

            ptot_sk_vals[k] /= (ptot_sk_hi - ptot_sk_lo)
            ptot_pf_vals[k] /= (ptot_pf_hi - ptot_pf_lo)
            ptot_ef_vals[k] /= (ptot_ef_hi - ptot_ef_lo)
        end

        ptot_sk_tot = sum(ptot_sk_vals, 1)
        ptot_pf_tot = sum(ptot_pf_vals, 1)
        ptot_ef_tot = sum(ptot_ef_vals, 1)
        ∑cθ_sk_weight = sum(cθ_sk_weight)
        ∑cθ_pf_weight = sum(cθ_pf_weight)
        ∑cθ_ef_weight = sum(cθ_ef_weight)
        #------------------------------------------------------------------------
        # Particles binned


        # Finally, copy pt**_vals and pt**_bins into the dNdp arrays; shift down by
        # one bin to make both thermal and cosmic ray dN/dp's start at bin 0.
        # Also, set zero values of dNdp_therm to a small but nonzero number.
        # Note on third dimension of dNdp_therm:
        #    1  -  Shock frame
        #    2  -  Plasma frame
        #    3  -  ISM frame
        for k in 1:num_hist_bins
            # Shock frame
            dNdp_therm[k-1, i, 1] = iszero(ptot_sk_vals[k]) ? 1e-99 : ptot_sk_vals[k]

            # Plasma frame
            dNdp_therm[k-1, i, 2] = iszero(ptot_pf_vals[k]) ? 1e-99 : ptot_pf_vals[k]

            # ISM frame
            dNdp_therm[k-1, i, 3] = iszero(ptot_ef_vals[k]) ? 1e-99 : ptot_ef_vals[k]

            # Set bin boundaries
            dNdp_therm_pvals[k-1, i, 1] = ptot_sk_bins[k]
            dNdp_therm_pvals[k-1, i, 2] = ptot_pf_bins[k]
            dNdp_therm_pvals[k-1, i, 3] = ptot_ef_bins[k]
        end
        dNdp_therm_pvals[num_hist_bins, i, 1] = ptot_sk_max
        dNdp_therm_pvals[num_hist_bins, i, 2] = ptot_pf_max
        dNdp_therm_pvals[num_hist_bins, i, 3] = ptot_ef_max

        # The next section is a set of tests and checks on thermal population.
        # Fit Maxwell-Boltzmann distribution to plasma frame dN/dp, and determine temperature
        # of distribution. Also, calculate pressure using integral formulation of
        # equation (17) from Ellison & Reynolds (1991) [1991ApJ...378..214E].
        #------------------------------------------------------------------------
        #num_skipped = 0
        #∑p²      = 0.0
        #∑p²_f    = 0.0
        #∑p²_lnp² = 0.0
        #∑pfth    = 0.0
        #∑f       = 0.0
        #∑lnp²    = 0.0
        #for k in 0:num_hist_bins-1 # now using shifted bins
        #    # Here, f is the natural logarithm of dN/dp because fitting occurs in log-log space
        #    if dNdp_therm(k, i, 2) > 0
        #        f = log(dNdp_therm[k, i, 2])
        #    else
        #        f = 0.0
        #        num_skipped += 1
        #        continue
        #    end
        #    p² = dNdp_therm_pvals[k,i,2]^2
        #    ∑p²      += p²
        #    ∑p²_f    += p² * f
        #    ∑p²_lnp² += p² * log(p²)
        #    ∑pfth    += p²^2
        #    ∑f       += f
        #    ∑lnp²    += log(p²)
        #end
        #num_bins = num_hist_bins - num_skipped
        #lnA    = (∑p²*∑p²_f - ∑p²*∑p²_lnp² - ∑pfth*∑f + ∑pfth*∑lnp²) / (∑p²^2 - num_bins*∑pfth)
        #expfac = (-num_bins*∑p²_f + ∑f*∑p² - ∑p²*∑lnp² + num_bins*∑p²_lnp²) / (∑p²^2 - num_bins*∑pfth)
        #temp   = - 1 / (2 * m_ion[numion] * kB * expfac)
        #n0     = exp(lnA) * (m_ion[numion]*kB*temp)^(3//2) * √(π/2)

        ## Generate values of fitted M-B distribution
        #mb_vals = zeros(num_hist_bins)
        #for k in 1:num_hist_bins
        #    p² = dNdp_therm_pvals[k-1,i,2]^2
        #    mb_vals[k] = exp(lnA + log(p²) + expfac*p²)
        #end
        #mb_vals ./= sum(mb_vals) # normalize

        ## Calculate pressure using integral formula.
        ## Plot dN/dp, fitted M-B distribution, and pressure
        #j_plot   = 0
        #pressure = 0.0
        #i_unit  = 700 + i
        #for k in 1:num_hist_bins
        #    ptot_pf = √(dNdp_therm_pvals[k-1,i,2] * dNdp_therm_pvals[k,i,2])
        #    Δptot_pf = dNdp_therm_pvals[k,i,2] - dNdp_therm_pvals[k-1,i,2]
        #
        #    γₚ_pf = √(1 + (ptot_pf*c/E₀)^2)
        #    pressure += (ptot_pf * ptot_pf/(m_ion[i_ion]*γₚ_pf) *
        #                 dNdp_therm[k-1,i,2] * Δptot_pf/3)
        #
        #    j_plot += 1
        #    write(i_unit, i, j_plot,
        #          log10(dNdp_therm_pvals[k-1,i,2]),    # 1
        #          dNdp_therm[k-1,i,2]/ptot_pf_tot,     # 2
        #          cθ_pf_bins[k],                       # 3
        #          cθ_pf_weight[k]/∑cθ_pf_weight,       # 4
        #          cθ_sk_bins[k],                       # 5
        #          cθ_sk_weight[k]/∑cθ_sk_weight,       # 6
        #          mb_vals[k],                          # 7
        #          pressure)                            # 8
        #    j_plot += 1
        #    if k < num_hist_bins
        #        write(i_unit, i, j_plot,
        #              log10(dNdp_therm_pvals[k,i,2]),  # 1
        #              dNdp_therm[k-1,i,2]/ptot_pf_tot,   # 2
        #              cθ_pf_bins[k+1],                 # 3
        #              cθ_pf_weight[k]/∑cθ_pf_weight,   # 4
        #              cθ_sk_bins[k+1],                 # 5
        #              cθ_sk_weight[k]/∑cθ_sk_weight,   # 6
        #              mb_vals[k],                      # 7
        #              pressure)                        # 8
        #    end
        #end # loop over num_hist_bins
        #print_plot_vals(i_unit)
        #------------------------------------------------------------------------
        # End of tests/checks section

    end # loop over grid zones
end # get_dNdp_therm

"""
    rebin_dNdp_therm(num_hist_bins, dNdp_therm, dNdp_therm_pvals, num_psd_mom_bins, psd_mom_bounds)

This subroutine takes the distribution of thermal particles and re-bins it in the bins used
for the cosmic rays.
Note that the thermal distribution is provided (and re-binned separately) in
all three frames: shock, plasma, and ISM.

### Arguments
- `num_hist_bins`: number of bins used for the thermal distribution
- `dNdp_therm`: distribution of thermal particles in all three frames
- `dNdp_therm_pvals`: arrays holding bin boundary values for the thermal distributions;
  as with the distributions themselves, they change based on grid location
- `num_psd_mom_bins`: number of bins used for the cosmic ray distribution
- `psd_mom_bounds`: boundaries of the bins for the cosmic ray distribution
  Unlike the thermal boundaries, these are constant everywhere on the grid.

### Returns
Re-binned histogram of the thermal distribution
"""
function rebin_dNdp_therm(num_hist_bins, dNdp_therm, dNdp_therm_pvals, num_psd_mom_bins, psd_mom_bounds)

    # Zero out the output array
    rebin_dNdp_therm = fill(1e-99, (0:psd_max,3))

    # Convert psd_mom_bounds from log space into linear space
    lin_bounds = exp10.(psd_mom_bounds)

    dN_therm       = OffsetVector{Float64}(undef, 0:psd_max)
    therm_pvals    = OffsetVector{Float64}(undef, 0:psd_max)
    dN_therm_rebin = OffsetVector{Float64}(undef, 0:psd_max)

    # Loop over the various frames.
    #     1  -  Shock
    #     2  -  Plasma
    #     3  -  ISM
    for m in 1:3

        # Convert the arrays from 2-D down to 1-D to make referring to them easier.
        # Also, turn dN/dp into dN, noting that dNdp_therm_pvals is already in cgs units.
        dN_therm    .= dNdp_therm[:,m]
        therm_pvals .= dNdp_therm_pvals[:,m]
        for i in 0:num_hist_bins - 1
            dN_therm[i] *= max(1e-99, (therm_pvals[i+1] - therm_pvals[i]))
        end

        # Find the lower bound on where the thermal distribution might fall, to save cycles later
        i_lo = findfirst(i -> lin_bounds[i] > therm_pvals[0], 1:num_psd_mom_bins) - 1

        fill!(dN_therm_rebin, 1e-99)
        # Now, re-bin the thermal distribution. Loop over the bins of dN_therm,
        # assigning each one to the CR distribution bin(s) where it belongs.
        for j in 0:num_hist_bins-1

            dN_therm[j] ≤ 1e-99 && continue # skip empty bins
            dN_remaining = dN_therm[j]

            p_lo = therm_pvals[j]
            p_hi = therm_pvals[j+1]

            p_bottom = p_lo
            p_length = p_hi - p_lo

            for i in i_lo:num_psd_mom_bins-1

                # Determine the fraction of p_length that the current cell represents.
                frac_p_length = (lin_bounds[i+1] - p_bottom) / p_length

                # If the top of the cell is above the top of the thermal bin, the remainder
                # of the thermal bin fits entirely in this cell of the CR dN/dp.
                # Exit the loop and move to the next thermal bin.
                if lin_bounds[i+1] > p_hi
                    dN_therm_rebin[i] += dN_remaining
                    break
                elseif frac_p_length < 0
                    # The current cell of the CR histogram is *entirely* below the bin
                    # of the thermal distribution. Skip and move on to the next cycle
                    continue
                end

                # If the code reached this point, the thermal bin extends beyond the
                # upper edge of the current cell. Add the fraction *within* the cell to
                # dN_therm_rebinned and reset for the next cycle
                dN_therm_rebin[i] += dN_therm[j] * frac_p_length
                dN_remaining      -= dN_therm[j] * frac_p_length
                p_bottom = lin_bounds[i+1]

            end # loop over cells in CR histogram

        end # loop over thermal bins

        # Convert the re-binned thermal distribution from dN into dN/dp.
        # Then copy it into the output array.
        for i in 0:num_psd_mom_bins-1
            if dN_therm_rebin[i] > 1e-90
                dN_therm_rebin[i] /= (lin_bounds[i+1]-lin_bounds[i])
            end
        end
        rebin_dNdp_therm[:,m] .= dN_therm_rebin

    end # loop over frames

    return rebin_dNdp_therm
end # rebin_dNdp_therm
end # module
