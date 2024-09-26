using .constants: c_cgs
using .parameters: na_photons
using .io: print_plot_vals

"""
This subroutine takes the output from the various emission subroutines, which are called one
time for each grid zone of emission, and combines all emission profiles into a single spectrum
for emission process in each spectral region.

For pion decay and inverse Compton (if used), if more than one species of nucleus or more
than one photon field was present, write the summed spectra (collected in pion_photon_sum
and ic_photon_sum) to the file and use those spectra instead.

It then combines all num_UpS_shells + num_DwS_shells spectra into a single energy spectrum
spanning all processes, summing the various emission spectra into a unified profile.

### Arguments
- n_shells: number of photon shells, i.e. number of spectra in file photon_***.dat
- n_spec_endpoints: grid boundaries between spectral regions
- n_photon_***: one more than number of energy bins used for emission spectrum
- photon_energy_***_MeV: min/max value of final (summed) emission spectrum's
- energy range, in units of MeV
- n_bins_per_decade: number of energy bins per decade of photon spectrum.
- photon_***_min_MeV: minimum value of specified emission spectrum's energy range,
  in units of MeV

### Returns
None

!!! warning
    This is currently set up to handle only single-iteration runs.
    If doing more than one iteration there will be more lines in the
    photon_***_grid.dat files than assumed by the current code.
"""
function get_summed_emission(
        n_shells, n_shell_endpoints,
        n_photon_pion, n_photon_synch, n_photon_IC, photon_energy_min_MeV, photon_energy_max_MeV,
        n_bins_per_decade, photon_pion_min_MeV, photon_synch_min_MeV, photon_ic_min_MeV,
        aa_ion, n_ions,
        γ_sf_grid, utot_grid, β_ef_grid, γ_ef_grid, n_grid,
        n_pion_specs, n_IC_specs, energy_pion_MeV, energy_IC_MeV, pion_photon_sum, ic_photon_sum,
    )

    # Check pion decay and inverse Compton for presence of multiple spectra due to more than
    # one nucleus species or more than one photon field. If needed, add an extra set of
    # spectra to the emission files and use the summed spectra directly for the rest of the
    # processing. This means setting a logical flag to be checked in the next section.
    #-------------------------------------------------------------------------
    # Pion decay first
    skip_pion_read = false
    i_grid_start = zeros(Int, n_ions)
    i_grid_end = zeros(Int, n_ions)
    if n_pion_specs > 1
        skip_pion_read = true
        pion_emission()
    end # check on n_pion_specs


    # Now do inverse Compton
    skip_IC_read = false
    if n_ic_specs > 1
        skip_IC_read = true
        inverse_compton_emission()
    end  # check on n_ic_specs
    #-------------------------------------------------------------------------
    # Summed spectra written if multiple spectra present


    energy_MeV_in     = zeros(0:na_photons)
    energy_MeV_pion   = zeros(0:na_photons)
    energy_MeV_synch  = zeros(0:na_photons)
    energy_MeV_IC     = zeros(0:na_photons)
    energy_mid        = zeros(0:na_photons)
    photon_flux_in    = Matrix{Float64}(undef, (na_photons, n_grid))
    photon_flux_pion  = Matrix{Float64}(undef, (na_photons, n_grid))
    photon_flux_synch = Matrix{Float64}(undef, (na_photons, n_grid))
    photon_flux_IC    = Matrix{Float64}(undef, (na_photons, n_grid))

    # Determine number of emission types based on whether electrons were used during this run.
    kmax = (count(aa_ion .< 1) < 1) ? 1 : 3

    # Read in data from all threee emission types: pion decay, synchrotron, and inverse Compton
    #-------------------------------------------------------------------------
    for k in 1:kmax

        # First, if the output files are still open, close and reopen them to start reading
        # at the beginning of the file. Also, get the zone number of the first block of
        # emission data (not necessarily zone 1 because of the location of the free escape
        # boundary during the MC iterations)
        if k == 1 # Pion decay
            inquire(file="./photon_pion_decay_grid.dat", opened=lopen, number=j_unit)
            lopen && close(unit=j_unit)

            open(newunit=grid_file_unit, status="unknown", file="./photon_pion_decay_grid.dat", position="append")
            read(grid_file_unit, i_grid_start[k])
            rewind(grid_file_unit)

        elseif k == 2 # Synchrotron
            inquire(file="./photon_synch_grid.dat",opened=lopen,number=j_unit)
            lopen && close(unit=j_unit)

            open(newunit=grid_file_unit, status="unknown", file="./photon_synch_grid.dat", position="append")
            read(grid_file_unit, i_grid_start[k])
            rewind(grid_file_unit)

        elseif k == 3 # Inverse Compton, using CMB (or other) photon field
            inquire(file="./photon_IC_grid.dat",opened=lopen,number=j_unit)
            lopen && close(unit=j_unit)

            open(newunit=grid_file_unit, status="unknown", file="./photon_IC_grid.dat", position="append")
            read(grid_file_unit, i_grid_start[k])
            rewind(grid_file_unit)
        end

        # Next, set the number of bins to be used in the summation. The "-1"s are
        # because the subroutines don't write out the final value of each spectrum
        if k == 1
            n_photon = n_photon_pion - 1
        elseif k == 2
            n_photon = n_photon_synch - 1
        elseif k == 3
            n_photon = n_photon_IC - 1
        end

        # Initialize array that will hold fluxes, knowing that it is in log space
        fill!(photon_flux_in, -99.0)

        # Loop over grid zones, building array photon_flux_in
        # WARNING: photon_flux_in is actually photon flux per unit ln(energy). To
        # get true photon flux, must multiply by width of energy bin in log space.
        # Additionally, photon_flux_in is logarithm of true photon number.
        #--------------------------------------------------------------------------
        # If more than one pion decay or IC spectrum was generated, handle this part differently.
        if k == 1 && skip_pion_read
            photon_flux_in[1:n_photon, i_grid_start[k]:n_grid] .= pion_photon_sum[1:n_photon, i_grid_start[k]:n_grid]
            energy_MeV_in[1:n_photon] .= log10.(energy_pion_MeV[1:n_photon])

        elseif k == 3 && skip_IC_read

            photon_flux_in[1:n_photon, i_grid_start[k]:n_grid] .= ic_photon_sum[1:n_photon, i_grid_start[k]:n_grid]
            energy_MeV_in[1:n_photon] .= log10.(energy_IC_MeV[1:n_photon])

        elseif k == 2 # Make sure different formats of synch. and other emission processes are managed

            # First spectrum read in gets slightly special treatment, in that its
            # energy bins are assigned directly to energy_γ_MeV_in without error checking.
            for j in 1:n_photon
                read(grid_file_unit,  # photon_***.dat
                     idum1, idum2,
                     xdum1,                             # 1 Log10(keV)
                     photon_flux_in[j,i_grid_start[k]], # 2 Log10(photons/(cm²⋅sec))
                     energy_MeV_in[j])                  # 3 Log10(MeV)
            end
            read(grid_file_unit,"(2I5,8ES12.3E2)") # Advance past output from print_plot_vals


            # Loop over remaining grid zones
            for i in i_grid_start[k]+1:n_grid
                # Each spectra read in after the first is given a simple error check -- do
                # the energy bins match previously read-in energy bins? -- before being
                # included in photon_flux_in
                for j in 1:n_photon
                    read(grid_file_unit, # photon_***.dat
                         i_grid_end[k], idum2,
                         xdum1,                     # 1 Log10(keV)
                         photon_flux_tmp,           # 2 Log10(photons/(cm²⋅sec))
                         energy_MeV_tmp,            # 3 Log10(MeV)
                         iostat=m)

                    # If using a downstream FEB, end-of-file will be reached before i = n_grid.
                    # In that case, exit the loop early.
                    m < 0 && break

                    energy_MeV_tmp != energy_MeV_in[j] && error(
                        "Energy bins for read-in spectra do not match. ",
                        "energy_MeV_in[$j] = ",energy_MeV_in[j],". energy_MeV_tmp[$j] = ",energy_MeV_tmp)

                    photon_flux_in[j,i] = photon_flux_tmp

                end

                # Again, exit the loop if DwS FEB means some grid zones don't have emission recorded
                m < 0 && break

                read(grid_file_unit,"(2I5,8ES12.3E2)") # Advance past output of print_plot_vals

            end # loop over grid zones

        else

            # First spectrum read in gets slightly special treatment, in that its energy
            # bins are assigned directly to energy_MeV_in without error checking.
            for j in 1:n_photon
                read(grid_file_unit,  # photon_***.dat
                     idum1, idum2,
                     xdum1,                                 # 1 nucleus species/photon field
                     xdum2,                                 # 2 Log10(keV)
                     photon_flux_in[j,i_grid_start[k]],     # 3 Log10(photons/(cm²⋅sec))
                     energy_MeV_in[j])                      # 4 Log10(MeV)
            end
            read(grid_file_unit,"(2I4,F5.1,48ES12.3E2)") # Advance past output of print_plot_vals


            # Loop over remaining grid zones
            for i in i_grid_start[k]+1:n_grid
                # Each spectra read in after the first is given a simple error check -- do
                # the energy bins match previously read-in energy bins? -- before being
                # included in photon_flux_in
                for j in 1:n_photon
                    read(grid_file_unit,        # photon_***.dat
                         i_grid_end[k], idum2,
                         xdum1,                 # 1 ion species/photon field
                         xdum2,                 # 2 Log10(keV)
                         photon_flux_tmp,       # 3 Log10(photons/(cm²⋅sec))
                         energy_MeV_tmp,        # 4 Log10(MeV)
                         iostat=m)

                    # If using a downstream FEB, end-of-file will be reached before i = n_grid.
                    # In that case, exit the loop early.
                    m < 0 && break

                    energy_MeV_tmp != energy_MeV_in[j] && error(
                        "Energy bins for read-in spectra do not match. ",
                        "energy_MeV_in[$j] = ",energy_MeV_in[j],". energy_MeV_tmp[$j] = ",energy_MeV_tmp)

                    photon_flux_in[j,i] = photon_flux_tmp
                end

                # Again, exit the loop if DwS FEB means some grid zones don't have emission recorded
                m < 0 && break

                read(grid_file_unit) # Advance past output of print_plot_vals
            end # loop over grid zones

        end # check on skip_***_read


        # Assign the temporary array to the emission-specific array for further handling.
        # Additionally, convert it from log space to linear space.
        if k == 1
            photon_flux_pion  .= exp10.(photon_flux_in)
            energy_MeV_pion   .= energy_MeV_in
        elseif k == 2
            photon_flux_synch .= exp10.(photon_flux_in)
            energy_MeV_synch  .= energy_MeV_in
        elseif k == 3
            photon_flux_IC    .= exp10.(photon_flux_in)
            energy_MeV_IC     .= energy_MeV_in
        end

        close(grid_file_unit)
    end
    #-------------------------------------------------------------------------
    # Emission data read into arrays


    # Convert pion and synchrotron emission from plasma frame to ISM frame.
    # Use basic Doppler shift formula, e.g. Longair, 3rd. ed., pg. 238
    #     E_new/E_old = γ (1 + β cos(θ))
    # WARNING: Current method assumes isotropic emission in both initial and final frame,
    # in that it strips away any angular dependence.
    #-------------------------------------------------------------------------
    # Determine number of emission types being processed. If electrons were included, kmax
    # is currently 3; however, since we don't need to deal with IC emission at this point
    # we'll skip it during this step.
    if kmax > 1
        kmax = 2
    end


    # Set number of angular bins to use, as well as the fractional area array. The
    # orientation is chosen to be consistent with particle momenta, so that boundary 0
    # corresponds to cos(θ) = -1, and is directed *upstream*. This choice of orientation
    # means the + sign in the Doppler shift formula must be changed to a - sign, as will be
    # noted again later.
    n_cos_bins = 180
    cosθ = range(-1.0, 1.0, length=n_cos_bins+1)
    # Fractional surface area on a sphere, i.e. ∫ sin(θ) dθ dφ, divided by 4π steradians
    frac_area = 1/n_cos_bins

    photon_flux_out = Matrix{Float64}(undef, (na_photons, n_grid))
    # Loop over the one (or two) plasma frame spectra to be Doppler shifted into the ISM frame.
    for k in 1:kmax

        if k == 1
            n_photon = n_photon_pion - 1
            energy_MeV_in[1:n_photon] .= exp10(energy_MeV_pion[1:n_photon])
            photon_flux_in  .= photon_flux_pion
            fill!(photon_flux_out, 1e-99)
        elseif k == 2
            n_photon = n_photon_synch - 1
            energy_MeV_in[1:n_photon] .= exp10(energy_MeV_synch[1:n_photon])
            photon_flux_in  .= photon_flux_synch
            fill!(photon_flux_out, 1e-99)
        end

        # Obtain energy widths of each bin now, since they're constant across all the grid zones.
        # Remember that photon_flux_in is number of photons per log energy, so multiplying by
        # the difference of the logarithms gives the area under the curve, or number of photons
        # in this bin. Also obtain the midpoint (in log space) of each energy bin.
        ΔlogE = energy_MeV_in[3] / energy_MeV_in[2] # position chosen arbitrarily
        for j in 1:n_photon-1
            energy_mid[j] = √(energy_MeV_in[j] * energy_MeV_in[j+1])
        end
        # Highest-energy zone handled separately, but correctly.
        energy_mid[n_photon] = energy_mid[n_photon-1] * exp10(ΔlogE)

        # Loop over grid positions with emission data
        for i in i_grid_start[k]:i_grid_end[k]

            # Obtain the current flow speed difference between the plasma and ISM frames
            γ_sk_pf  = γ_sf_grid[i]
            β_sk_pf  = utot_grid[i] / c_cgs
            β_pf_ISM = β_ef_grid[i]
            γ_pf_ISM = γ_ef_grid[i]

            # Convert photon fluxes from number per log(energy) to pure number
            for j in 1:n_photon
                if photon_flux_in[j,i] > 1e-90
                    photon_flux_in[j,i] *= ΔlogE
                end
            end

            # Assume emission was isotropic, so that emission in each slice of a sphere is
            # proportional to the fractional area of that slice. Now use the Doppler shift
            # formula to transform the energy of the *central* ray in each energy/angle bin,
            # which we will use to re-bin in the ISM frame.
            for l in 1:n_cos_bins, j in 1:n_photon

                xnum_trans = photon_flux_in[j,i] * frac_area

                # With the energy of the transformed photons, as well as the number of
                # photons coming from this angle, re-bin in terms of energy
                if xnum_trans > 1e-99

                    # Determine the two energy ratios at the angular boundaries using the
                    # Doppler shift formula. Recall that since a cosine of -1 points upstream,
                    # i.e. towards the observer, the + sign in the Doppler shift formula must
                    # be changed to a - sign.
                    dimless_avg = √((1 - β_pf_ISM*cosθ[l]) * (1 - β_pf_ISM*cosθ[l+1]))
                    energy_trans = energy_mid[j] * γ_pf_ISM * dimless_avg

                    m = findnext(≥(energy_trans), energy_MeV_in, 2)
                    photon_flux_out[m-1,i] += xnum_trans
                end

            end # loop over incoming photon energies and angles

        end # loop over grid positions


        # Now that we have the transformed emission spectra, put them back into the
        # emission-specific arrays. Because these arrays were shifted from plasma to shock
        # frame, three factors of γ_pf_ISM need to be included: two from relativistic
        # beaming, and one from time dilation.
        mask = photon_flux_out .> 1e-95
        if k == 1
            photon_flux_pion[  mask]  .= photon_flux_out[mask] * γ_pf_ISM^3
            photon_flux_pion[.!mask]  .= 1e-99
        elseif k == 2
            photon_flux_synch[  mask] .= photon_flux_out[mask] * γ_pf_ISM^3
            photon_flux_synch[.!mask] .= 1e-99
        end

    end # loop over emission types
    #-------------------------------------------------------------------------
    # Pion and synchrotron emission converted from plasma frame to ISM frame


    # Loop over all emission types and sum them into spectral regions, still separated by emission process.
    #--------------------------------------------------------------------------
    # Determine number of emission types being processed. If electrons were
    # included, kmax is currently 2 due to the previous section of code
    if kmax > 1
        kmax = 3
    end


    # Now loop over emission types to get the summed emission for each type in each spectral region
    for k in 1:kmax

        if k == 1
            n_photon = n_photon_pion - 1
            photon_flux_in .= photon_flux_pion
        elseif k == 2
            n_photon = n_photon_synch - 1
            photon_flux_in .= photon_flux_synch
        elseif k == 3
            n_photon = n_photon_IC - 1
            photon_flux_in .= photon_flux_IC
        end

        # Loop over spectral regions; summation will be quite easy because we
        # previously converted photon fluxes from log space to linear space.
        for n in 1:n_shells

            nn_lo = n_shell_endpoints[n]
            nn_hi = n_shell_endpoints[n+1]-1

            # Sum over second dimension in photon_flux_in, from nn_lo to nn_hi
            photon_flux_out[1:n_photon,n] .= sum(photon_flux_in[1:n_photon, nn_lo:nn_hi], dims=2)

            # Quick check to eliminate sums of "zeros"
            for j in 1:n_photon
                if photon_flux_out[j,n] < 1e-95
                    photon_flux_out[j,n] = 1e-99
                end
            end
        end

        # Fill the correct emission array with the summed emission; fill the rest with zero emission
        if k == 1
            photon_flux_pion .= photon_flux_out
            photon_flux_pion[:,n_shells+1:end] .= 1e-99
        elseif k == 2
            photon_flux_synch .= photon_flux_out
            photon_flux_synch[:,n_shells+1:end] .= 1e-99
        elseif k == 3
            photon_flux_IC .= photon_flux_out
            photon_flux_IC[:,n_shells+1:end] .= 1e-99
        end

    end
    #-------------------------------------------------------------------------
    # Summed spectra in regions found


    # Create a single emission spectrum accounting for all emission processes. During this
    # process, return all photon flux arrays to log space and divide by d(logE) to restore
    # the units of the original spectra as calculated in the photon production subroutines.
    #-------------------------------------------------------------------------
    # Zero out array that will hold totalled spectra, set total number of points in those spectra,
    # and the energy range of each bin (the same for all three emission processes by design)
    n_pts_sum = trunc(Int, (log10(photon_energy_max_MeV) - log10(photon_energy_min_MeV)) * n_bins_per_decade )
    ΔlogE = 1 / n_bins_per_decade

    # All three emission spectra use logarithmic energy scales with the same number of bins
    # per decade in energy. So determine where the first bin of each spectrum falls in the
    # summed spectrum
    n_pion_start  = trunc(Int, (log10(photon_pion_min_MeV)  - log10(photon_energy_min_MeV)) * n_bins_per_decade)
    n_synch_start = trunc(Int, (log10(photon_synch_min_MeV) - log10(photon_energy_min_MeV)) * n_bins_per_decade)
    n_ic_start    = trunc(Int, (log10(photon_ic_min_MeV)    - log10(photon_energy_min_MeV)) * n_bins_per_decade)


    # Loop over spectral shells, creating a total spectrum for each
    fill!(photon_flux_out, 1e-99) # "zero" out array
    for n in 1:n_shells

        # Add pion spectrum to total
        n_photon = n_photon_pion - 1
        n_start  = n_pion_start

        for j in 1:n_photon
            photon_flux_out[n_start+j,n] = photon_flux_pion[j,n]

            if photon_flux_pion[j,n] > 1e-99
                photon_flux_pion[j,n] = log10(photon_flux_pion[j,n] / ΔlogE)
            else
                photon_flux_pion[j,n] = -99.0
            end
        end


        # If electrons were included, add synchrotron and IC spectra as well
        if kmax > 1
            n_photon = n_photon_synch - 1
            n_start  = n_synch_start

            for j in 1:n_photon
                photon_flux_out[n_start+j,n] += photon_flux_synch[j,n]

                if photon_flux_synch[j,n] > 1e-99
                    photon_flux_synch[j,n] = log10(photon_flux_synch[j,n] / ΔlogE)
                else
                    photon_flux_synch[j,n] = -99.0
                end
            end

            n_photon = n_photon_IC - 1
            n_start  = n_ic_start

            for j in 1:n_photon
                photon_flux_out[n_start+j,n] += photon_flux_IC[j,n]

                if photon_flux_IC[j,n] > 1e-99
                    photon_flux_IC[j,n] = log10(photon_flux_IC[j,n] / ΔlogE)
                else
                    photon_flux_IC[j,n] = -99.0
                end
            end
        end

    end  # loop over shells


    # Sum all spectral shells into a single, total, spectrum
    photon_flux_tot = sum(photon_flux_out[1:n_pts_sum,1:n_shells], dims=2)

    # Return spectral regions to log space, flux per d(logE).
    for j in 1:n_pts_sum
        for n in 1:n_shells
            if photon_flux_out[j,n] > 1e-96
                photon_flux_out[j,n] = log10(photon_flux_out[j,n] / ΔlogE)
            else
                photon_flux_out[j,n] = -99.0
            end
        end

        if photon_flux_tot[j] > 1e-96
            photon_flux_tot[j] = log10(photon_flux_tot[j] / ΔlogE)
        else
            photon_flux_tot[j] = -99.0
        end
    end
    #-------------------------------------------------------------------------
    # Completely summed spectrum calculated


    # Generate output for the spectral regions of the summed emission, each emission process,
    # and for the total spectrum. Because all fluxes are back in log space, generating
    # *energy* fluxes (rather than number fluxes) takes slightly different formula.
    #-------------------------------------------------------------------------
    energy_MeV_out = zeros(0:na_photons)
    energy_MeV_out[1:n_pts_sum] .= range(start  = log10(photon_energy_min_MeV),
                                         step   = ΔlogE,
                                         length = n_pts_sum)

    # Create histograms for summed emission
    create_histograms()

    # Now create the histograms for the various emission processes' spectral regions
    for k in 1:kmax

        if k == 1
            n_photon = n_photon_pion - 1

            # energy_MeV_*** in log space already
            energy_MeV_out[1:n_photon] .= energy_MeV_pion[1:n_photon]
            photon_flux_out[1:n_photon, 1:n_shells] .= photon_flux_pion[1:n_photon, 1:n_shells]

            open(newunit=j_unit, status="unknown", file="./photon_pion_summed.dat")

        elseif k == 2
            n_photon = n_photon_synch - 1

            # energy_MeV_*** in log space already
            energy_MeV_out[1:n_photon] .= energy_MeV_synch[1:n_photon]
            photon_flux_out[1:n_photon, 1:n_shells] .= photon_flux_synch[1:n_photon, 1:n_shells]

            open(newunit=j_unit, status="unknown", file="./photon_synch_summed.dat")

        elseif k == 3
            n_photon = n_photon_IC - 1

            # energy_MeV_*** in log space already
            energy_MeV_out[1:n_photon] .= energy_MeV_IC[1:n_photon]
            photon_flux_out[1:n_photon, 1:n_shells] .= photon_flux_IC[1:n_photon, 1:n_shells]

            open(newunit=j_unit, status="unknown", file="./photon_IC_summed.dat")

        end

        for n in 1:n_shells
            j_plot = 1

            for j in 1:n_photon

                energy_keV  = energy_MeV_out[j] + 3
                photon_flux = photon_flux_out[j,n]
                energy_MeV  = energy_MeV_out[j]
                if photon_flux > -98
                    energy_flux_MeV = photon_flux + energy_MeV
                    energy_flux_keV = photon_flux + energy_keV
                else
                    energy_flux_MeV = -99.0
                    energy_flux_keV = -99.0
                end

                write(j_unit, # photon_***_summed.dat
                      n, j_plot,
                      energy_keV,           # 1 Log10(keV)
                      photon_flux,          # 2 Log10(photons/(cm²⋅sec))
                      energy_MeV,           # 3 Log10(MeV)
                      energy_flux_MeV,      # 4 Log10(MeV/(cm²⋅sec))
                      energy_flux_keV)      # 5 Log10(keV/(cm²⋅sec))

                j_plot += 1
                energy_keV = energy_MeV_out[j+1] + 3
                energy_MeV = energy_MeV_out[j+1]

                if j != n_photon
                    write(j_unit, # photon_***_summed.dat
                          n, j_plot,
                          energy_keV,       # 1 Log10(keV)
                          photon_flux,      # 2 Log10(photons/(cm²⋅sec))
                          energy_MeV,       # 3 Log10(MeV)
                          energy_flux_MeV,  # 4 Log10(MeV/(cm²⋅sec))
                          energy_flux_keV)  # 5 Log10(keV/(cm²⋅sec))
                end

                j_plot += 1

            end # loop over photon energies

            print_plot_vals(j_unit)

        end # loop over spectral locations

        close(j_unit)

    end # loop over emission types


    # Create histogram for total emission
    create_total_emission_histogram()

    #-------------------------------------------------------------------------
    # Output written

    #deallocate(photon_flux_in, photon_flux_out)
    #deallocate(photon_flux_pion, photon_flux_synch, photon_flux_IC)
end

function pion_emission()
    pion_photon_sum .= log10.(pion_photon_sum)

    # Make sure we're writing to photon_pion_decay_grid.dat, and at the end of the file.
    grid_file_unit = open("./photon_pion_decay_grid.dat", "a")

    # Loop over grid positions, skipping any where there was no photon
    # production due to pion decay
    i_grid_start[1] = n_grid + 1
    for i in 1:n_grid
        if count(pion_photon_sum[:,i] .> -90) ≥ 1

            # Find first/last position with photon production, since that task
            # would normally be handled during read-in from files
            i_grid_start[1] = min(i_grid_start[1], i)
            i_grid_end[1]   = i

            j_plot = 0
            for j in 1:n_photon_pion-1
                j_plot += 1

                energy_γ_MeV = log10(energy_pion_MeV[j])
                photon_flux  = pion_photon_sum[j,i]
                if photon_flux > -99
                    emis_γ_MeV = photon_flux + energy_γ_MeV
                else
                    emis_γ_MeV = -99.0
                end

                energy_γ_keV = energy_γ_MeV + 3
                if emis_γ_MeV > -96
                    emis_γ_keV = emis_γ_MeV + 3
                else
                    emis_γ_keV = -99.0
                end

                write(grid_file_unit, # photon_pion_decay_grid.dat
                      i, j_plot,
                      energy_γ_keV,             # 1 Log10(keV)
                      photon_flux,              # 2 Log10(photons/(cm²⋅sec))
                      energy_γ_MeV,             # 3 Log10(MeV)
                      emis_γ_MeV,               # 4 Log10[MeV/(cm²⋅sec)] at Earth
                      emis_γ_keV,               # 5 Log10[keV/(cm²⋅sec)] at Earth
                      photon_flux-energy_γ_keV) # 6 Log10[photons/(cm²⋅sec-keV)]
            end
            print_plot_vals(grid_file_unit)

        end
    end

    close(grid_file_unit)
end

function inverse_compton_emission()
    ic_photon_sum .= log10(ic_photon_sum[:,:])

    # Make sure we're writing to photon_IC_grid.dat, and at the end of the file.
    # Note that since IC emission is separated by field, this file should not exist
    # before this part of the code.
    grid_file_unit = open("./photon_IC_grid.dat", "a")

    # Loop over grid positions, skipping any where there was no photon
    # production due to pion decay
    i_grid_start[3] = n_grid + 1
    for i in 1:n_grid
        if count(ic_photon_sum[:,i] > -90) ≥ 1

            # Find first/last position with photon production, since that task would
            # normally be handled during read-in from files
            i_grid_start[3] = min(i_grid_start[3], i)
            i_grid_end[3]   = i

            j_plot = 0
            for j in 1:n_photon_IC-1
                j_plot += 1

                energy_γ_MeV = log10(energy_IC_MeV[j])
                photon_flux  = ic_photon_sum[j,i]
                if photon_flux > -99
                    emis_γ_MeV = photon_flux + energy_γ_MeV
                else
                    emis_γ_MeV = -99.0
                end

                energy_γ_keV = energy_γ_MeV + 3
                if emis_γ_MeV > -96
                    emis_γ_keV = emis_γ_MeV + 3
                else
                    emis_γ_keV = -99.0
                end

                write(grid_file_unit, # photon_IC_grid.dat
                      i, j_plot,
                      energy_γ_keV,             # 1 Log10(keV)
                      photon_flux,              # 2 Log10(photons/(cm²⋅sec))
                      energy_γ_MeV,             # 3 Log10(MeV)
                      emis_γ_MeV,               # 4 Log10[MeV/(cm²⋅sec)] at Earth
                      emis_γ_keV,               # 5 Log10[keV/(cm²⋅sec)] at Earth
                      photon_flux-energy_γ_keV) # 6 Log10[photons/(cm²⋅sec-keV)]
            end
            print_plot_vals(grid_file_unit)

        end
    end

    close(grid_file_unit)
end

function create_histograms()
    open(newunit=j_unit, status="unknown", file="./photon_tot_summed.dat")
    for n in 1:n_shells
        j_plot = 1

        for j in 1:n_pts_sum

            energy_keV  = energy_MeV_out[j] + 3
            photon_flux = photon_flux_out[j,n]
            energy_MeV  = energy_MeV_out[j]
            if photon_flux > -98
                energy_flux_MeV = photon_flux + energy_MeV
                energy_flux_keV = photon_flux + energy_keV
            else
                energy_flux_MeV = -99.0
                energy_flux_keV = -99.0
            end

            write(j_unit,             # photon_tot_summed.dat
                  n, j_plot,
                  energy_keV,           # 1 Log10(keV)
                  photon_flux,          # 2 Log10(photons/(cm²⋅sec))
                  energy_MeV,           # 3 Log10(MeV)
                  energy_flux_MeV,      # 4 Log10(MeV/(cm²⋅sec))
                  energy_flux_keV)      # 5 Log10(keV/(cm²⋅sec))

            j_plot += 1
            energy_keV = energy_MeV_out[j+1] + 3
            energy_MeV = energy_MeV_out[j+1]

            if j != n_pts_sum
                write(j_unit,         # photon_tot_summed.dat
                      n, j_plot,
                      energy_keV,       # 1 Log10(keV)
                      photon_flux,      # 2 Log10(photons/(cm²⋅sec))
                      energy_MeV,       # 3 Log10(MeV)
                      energy_flux_MeV,  # 4 Log10(MeV/(cm²⋅sec))
                      energy_flux_keV)  # 5 Log10(keV/(cm²⋅sec))
            end

            j_plot += 1

        end # loop over photon energies

        print_plot_vals(j_unit)

    end # loop over spectral locations
    close(j_unit)
end

function create_total_emission_histogram()
    energy_MeV_out[1:n_pts_sum] .= range(start = log10(photon_energy_min_MeV),
                                         step = ΔlogE,
                                         length = n_pts_sum)

    open(newunit=j_unit, status="unknown", file="./photon_tot.dat")

    j_plot = 1

    for j in 1:n_pts_sum

        energy_keV  = energy_MeV_out[j] + 3
        photon_flux = photon_flux_tot[j]
        energy_MeV  = energy_MeV_out[j]
        if photon_flux > -98.0
            energy_flux_MeV = photon_flux + energy_MeV
            energy_flux_keV = photon_flux + energy_keV
        else
            energy_flux_MeV = -99.0
            energy_flux_keV = -99.0
        end

        write(j_unit,   # photon_tot.dat
              n, j_plot,
              energy_keV,                       # 1 Log10(keV)
              photon_flux,                      # 2 Log10(photons/(cm²⋅sec))
              energy_MeV,                       # 3 Log10(MeV)
              energy_flux_MeV,                  # 4 Log10(MeV/(cm²⋅sec))
              energy_flux_keV)                  # 5 Log10(keV/(cm²⋅sec))

        j_plot += 1
        energy_keV = energy_MeV_out[j+1] + 3
        energy_MeV = energy_MeV_out[j+1]

        if j != n_pts_sum
            write(j_unit, # photon_tot.dat
                  n, j_plot,
                  energy_keV,                   # 1 Log10(keV)
                  photon_flux,                  # 2 Log10(photons/(cm²⋅sec))
                  energy_MeV,                   # 3 Log10(MeV)
                  energy_flux_MeV,              # 4 Log10(MeV/(cm²⋅sec))
                  energy_flux_keV)              # 5 Log10(keV/(cm²⋅sec))
        end

        j_plot += 1

    end # loop over photon energies

    print_plot_vals(j_unit)
    close(j_unit)
end
