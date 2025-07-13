using Unitful, UnitfulAstro
using .constants: h

"""
    photon_pion_decay(...)

Calculates photon production by the particle distribution due to pion decay emission process.

TODO: Incorporate possibility that target particles aren't at rest

### Arguments
- `n_grid`: current grid zone. Used to get plasma-frame density and for
  tracking emission output
- `num_hist_bins`: number of momentum bins in the distribution of thermal particles
- `p_pf_therm`: momentum boundary values, cgs units, of thermal distribution histogram
- `dNdp_pf_therm`: thermal particle distribution. Calculated at end of each ion species, it is
  number of particles per dp (NOT #/cm³/dp)
- `num_psd_mom_bins`: number of momentum bins in the distribution of accelerated particles
- `p_pf_cr`: momentum boundary values, cgs units, of cosmic ray distribution histogram
- `dNdp_pf_cr`: cosmic ray distribution. Calculated at end of each ion species, it is
  number of particles per dp (NOT #/cm³/dp)
- `n_photon_pion`: number of energy bins to use for photon production
- `photon_pion_min`: minimum photon energy to use for pion decay spectrum.
  Passed to `pion_karlsson()`.
- `bins_per_dec_photon`: number of energy bins per decade of photon spectrum.
  Passed to pion_karlsson().
- `dist_lum`: luminosity distance (i.e. including redshift correction) to source
- `redshift`: redshift of source, used to adjust photon energies

### Returns
TODO

### Modifies
TODO
- `energy_pion`
- `pion_photon_sum`
"""
function photon_pion_decay(
        n_grid, num_hist_bins, p_pf_therm,
        dNdp_pf_therm, num_psd_mom_bins, p_pf_cr, dNdp_pf_cr,
        n_photon_pion, photon_pion_min, bins_per_dec_photon, dist_lum, redshift,
        aa_ion, n₀_ion, n_ions, β₀, γ₀,
        # from grid_vars
        γ_sf_grid,
        # from species_vars
        i_ion,
        # from photons
        energy_pion, pion_photon_sum,
    )

    pion_emis = zeros(n_photon_pion)
    energy_γ = zeros(typeof(1.0MeV), n_photon_pion)
    emis_γ = zeros(n_photon_pion)
    dN_therm = zeros(0:num_hist_bins+1)
    dN_cr = zeros(0:num_psd_mom_bins+1)


    # Set a handful of constants
    γ_β_loc = √(γ_sf_grid[n_grid]^2 - 1)
    target_density = n₀_ion[1] * (γ₀*β₀)/γ_β_loc

    aa = aa_ion[i_ion]

    # Set number of pion decay spectra. Will be one for each nucleus type.
    # Also, zero out pion_photon_sum.
    n_pion_specs = count(aa_ion[1:n_ions] .≥ 1)
    if i_ion == 1
        pion_photon_sum .= 1e-99
        energy_pion .= 1e-99MeV
    end

    # Our distribution functions have already been normalized to the total number of
    # emitting particles, so no additional scaling is needed. However, they must be
    # converted from particles per momentum into pure particle count
    for i in 0:num_hist_bins-1
        if dNdp_pf_therm[i] ≤ 1e-99
            dN_therm[i] = 1e-99
        else
            dN_therm[i] = dNdp_pf_therm[i] * (p_pf_therm[i+1] - p_pf_therm[i])
        end
    end
    for i in 0:num_psd_mom_bins
        if dNdp_pf_cr[i] ≤ 1e-99
            dN_cr[i] = 1e-99
        else
            dN_cr[i] = dNdp_pf_cr[i] * (p_pf_cr[i+1] - p_pf_cr[i])
        end
    end


    # Now loop over all secondary species: photons, e-/e+, (anti-)neutrinos.
    # For now, just limit to photons
    for ID in 1:1# 7

        # Open the files to which we will write our spectral data
        if ID == 1
            inquire(file="./photon_pion_decay_grid.dat", opened=lopen)
            if !lopen
                open(newunit=j_unit, status="unknown", position="append", file="./photon_pion_decay_grid.dat")
            end
        end


        # Generate both the energy bins for pion decay photon emission and the
        # amount of emission in those bins
        pion_kafexhiu(num_hist_bins, p_pf_therm, dN_therm,
                      num_psd_mom_bins, p_pf_cr, dN_cr, n_photon_pion, energy_γ,
                      pion_emis, target_density, ID, photon_pion_min, bins_per_dec_photon, aa)


        # Convert units of energy_γ and pion_emis
        # Note that pion_emis is energy radiated per second per logarithmic
        # energy bin, i.e. dP/d(lnE). Its units are [erg/sec].
        for i in 1:n_photon_pion
            emis_γ[i]   = pion_emis[i] / (4π*dist_lum^2)
            if emis_γ[i] < 1e-99
                emis_γ[i]  = 1e-99
            else
                # Add the current photon flux to the running total over all species
                pion_photon_sum[i,n_grid] += emis_γ[i] / energy_γ[i]
            end
            energy_pion[i] = energy_γ[i]
        end


        # Don't write out anything if emis_γ is empty
        count(emis_γ .> 1e-99) < 1 && continue


        # Do necessary unit conversions and write out results to file
        for i in 1:n_photon_pion

            # This is energy flux [MeV/(cm²⋅sec)] per log energy bin d(lnE) = dE/E.
            # Energy in spectrum is area under curve when plotted with a logarithmic
            # energy axis, i.e. dΦ/d(lnE) ⋅ d(lnE).
            if emis_γ[i] > 1e-99
                emis_γ = emis_γ[i]  # MeV/(cm²⋅sec) at earth
            else
                emis_γ = 1e-99MeV
            end

            ν_γ        = energy_γ[i]/h # frequency (ν)
            #ω_γ       = energy_γ[i]/ħ # angular frequency (ω)
            #f_jansky  = max(1e-99Jy, uconvert(Jy, emis_γ[i]/ν_γ*(erg/cm^2)))
            xMeV_log   = log10(ustrip(MeV, energy_γ[i]))

            # This is photon flux [#/(cm²⋅sec)] per log energy bin d(lnE) = dE/E.
            # Number of photons in spectrum is area under curve when plotted with
            # a logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if emis_γ ≤ 1e-99MeV
                photon_flux = 1e-99                  # "zero" emission
            else
                photon_flux = emis_γ/energy_γ[i]
            end

            # Don't write out the last data point
            i == n_photon_pion && continue

            ##TODO: incorporate redshift into this writeout; right now,
            # everything is handled in time_seq_photons so this section is irrelevant
            write(j_unit, # photon_pion_decay_grid.dat
                  n_grid, i,
                  i_ion,                          # 1 nucleus species
                  log10(photon_flux),             # 3 Log10(photons/(cm²⋅sec))
                  xMeV_log,                       # 3 Log10(MeV)
                  log10(emis_γ),                  # 4 Log10[MeV/(cm²⋅sec)] at earth
                  log10(photon_flux/energy_γ))    # 5 Log10[photons/(cm²⋅sec⋅MeV)]

        end # Loop over photon energies

        print_plot_vals(j_unit)

    end # Loop over secondary particle types

    close(j_unit)

    return n_pion_specs
end
