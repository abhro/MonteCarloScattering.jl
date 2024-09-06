using Unitful, UnitfulAstro
using .constants: h_cgs#, ħ_cgs

"""
Calculates photon production by the particle distribution due to pion decay emission process.

TODO: Incorporate possibility that target particles aren't at rest

### Arguments
- n_grid: current grid zone. Used to get plasma-frame density and for
  tracking emission output
- num_hist_bins: number of momentum bins in the distribution of thermal particles
- p_pf_cgs_therm: momentum boundary values, cgs units, of thermal distribution histogram
- dNdp_pf_therm: thermal particle distribution. Calculated at end of each ion species, it is
  number of particles per dp (NOT #/cm³/dp)
- num_psd_mom_bins: number of momentum bins in the distribution of accelerated particles
- p_pf_cgs_cr: momentum boundary values, cgs units, of cosmic ray distribution histogram
- dNdp_pf_cr: cosmic ray distribution. Calculated at end of each ion species, it is
  number of particles per dp (NOT #/cm³/dp)
- n_photon_pion: number of energy bins to use for photon production
- photon_pion_min_MeV: minimum photon energy, in MeV, to use for pion decay spectrum.
  Passed to pion_karlsson().
- bins_per_dec_photon: number of energy bins per decade of photon spectrum.
  Passed to pion_karlsson().
- dist_lum: luminosity distance (i.e. including redshift correction) to source
- redshift: redshift of source, used to adjust photon energies

### Returns
TODO

### Modifies
TODO
- energy_pion_MeV
- pion_photon_sum
"""
function photon_pion_decay(
        n_grid, num_hist_bins, p_pf_cgs_therm,
        dNdp_pf_therm, num_psd_mom_bins, p_pf_cgs_cr, dNdp_pf_cr,
        n_photon_pion, photon_pion_min_MeV, bins_per_dec_photon, dist_lum, redshift,
        aa_ion, ρ_N₀_ion, n_ions, β₀, γ₀,
        # from grid_vars
        γ_sf_grid,
        # from species_vars
        i_ion,
        # from photons
        energy_pion_MeV, pion_photon_sum,
    )

    energy_γ_cgs = zeros(n_photon_pion)
    pion_emis = zeros(n_photon_pion)
    energy_γ_MeV = zeros(n_photon_pion)
    emis_γ = zeros(n_photon_pion)
    dN_therm = zeros(0:num_hist_bins+1)
    dN_cr = zeros(0:num_psd_mom_bins+1)


    # Set a handful of constants
    γ_β_loc = √(γ_sf_grid[n_grid]^2 - 1)
    target_density = ρ_N₀_ion[1] * (γ₀*β₀)/γ_β_loc

    aa = aa_ion[i_ion]

    # Set number of pion decay spectra. Will be one for each nucleus type.
    # Also, zero out pion_photon_sum.
    n_pion_specs = count(aa_ion[1:n_ions] .≥ 1)
    if i_ion == 1
        pion_photon_sum .= 1e-99
        energy_pion_MeV .= 1e-99
    end

    # Our distribution functions have already been normalized to the total number of
    # emitting particles, so no additional scaling is needed. However, they must be
    # converted from particles per momentum into pure particle count
    for i in 0:num_hist_bins-1
        if dNdp_pf_therm[i] ≤ 1e-99
            dN_therm[i] = 1e-99
        else
            dN_therm[i] = dNdp_pf_therm[i] * (p_pf_cgs_therm[i+1] - p_pf_cgs_therm[i])
        end
    end
    for i in 0:num_psd_mom_bins
        if dNdp_pf_cr[i] ≤ 1e-99
            dN_cr[i] = 1e-99
        else
            dN_cr[i] = dNdp_pf_cr[i] * (p_pf_cgs_cr[i+1] - p_pf_cgs_cr[i])
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
        pion_kafexhiu(num_hist_bins, p_pf_cgs_therm, dN_therm,
                      num_psd_mom_bins, p_pf_cgs_cr, dN_cr, n_photon_pion, energy_γ_cgs,
                      pion_emis, target_density, ID, photon_pion_min_MeV, bins_per_dec_photon, aa)


        # Convert units of energy_γ_cgs and pion_emis
        # Note that pion_emis is energy radiated per second per logarithmic
        # energy bin, i.e. dP/d(lnE). Its units are [erg/sec].
        for i in 1:n_photon_pion
            energy_γ_MeV[i] = ustrip(u"MeV", energy_γ_cgs[i]*u"erg")
            emis_γ[i]       = pion_emis[i] / (4π*dist_lum^2)
            if emis_γ[i] < 1e-99
                emis_γ[i]  = 1e-99
            else
                # Add the current photon flux to the running total over all species
                pion_photon_sum[i,n_grid] += emis_γ[i] / energy_γ_cgs[i]
            end
            energy_pion_MeV[i] = energy_γ_MeV[i]
        end


        # Don't write out anything if emis_γ is empty
        count(emis_γ .> 1e-99) < 1 && continue


        # Do necessary unit conversions and write out results to file
        for i in 1:n_photon_pion

            # This is energy flux [MeV/(cm²⋅sec)] per log energy bin d(lnE) = dE/E.
            # Energy in spectrum is area under curve when plotted with a logarithmic
            # energy axis, i.e. dΦ/d(lnE) ⋅ d(lnE).
            if emis_γ[i] > 1e-99
                emis_γ_MeV = ustrip(u"MeV", emis_γ[i]*u"erg")  # MeV/(cm²⋅sec) at earth
            else
                emis_γ_MeV = 1e-99
            end

            if emis_γ_MeV ≤ 1e-99
                emis_γ_keV = 1e-99            # "zero" emission
            else
                emis_γ_keV = emis_γ_MeV * 1e3 # keV/(cm²⋅sec) at earth
            end

            ν_γ        = energy_γ_cgs[i]/h_cgs # frequency (ν)
            #ω_γ       = energy_γ_cgs[i]/ħ_cgs # angular frequency (ω)
            #f_jansky  = max(1e-99, ustrip(u"Jy", emis_γ[i]/ν_γ*u"erg/cm^2"))
            xMeV_log   = log10(energy_γ_MeV[i])
            energy_keV = energy_γ_MeV[i]*1000
            #xkeV_log  = xMeV_log + 3

            # This is photon flux [#/(cm²⋅sec)] per log energy bin d(lnE) = dE/E.
            # Number of photons in spectrum is area under curve when plotted with
            # a logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if emis_γ_MeV ≤ 1e-99
                photon_flux = 1e-99                  # "zero" emission
            else
                photon_flux = emis_γ_MeV/energy_γ_MeV[i]
            end

            # Don't write out the last data point
            i == n_photon_pion && continue

            ##TODO: incorporate redshift into this writeout; right now,
            # everything is handled in time_seq_photons so this section is irrelevant
            write(j_unit, # photon_pion_decay_grid.dat
                  n_grid, i,
                  i_ion,                          # 1 nucleus species
                  (xMeV_log + 3),                 # 2 Log10(keV)
                  log10(photon_flux),             # 3 Log10(photons/(cm²⋅sec))
                  xMeV_log,                       # 4 Log10(MeV)
                  log10(emis_γ_MeV),              # 5 Log10[MeV/(cm²⋅sec)] at earth
                  log10(emis_γ_keV),              # 6 Log10[keV/(cm²⋅sec)] at earth
                  log10(photon_flux/energy_keV))  # 7 Log10[photons/(cm²⋅sec⋅keV)]

        end # Loop over photon energies

        print_plot_vals(j_unit)

    end # Loop over secondary particle types

    close(j_unit)

    return n_pion_specs
end
