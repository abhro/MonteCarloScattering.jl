#using .constants: h_cgs, ħ_cgs
using .parameters: psd_max, na_photons
using .io: print_plot_vals

"""
Calculates photon production by an electron distribution due to synchtrotron emission.

### Arguments
- n_grid: current grid zone. Used to get plasma-frame density and for tracking emission output
- num_hist_bins: number of momentum bins in the distribution of thermal particles
- p_pf_cgs_therm: momentum boundary values, cgs units, of thermal distribution histogram
- dNdp_pf_therm: thermal particle distribution. Calculated at end of each ion species, it is
  number of particles per dp (NOT #/cm³/dp)
- num_psd_mom_bins: number of momentum bins in the distribution of accelerated particles
- p_pf_cgs_cr: momentum boundary values, cgs units, of cosmic ray distribution histogram
- dNdp_pf_cr: cosmic ray distribution. Calculated at end of each ion species, it is number
  of particles per dp (NOT #/cm³/dp)
- n_photon_synch: number of energy bins to use for photon production
- photon_synch_min_MeV: minimum photon energy, in MeV, to use for synchrotron spectrum.
  Passed to function synch_emission.
- bins_per_dec_photon: number of energy bins per decade of photon spectrum. Passed to
  function synch_emission.
- dist_lum: luminosity distance (i.e. including redshift correction) to source
- redshift: redshift of source, used to adjust photon energies
"""
function photon_synch(
        n_grid, num_hist_bins, p_pf_cgs_therm,
        dNdp_pf_therm, num_psd_mom_bins, p_pf_cgs_cr, dNdp_pf_cr, n_photon_synch,
        photon_synch_min_MeV, bins_per_dec_photon, dist_lum, redshift)

    # Our distribution function has already been normalized to the total number of emitting
    # particles, so no additional scaling is needed. However, it must be converted from
    # particles per momentum into pure particle count
    dN_therm = OffsetVector{Float64}(undef, 0:psd_max)
    for i in 0:num_hist_bins-1
        if dNdp_pf_therm[i] ≤ 1e-99
            dN_therm[i] = 1e-99
        else
            dN_therm[i] = dNdp_pf_therm[i] * (p_pf_cgs_therm[i+1] - p_pf_cgs_therm[i])
        end
    end
    dN_cr = OffsetVector{Float64}(undef, 0:psd_max)
    for i in 0:num_psd_mom_bins
        if dNdp_pf_cr[i] ≤ 1e-99
            dN_cr[i] = 1e-99
        else
            dN_cr[i] = dNdp_pf_cr[i] * (p_pf_cgs_cr[i+1] - p_pf_cgs_cr[i])
        end
    end


    # Generate both the energy bins for synchtrotron emission and the amount of emission
    energy_γ_cgs, synch_emis = synch_emission(
        n_grid, num_hist_bins, p_pf_cgs_therm, dN_therm,
        num_psd_mom_bins, p_pf_cgs_cr, dN_cr, n_photon_synch,
        photon_synch_min_MeV, bins_per_dec_photon,
        n_ions, aa_ion, ρ_N₀_ion, γ₀, u₀, flux_pₓ_UpS, flux_energy_UpS, u₂,
        n_grid, btot_grid,
        i_ion, mc)

    # Convert units of energy_γ_cgs and synch_emis
    # Note that synch_emis is energy radiated per second per logarithmic energy bin, i.e.,
    # dP/d(lnE). Its units are [erg/sec].
    energy_γ_MeV = ustrip(u"MeV", energy_γ_cgs*u"erg")
    emis_γ       = max.(synch_emis / (4π*dist_lum^2), 1e-99)

    # Don't write out anything if emis_γ is empty; different structure compared to pion
    # and IC subroutines because synchrotron subroutine doesn't have an internal loop.
    do_write = (count(emis_γ[1:n_photon_synch] > 1e-99) ≥ 1)


    # Do necessary unit conversions and write out results to file
    # TODO: maybe hold off on writing to file until IC has been done so that
    # we can take the SSC effect on the synchrotron spectrum into account
    iplot = 0
    if do_write

        for i in 1:n_photon_synch
            iplot += 1

            # This is energy flux [MeV/(cm²⋅sec)] per log energy bin d(lnE) = dE/E.
            # Energy in spectrum is area under curve when plotted with a
            # logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if emis_γ[i] > 1e-99
                emis_γ_MeV = ustrip(u"MeV", emis_γ[i]*u"erg")  # MeV/(cm²⋅s) at earth
            else
                emis_γ_MeV = 1e-99
            end

            if emis_γ_MeV ≤ 1e-99
                emis_γ_keV = 1e-99               # "zero" emission
            else
                emis_γ_keV = emis_γ_MeV * 1e3 # keV/(cm²⋅s) at earth
            end

            #ν_γ = energy_γ_cgs[i]/h_cgs # frequency (ν)
            #ω_γ = energy_γ_cgs[i]/ħ_cgs # ω

            #f_jansky = max(ustrip(u"Jy", emis_γ[i]/ν_γ * u"erg/cm^2"), 1e-99)

            xMeV_log   = log10(energy_γ_MeV[i])
            xkeV_log   = xMeV_log + 3
            energy_keV = energy_γ_MeV[i]*1000

            # This is photon flux [#/(cm²⋅s)] per log energy bin d(lnE) = dE/E.
            # Number of photons in spectrum is area under curve when plotted
            # with a logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            #                                  ↓ "zero" emission
            photon_flux = emis_γ_MeV ≤ 1e-99 ? 1e-99 : emis_γ_MeV/energy_γ_MeV[i]


            # Open the file and write the spectral data
            inquire(file="./photon_synch_grid.dat", opened=lopen)
            lopen || open(newunit=j_unit, status="unknown", file="./photon_synch_grid.dat")

            # Don't write out the last data point
            i == n_photon_synch && continue

            # TODO: incorporate redshift into this writeout; right now,
            # everything is handled in time_seq_photons so this section is irrelevant
            write(j_unit,           # photon_synch_grid.dat
                  n_grid, iplot,
                  xkeV_log,                     #1 Log10(keV)
                  log10(photon_flux),           #2 Log10(photons/(cm²⋅sec))
                  xMeV_log,                     #3 Log10(MeV)
                  log10(emis_γ_MeV),            #4 Log10[MeV/(cm²⋅sec)] at earth
                  log10(emis_γ_keV),            #5 Log10[keV/(cm²⋅sec)] at earth
                  log10(photon_flux/energy_keV))#6 Log10[photons/(cm²⋅sec⋅keV)]

        end # loop over n_photon_synch

        print_plot_vals(j_unit)

    end # check on do_write

    close(j_unit)
end
