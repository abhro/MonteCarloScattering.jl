using .constants: erg2MeV, MeV2erg
using .parameters: psd_max, na_photons
using .io: print_plot_vals

"""
Calculates photon production by inverse Compton emission for a provided distribution of electrons.

While the subroutine can handle any number of photon distributions, and can readily be
changed to read in an arbitrary distribution (e.g. the photon fields calculated in
photon_pion_decay and photon_synch), right now the subroutine is set up to only work with
one photon field: the CMB.

### Arguments
- n_grid: current grid zone. Used for tracking emission output
- p_pf_cgs_cr: momentum boundary values, cgs units, of cosmic ray distribution histogram
- num_psd_mom_bins: number of momentum bins in the distribution of accelerated particles
- cos_bounds: boundaries of the angular bins of PSD, but in cosine form rather than mixed
  cosine and θ form (as in psd_θ_bounds)
- num_psd_θ_bins: number of angular bins in the distribution array
- d2Ndp_slice: particle distribution array. It is number of particles per dp (NOT #/cm^3/dp),
  split by pitch angle (NOT per pitch angle)
- n_photon_IC: number of energy bins to use for photon production
- photon_ic_min_MeV: minimum photon energy, in MeV, to use for inverse Compton spectrum.
  Passed to IC_emission_FCJ.
- bins_per_dec_photon: number of energy bins per decade of photon spectrum. Passed to IC_emission_FCJ.
- dist_lum: luminosity distance (i.e. including redshift correction) to source; passed to IC_emission_FCJ.
- redshift: redshift of source, used to adjust photon energies; passed to IC_emission_FCJ.
"""
function photon_IC(
        n_grid, p_pf_cgs_cr, num_psd_mom_bins, cos_bounds,
        num_psd_θ_bins, d2Ndp_slice, n_photon_IC, photon_ic_min_MeV,
        bins_per_dec_photon, dist_lum, redshift,
        n_IC_specs, energy_IC_MeV, ic_photon_sum,
    )

    #energy_γ_cgs = zeros(na_photons)
    #ic_emis = zeros(na_photons)
    energy_γ_MeV = zeros(na_photons)

    # Our distribution function has already been normalized to the total
    # number of emitting particules, so no additional scaling is needed.
    # However, it must be converted from particles per momentum into pure particle count.
    # Note the ordering of the dimensions in d2N_slice:
    #    1  -  angle          2  -  momentum
    d2N_slice = zeros(0:psd_max, 0:psd_max)
    for i in 0:num_psd_mom_bins, j in 0:num_psd_θ_bins
        if d2Ndp_slice[j,i] ≤ 1e-99
            d2N_slice[j,i] = 1e-99
        else
            d2N_slice[j,i] = d2Ndp_slice[j,i] * ( p_pf_cgs_cr[i+1] -  p_pf_cgs_cr[i])
        end
    end


    # Loop over the photon fields, creating an output file for each field's
    # contribution to IC emission.
    # Right now, there's just one photon field: the CMB.
    # TODO: have access to synchrotron photon field. Use it.
    # TODO: investigate interplay between j3 and n_IC_specs if more than
    #  one field is used. not sure current code is right
    for j3 in 1:1# i_photon_fields

        # Set number of inverse Compton spectra here. Will be 1 unless more than one photon
        # field was used. Also, zero out the column of ic_photon_sum that would hold the
        # summed emission if more than one photon field is used.
        n_IC_specs = j3
        if j3 == 1
            ic_photon_sum[:,n_grid] = 1e-99
            energy_IC_MeV .= 1e-99
        end

        # Compute the spectra
        energy_γ_cgs, ic_emis = IC_emission_FCJ(
            num_psd_mom_bins, p_pf_cgs_cr,
            num_psd_θ_bins, cos_bounds, d2N_slice, n_photon_IC, j3,
            photon_ic_min_MeV, bins_per_dec_photon, dist_lum, redshift,
            jet_sph_frac, o_o_mc)

        # Convert units of energy_γ_cgs and ic_emis
        # NOTE: unlike the other two photon production subroutines, ic_emis comes out of
        # IC_emission_FCJ already in the form of observed energy flux at Earth per
        # log energy bin, i.e. dP/(d(lnE)-dA). Pion production and synchrotron both require
        # processing from dP/d(lnE). The units of ic_emis are [erg/sec-cm^2].
        for i in 1:n_photon_IC
            energy_γ_MeV[i] = energy_γ_cgs[i]*erg2MeV

            # Add the current photon flux to the running total over all fields
            if ic_emis[i] > 1e-99
                ic_photon_sum[i,n_grid] += ic_emis[i] / energy_γ_cgs[i]
            end
            energy_IC_MeV[i] = energy_γ_MeV[i]

        end   # units on emis_γ: erg/(cm^2-sec)

        # Don't write out anything if ic_emis is empty
        count(ic_emis[1:n_photon_IC] .> 1e-99) < 1 && continue

        # Do necessary unit conversion and write out results to file
        iplot = 0
        for i in 1:n_photon_IC
            iplot += 1

            # Above is energy flux [MeV/(cm^2-sec)] per log energy bin d(lnE) = dE/E.
            # Energy in spectrum is area under curve when plotted with a
            # logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if ic_emis[i] > 1e-99
                emis_γ_MeV = ic_emis[i]/MeV2erg  # MeV/(cm^2-sec) at earth
            else
                emis_γ_MeV = 1e-99
            end


            if emis_γ_MeV ≤ 1e-99
                emis_γ_keV = 1e-99               # "zero" emission
            else
                emis_γ_keV = emis_γ_MeV * 1e3 # keV/(cm^2-sec) at earth
            end

            xMeV_log   = log10(energy_γ_MeV[i])
            energy_keV = energy_γ_MeV[i]*1000
            xkeV_log   = xMeV_log + 3

            # This is photon flux [#/(cm^2-sec)] per log energy bin d(lnE) = dE/E.
            # Number of photons in spectrum is area under curve when plotted
            # with a logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if emis_γ_MeV ≤ 1e-99
                photon_flux = 1e-99                  # "zero" emission
            else
                photon_flux = emis_γ_MeV/energy_γ_MeV[i]
            end

            # Open the file to which we will write the spectral data.
            if j3 == 1
                inquire(file="./photon_IC_grid.dat", opened=lopen)
                lopen || open(newunit=j_unit, status="unknown", file="./photon_IC_grid.dat")
            end

            # Don't write out the last data point
            i == n_photon_IC && continue

            ##TODO: incorporate redshift into this writeout;
            # right now, everything is handled in time_seq_photons so this section is irrelevant
            write(j_unit, n_grid, iplot, # photon_IC_grid.dat
                  j3,                           # 1 photon source
                  (xMeV_log + 3),               # 2 Log10(keV)
                  log10(photon_flux),           # 3 Log10(photons/(cm^2-sec))
                  xMeV_log,                     # 4 Log10(MeV)
                  log10(emis_γ_MeV),            # 5 Log[MeV/(cm^2-sec)] at earth
                  log10(emis_γ_keV),            # 6 Log[keV/(cm^2-sec)] at earth
                  log10(photon_flux/energy_keV))# 7 Log10[photons/(cm^2-sec-keV)]

        end # loop over n_photon_IC

        print_plot_vals(j_unit)

    end # Loop over photon fields

    close(j_unit)

    #deallocate(d2N_slice)
end
