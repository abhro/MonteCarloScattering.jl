using Unitful, UnitfulAstro
using Unitful: cm
using UnitfulAstro: kpc

using .parameters: psd_max, na_photons

# To make the resultant spectra easier to add, set min and max energies
# for the overall spectrum and the individual emission spectra.
const E_photon_min = 1e-13MeV
const E_photon_max = 1e+12MeV
const bins_per_dec_photon   = 10

const E_photon_pion_min  = 1.0MeV
const E_photon_pion_max  = E_photon_max
const E_photon_synch_min = E_photon_min
const E_photon_synch_max = 1e5MeV
const E_photon_ic_min    = 1e-2MeV
const E_photon_ic_max    = E_photon_max

"""
    photon_calcs(...)

### Arguments
TODO
"""
function photon_calcs(
        dNdp_therm_pvals, dNdp_therm, dNdp_cr, aa,
        n_shell_endpoints, d²N_dpdcos_ef,
        jet_dist_kpc, redshift, psd_lin_cos_bins, n_ions,
        num_upstream_shells, num_downstream_shells,
        psd_mom_bounds, psd_θ_bounds, num_psd_mom_bins, num_psd_θ_bins,
        i_ion,
    )

    # Before doing any photon calculations, calculate the luminosity distance
    # (which includes the correction due to redshift)
    dist_lum = ustrip(cm, jet_dist_kpc*(1 + redshift) * kpc)

    # Check to make sure that the arrays can store the emission spectra
    if log10(E_photon_max / E_photon_min) * bins_per_dec_photon > na_photons
        error("Photon arrays will not be big enough. Check na_photons, compare against photon_energy_***_MeV.")
    end

    # Set the number of bins used in the individual emission spectra
    n_photon_pion  = trunc(Int, log10(E_photon_pion_max  / E_photon_pion_min)  * bins_per_dec_photon)
    n_photon_synch = trunc(Int, log10(E_photon_synch_max / E_photon_synch_min) * bins_per_dec_photon)
    n_photon_IC    = trunc(Int, log10(E_photon_ic_max    / E_photon_ic_min)    * bins_per_dec_photon)

    n_shells = num_upstream_shells + num_downstream_shells

    num_hist_bins = psd_max
    p_pf_therm = zeros(0:num_hist_bins)
    dNdp_pf_therm  = zeros(0:num_hist_bins)
    p_pf_cr    = zeros(0:num_psd_mom_bins+1)
    dNdp_pf_cr     = zeros(0:num_psd_mom_bins+1)
    cos_bounds     = zeros(0:num_psd_θ_bins+1)
    d²Ndp_slice    = zeros(0:psd_max, 0:psd_max)
    # Main loop over photon emission shells and zones in this spectral region
    for i in 1:n_shells,
        n in n_shell_endpoints[i]:n_shell_endpoints[i+1]-1

        if aa ≥ 1 # Pion emission from nuclei

            # Take info from dN/dp and prepare it for photon_pion_decay
            # For pion decay, need plasma frame distribution, which is index
            # 2 in the third dimension (1 is shock frame, 3 is ISM frame)
            p_pf_therm .= dNdp_therm_pvals[:,n,2]
            dNdp_pf_therm  .= dNdp_therm[:,n,2]

            # For pion decay, need plasma frame distribution, which is index
            # 2 in the third dimension (1 is shock frame, 3 is ISM frame)
            p_pf_cr .= exp10.(psd_mom_bounds)
            dNdp_pf_cr  .= dNdp_cr[:,n,2]

            # Don't bother calculating photon production unless at least one
            # particle was present in this grid zone
            if (count(>(1e-99), dNdp_pf_therm) + count(>(1e-99), dNdp_pf_cr)) ≥ 1
                photon_pion_decay(n, num_hist_bins, p_pf_therm,
                                  dNdp_pf_therm, num_psd_mom_bins, p_pf_cr, dNdp_pf_cr,
                                  n_photon_pion, E_photon_pion_min, bins_per_dec_photon, dist_lum, redshift)
            end

        else # Synchrotron and IC emission from electrons

            # Take info from dN/dp and prepare it for photon_synch
            num_hist_bins = psd_max  # dNdp_**** ranges from 0 to psd_max
            # For synchrotron, need plasma frame distribution, which is index
            # 2 in the third dimension (1 is shock frame, 3 is ISM frame)
            p_pf_therm .= dNdp_therm_pvals[:,n,2]
            dNdp_pf_therm  .= dNdp_therm[:,n,2]
            # For synchrotron, need plasma frame distribution, which is index
            # 2 in the third dimension (1 is shock frame, 3 is ISM frame)
            p_pf_cr .= exp10.(psd_mom_bounds)
            dNdp_pf_cr  .= dNdp_cr[:,n,2]

            # Don't bother calculating photon production unless at least one
            # particle was present in this grid zone
            if (count(>(1e-99), dNdp_pf_therm) + count(>(1e-99), dNdp_pf_cr)) ≥ 1
                photon_synch(n, num_hist_bins, p_pf_therm, dNdp_pf_therm,
                             num_psd_mom_bins, p_pf_cr, dNdp_pf_cr, n_photon_synch,
                             E_photon_synch_min, bins_per_dec_photon, dist_lum, redshift)
            end

            # Take info from d²N_dpdcos_ef and prepare it for photon_IC
            p_pf_cr .= exp10.(psd_mom_bounds)
            for j in 0:num_psd_θ_bins+1
                # Determine current cosine, remembering that psd_θ_bounds has both a
                # linearly-spaced region in cosine and a logarithmically-spaced region in θ.
                # Also need to remember that the most finely spaced bins should occur in the
                # upstream-pointing direction, so need to negate psd_θ_bounds to get true cosine value.
                if j > (num_psd_θ_bins - psd_lin_cos_bins)
                    cos_bounds[j] = -psd_θ_bounds[j]
                else
                    cos_bounds[j] = -cos(psd_θ_bounds[j])
                end
            end
            d²Ndp_slice .= d²N_dpdcos_ef[:,:,n]

            # Don't bother calculating photon production unless at least one
            # particle was present in this grid zone
            if count(d²Ndp_slice .> 1e-99) ≥ 1
                photon_IC(n, p_pf_cr, num_psd_mom_bins, cos_bounds, num_psd_θ_bins, d²Ndp_slice,
                          n_photon_IC, E_photon_ic_min, bins_per_dec_photon, dist_lum, redshift)
            end

        end  # check on aa
    end
    #-------------------------------------------------------------------------
    # Loop over grid zones in this shell finished and loop over photon shells complete


    # Sum emission from the shells, and do additional processing to get fluxes that would be
    # observed at Earth.
    if i_ion == n_ions
        # The fluxes that come out of the emission subroutines are specific to the
        # individual zones of emission. They are also in different reference frames (IC was
        # calculated in ISM frame, everything else in plasma). They must be processed and
        # combined into emission viewed in a single reference frame.
        get_summed_emission(n_shells, n_shell_endpoints, n_photon_pion, n_photon_synch, n_photon_IC,
                            E_photon_min, E_photon_max,
                            bins_per_dec_photon, E_photon_pion_min,
                            E_photon_synch_min, E_photon_ic_min)
    end
end
