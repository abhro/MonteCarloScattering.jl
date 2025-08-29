function ion_finalize(
    outfile, nc_unit,
    aa, esc_psd_feb_upstream, esc_psd_feb_downstream,
    jet_rad_pc, jet_sph_frac, m_ion, aa_ion, zz_ion, T₀_ion, n₀_ion,
    (u₀, β₀, γ₀), u₂,
    n_ions,
    do_multi_dNdps, do_photons,
    mc,
    n_grid, x_grid_cm, uₓ_sk_grid,
    i_iter,
    i_ion,
    γ_sf_grid,
    therm_grid, therm_pₓ_sk, therm_ptot_sk, therm_weight, num_crossings, n_cr_count,
    num_psd_mom_bins, psd_mom_bounds,
    psd, psd_lin_cos_bins, num_psd_θ_bins, psd_θ_bounds,
    psd_bins_per_dec_mom, psd_mom_min, psd_bins_per_dec_θ, psd_cos_fine, Δcos, psd_θ_min,
    n_shell_endpoints,
    flux_px_upstream, flux_energy_upstream, btot_grid,
    zone_vol, therm_energy_density, energy_density,
)
    println()
    println(outfile)

    # Obtain the dN/dp arrays for the current species: one 2-D array for each grid zone.
    dNdp_therm, dNdp_therm_pvals, dNdp_cr, zone_pop = get_normalized_dNdp(
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

    # Get pressure (both components) and kinetic energy density everywhere in the shock profile
    pressure_psd_par, pressure_psd_perp, energy_density_psd = thermo_calcs(
        num_crossings, n_cr_count, therm_grid, therm_pₓ_sk, therm_ptot_sk,
        therm_weight, nc_unit, psd, zone_pop, aa_ion, zz_ion, T₀_ion, n₀_ion,
        n_grid, γ_sf_grid, uₓ_sk_grid, i_ion, mc,
        num_psd_mom_bins, num_psd_θ_bins, psd_max, psd_θ_bounds, psd_mom_bounds,
        psd_bins_per_dec_mom, psd_mom_min, psd_bins_per_dec_θ, psd_cos_fine, Δcos, psd_θ_min,
    )

    # Transform just the electron PSD into the explosion frame, since we
    # will need it shortly to calculate inverse Compton emission
    d²N_dpdcos_ef = get_dNdp_2D(
        nc_unit, zone_pop, m_ion, n_ions, n₀_ion,
        psd_lin_cos_bins, γ₀, β₀, psd, num_psd_θ_bins,
        psd_θ_bounds, num_psd_mom_bins, psd_mom_bounds, n_grid,
        γ_sf_grid, i_ion, num_crossings, n_cr_count, therm_grid,
        therm_pₓ_sk, therm_ptot_sk, therm_weight,
        psd_bins_per_dec_mom, psd_mom_min, psd_bins_per_dec_θ, psd_cos_fine,
        Δcos, psd_θ_min)


    # Output the spectra associated with this ion species
    #TODO: spectrum_plot

    # Print out escaping particle population for this species
    print_dNdp_esc(esc_psd_feb_upstream, esc_psd_feb_downstream)

    # Handle photon calculations
    if do_photons
        photon_calcs(
            dNdp_therm_pvals, dNdp_therm, dNdp_cr, aa, n_shell_endpoints, d²N_dpdcos_ef,
            i_ion, mc,
            aa_ion, n₀_ion, u₀, γ₀, u₂,
            flux_px_upstream, flux_energy_upstream, btot_grid,
        )
    end

    return (;
        dNdp_therm, dNdp_therm_pvals, dNdp_cr, zone_pop,
        pressure_psd_par, pressure_psd_perp, energy_density_psd,
    )
end
