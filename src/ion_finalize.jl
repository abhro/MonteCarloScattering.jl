function ion_finalize()
    println()
    println(outfile)

    # Obtain the dN/dp arrays for the current species: one 2-D array for each grid zone.
    get_normalized_dNdp(dNdp_therm, dNdp_therm_pvals, dNdp_cr, nc_unit, zone_pop)

    # Get pressure (both components) and kinetic energy density everywhere in the shock profile
    thermo_calcs(num_crossings, n_cr_count, therm_grid, therm_pₓ_sk, therm_ptot_sk,
                 therm_weight, nc_unit, psd, zone_pop,
                 pressure_psd_par, pressure_psd_perp, energy_density_psd)

    # Transform just the electron PSD into the explosion frame, since we
    # will need it shortly to calculate inverse Compton emission
    get_dNdp_2D(nc_unit, zone_pop, d²N_dpdcos_ef)


    # Output the spectra associated with this ion species
    #TODO: spectrum_plot

    # Print out escaping particle population for this species
    print_dNdp_esc(esc_psd_feb_UpS, esc_psd_feb_DwS)

    # Handle photon calculations
    if do_photons
      photon_calcs(dNdp_therm_pvals, dNdp_therm, dNdp_cr, aa, n_shell_endpoints, d²N_dpdcos_ef)
    end
end
