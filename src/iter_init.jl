function iter_init()
    # Zero out numerous quantities that will be modified over the course of this iteration.
    # Minimally positive number is used to prevent errors when taking logarithms later
    fill!(pxx_flux, 1e-99)
    fill!(pxz_flux, 1e-99)
    fill!(energy_flux, 1e-99)

    fill!(esc_spectra_feb_UpS, 1e-99)
    fill!(esc_spectra_feb_DwS, 1e-99)

    fill!(pressure_psd_par, 1e-99)
    fill!(pressure_psd_perp, 1e-99)
    fill!(energy_density_psd, 1e-99)

    fill!(weight_coupled, 1e-99)


    # Additionally, set/reset scalar quantities that will change
    ∑P_DwS  = 1e-99 # downstream pressure
    sum_KEdensity_DwS = 1e-99

    energy_esc_UpS    = 1e-99
    px_esc_UpS        = 1e-99


    # To facilitate energy transfer from ions to electrons, calculate here the target energy
    # density fraction for electrons at each grid zone, and zero out the pool of
    # plasma-frame energy that will be taken from ions and donated to electrons
    # Per Ardaneh+ (2015) [10.1088/0004-637X/811/1/57], ε_electron is proportional to √ε_b.
    # ε_b is itself roughly proportional to density² -- B² is proportional to z² (z being the
    # density compression factor) -- so ε_electron should vary roughly linearly with density.
    z_max = γ₀ * β₀ / (γ₂ * β₂)
    for i in eachindex(ε_target)
        if ux_sk_grid[i] != u₀
            z_curr = γ₀ * u₀ / (γ_sf_grid[i] * ux_sk_grid[i])
            ε_target[i] = energy_transfer_frac * (z_curr - 1) / (z_max - 1)
        end
    end

    fill!(energy_transfer_pool, 0.0)
    fill!(energy_recv_pool, 0.0)
    fill!(energy_density, 0.0)
    fill!(therm_energy_density, 0.0)
end
