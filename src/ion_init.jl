function ion_init(i_iter, i_ion, species)
    zz = charge(species[i_ion])
    m = mass(species[i_ion])
    aa = m/u"mp" |> NoUnits

    mc = aa*mₚ_cgs * c_cgs

    # At the start of each ion, print a glyph to the screen
    @info("Starting species iteration", i_iter, i_ion)

    # Determine the pcut at which to switch from low-E particle counts to high-E
    # particle counts. Recall that energy_pcut_hi has units of keV per aa, so when
    # dividing by particle mass the factor of aa is already present in the denominator.
    # Also set the maximum momentum cutoff based on the values given in keyword "ENMAX"
    E_pcut_hi_rmproton = energy_pcut_hi*u"keV" / E₀_proton # FIXME pick better name
    if E_pcut_hi_rmproton < energy_rel_pt
        p_pcut_hi = √(2 * E_pcut_hi_rmproton)
    else
        p_pcut_hi = aa * mₚ_cgs * c_cgs * √((E_pcut_hi_rmproton + 1)^2 - 1)
    end

    if Emax_keV > 0u"keV"
        γ = 1 + Emax_keV/(aa*E₀_proton)
        pmax_cutoff = aa*mₚ_cgs * c_cgs * √(γ^2 - 1)
    elseif Emax_keV_per_aa > 0u"keV"
        γ = 1 + Emax_keV_per_aa/E₀_proton
        pmax_cutoff = aa*mₚ_cgs * c_cgs * √(γ^2 - 1)
    elseif pmax_cgs > 0u"mp*c"
        pmax_cutoff = pmax_cgs
    else
        # Something has gone very wrong.
        error("Max CR energy not set in data_input, so can't set pmax_cutoff.")
    end

    # Zero out the phase space distributions and set variables related to
    # tracking thermal particles
    n_cr_count = clear_psd!()

    # In addition to initializing the phase space distribution, open the scratch (i.e.
    # temporary) file to which we will write information about thermal particle grid crossings
    nc_unit = open("mc_crossings.dat", "a")

    # To maintain identical results between OpenMP and serial versions,
    # set RNG seed based on current iteration/ion/pcut/particle number
    iseed_mod = (i_iter - 1)*n_ions + (i_ion - 1)
    Random.seed!(iseed_mod)


    # Initialize the particle populations that will be propagated through
    # the shock structure
    global n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in, x_PT_cm_in, globals... = init_pop(
        do_fast_push, inp_distr, i_ion, m,
        temperature.(species), energy_inj, inj_weight, n_pts_inj,
        density.(species), x_grid_start, rg₀, η_mfp, x_fast_stop_rg,
        β₀, γ₀, u₀, n_ions, mass.(species),
        n_grid, x_grid_rg, uₓ_sk_grid, γ_sf_grid,
        ptot_inj, weight_inj, n_pts_MB,
    )
    global pₓₓ_flux = globals[1]
    global pxz_flux = globals[2]
    global energy_flux = globals[3]
    @info("Finished init_pop on",
          i_iter, i_ion,
          n_pts_use, i_grid_in, weight_in, ptot_pf_in,
          pb_pf_in, x_PT_cm_in, pₓₓ_flux, pxz_flux, energy_flux)

    # Assign the various particle properties to the population
    assign_particle_properties_to_population!(n_pts_use, xn_per_fine, x_grid_stop)

    # Weight of remaining particles, printed after each pcut; note that this will not be
    # correct for all particles if they were originally created so each thermal bin
    # would have equal weight
    global weight_running = weight_in[1]

    # When using OpenMP, the array energy_transfer_pool can't be conditionally assigned
    # shared or reduction status, so it can't be used for both the ion and electron loops.
    # To get around this, use one array to hold the donated energy, and another to hold
    # the received energy.
    energy_recv_pool .= energy_transfer_pool


    # The array of pcuts read in by data_input has units momentum/mc.
    # Convert to momentum for this species
    pcuts_use[1:n_pcuts] .= pcuts_in[1:n_pcuts] * aa*mₚ_cgs*c_cgs

    return (
            aa, zz, m, mc,
            nc_unit,
           )
end

function clear_psd!()
    fill!(num_crossings,   0)
    fill!(therm_grid,      0)
    fill!(psd,             1e-99)
    fill!(esc_psd_feb_UpS, 1e-99)
    fill!(esc_psd_feb_DwS, 1e-99)
    fill!(therm_pₓ_sk,     0.0)
    fill!(therm_ptot_sk,   0.0)
    fill!(therm_weight,    0.0)
    n_cr_count = 0
    return n_cr_count
end
function assign_particle_properties_to_population!(n_pts_use, xn_per_fine, x_grid_stop)
    weight_new[1:n_pts_use]  .= weight_in[1:n_pts_use]
    ptot_pf_new[1:n_pts_use] .= ptot_pf_in[1:n_pts_use]
    pb_pf_new[1:n_pts_use]   .= pb_pf_in[1:n_pts_use]
    x_PT_cm_new[1:n_pts_use] .= x_PT_cm_in[1:n_pts_use]
    grid_new[1:n_pts_use]    .= i_grid_in[1:n_pts_use]

    DwS_new[1:n_pts_use]         .= false
    inj_new[1:n_pts_use]         .= false
    xn_per_new[1:n_pts_use]      .= xn_per_fine
    prp_x_cm_new[1:n_pts_use]    .= x_grid_stop
    acctime_sec_new[1:n_pts_use] .= 0.0
    tcut_new[1:n_pts_use]        .= 1

    φ_rad_new[1:n_pts_use] .= 2π*Random.rand(n_pts_use)
end
