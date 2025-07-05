module smoothers

using Roots

using ..constants: kB, E₀_proton
using ..parameters: β_rel_fl
using ..io: print_plot_vals

export smooth_grid_part, smooth_profile!

const grid_smoothing_maxitrs = 10_000
const grid_smoothing_target_err = 1e-6

"""
    smooth_grid_par(...)

Uses tracked fluxes of particles, and the Rankine-Hugoniot jump conditions, to determine the
smoothed shock profile for the next iteration of the code. Only valid in parallel case,
which makes subroutine, inputs and equations simpler than would be required for oblique case.

### Arguments
FIXME
- `i_iter`: current iteration number
- `n_grid`: number of grid zones
- `x_grid_rg`: locations[rg₀] of grid zone boundaries
- `uz_sk_grid`: z-component of bulk fluid velocity, in shock frame. This is a parallel shock,
  so `uz_sk_grid = 0` identically; just used for output
- `θ_grid`: angle[rad] between mean magnetic field and shock normal.
  This is a parallel shock, so θ_grid = 0 identically; just used for output
- `pxx_flux`: momentum flux of particles across grid zone boundaries
- `energy_flux`: energy flux of particles across grid zone boundaries
- `Γ₂`: downstream adiabatic index
- `flux_px_UpS`: far upstream momentum flux, calculated in upstream_fluxes
- `flux_energy_UpS`: far upstream energy flux, calculated in upstream_fluxes

### Returns

- `utot_grid`: total bulk flow speed of fluid in shock frame (is a pure output
  because `uₓ_sk_grid` is used in the calculations instead)
- `γ_ef_grid`: Lorentz factor of flow relative to far upstream plasma (i.e. in the explosion frame)
- `β_ef_grid`: bulk flow speed associated with `γ_ef_grid`

### Modifies

- `uₓ_sk_grid`: x-component of bulk fluid velocity, in shock frame
- `γ_sf_grid`: Lorentz factor associated with `utot_grid`, but since this is a parallel shock
  it's the Lorentz factor associated with `uₓ_sk_grid`
- `btot_grid`: magnetic field strength[G] in each grid zone, including any compression of
  turbulence or additional amplification
"""
function smooth_grid_par(
        i_iter, i_shock, n_grid, x_grid_rg, x_grid_cm,
        Γ_grid, uz_sk_grid, θ_grid,
        pressure_psd_par, pressure_psd_perp,
        flux_px_UpS, flux_energy_UpS, Γ₂, q_esc_cal_pₓ, q_esc_cal_energy,
        pxx_flux, energy_flux, uₓ_sk_grid, γ_sf_grid, btot_grid, utot_grid,
        γ_ef_grid, β_ef_grid, εB_grid,
        n_ions, aa_ion, zz_ion, T₀_ion, n₀_ion,
        rg₀, do_prof_fac_damp, prof_weight_fac, γ₀, u₀, β₀,
        γ₂, β₂, u₂, do_smoothing, smooth_mom_energy_fac,
        # ω = smooth_pressure_flux_psd_fac
        ω, bturb_comp_frac, bfield_amp, bmag₀,
        x_art_start_rg, use_custom_εB)

    pxx_norm = zeros(n_grid)
    energy_norm = zeros(n_grid)
    pxz_tot = zeros(n_grid)
    pxz_norm = zeros(n_grid)
    pressure_tot_MC = zeros(n_grid)
    uₓ_new_pₓ = zeros(n_grid)
    uₓ_new_energy = zeros(n_grid)
    uₓ_new = zeros(n_grid)


    # Set constants
    # First, density in units of proton rest mass, pressure, and number density
    n₀ = dot(n₀_ion, aa_ion)
    P₀   = dot(n₀_ion, T₀_ion) * kB

    #DEBUGLINE (for now)
    # Calculate the far upstream magnetization -- the ratio of the energy fluxes in EM fields
    # and particles. This will be used to scale the downstream decay of magnetic field,
    # linking rg₀ to the ion skin depth used in PIC sims.
    ##TODO: this uses γ₀ for the kinetic energy, rather than γ₀ - 1. Okay for ultra-relativistic
    # shocks, but badly mistaken in trans-relativistic limit. Does this affect the results?
    #σ₀ = bmag₀^2 / (4π * γ₀ * n₀ * E₀_proton)

    # Determine weighting factor for profile averaging
    if do_prof_fac_damp && i_iter != 1
        prof_weight_fac *= i_iter < 6 ? 1.15 : 1.50
        prof_weight_fac = max(10.0, prof_weight_fac)
    end

    x_grid_log = zeros(n_grid)
    x_grid_log_cm = zeros(n_grid)
    pxx_tot = zeros(n_grid)
    energy_tot = zeros(n_grid)
    # Compute a bunch of stuff about the current shock profile and print it
    # to file; loop 4111 in old code
    #-------------------------------------------------------------------------
    for i in 1:n_grid

        # Grid coordinates in log space
        if x_grid_rg[i] < -1
            x_grid_log[i] = -log10(-x_grid_rg[i])
        elseif x_grid_rg[i] > 1
            x_grid_log[i] =  log10( x_grid_rg[i])
        else
            x_grid_log[i] = 0.0
        end

        if x_grid_rg[i] <  0
            x_grid_log_cm[i] = -log10(-x_grid_rg[i] * rg₀)
        elseif x_grid_rg[i] > 0
            x_grid_log_cm[i] =  log10( x_grid_rg[i] * rg₀)
        else
            x_grid_log_cm[i] = 0.0
        end

        # Pull from ***_grid arrays into easier-to-use variables
        uₓ     = uₓ_sk_grid[i]
        utot   = utot_grid[i]
        β_uₓ   = uₓ            / c
        β_uz   = uz_sk_grid[i] / c
        γᵤ_sf = γ_sf_grid[i]
        γᵤ_ef = γ_ef_grid[i]
        βᵤ_ef = β_ef_grid[i]
        θ_deg  = rad2deg(θ_grid[i])
        bmag   = btot_grid[i]

        # "Pre" and "Post" refer respectively to the adiabatic index calculated before the
        # particles have propagated through the profile for this iteration, and the
        # adiabatic index calculated using the thermal crossing and PSD info of the most
        # recent iteration
        Γ_pre  = Γ_grid[i,1]
        Γ_post = Γ_grid[i,2]


        # Basic calculations using those variables
        uₓ_norm = uₓ / uₓ_sk_grid[1]
        uz_norm = 1e-99  # parallel shock; set to 0

        γ² = γᵤ_sf^2
        γβ = γᵤ_sf * β_uₓ

        density_ratio = γ₀ * β₀ / γβ

        # Magnetic field components, and associated fluxes according to
        # Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D]
        #TODO: this assumes a mean field, not turbulence. How do the
        # equations change when there's turbulence?
        B_x = bmag * cos(θ_grid[i])
        B_z = bmag * sin(θ_grid[i])

        pxx_EM = γβ^2 / 8π * bmag^2 + γ² / 8π * (B_z^2 - B_x^2) - (γ² - γᵤ_sf) / 2π * (β_uz/β_uₓ) * B_x * B_z

        energy_EM  = γᵤ_sf^2 / 4π * β_uₓ * B_z^2 - (2γ_sq - γᵤ_sf) / 4π * β_uz * B_x * B_z

        # Total momentum/energy fluxes, including electrons (if needed) and EM.
        # Also normalized against far upstream values and in log space for plotting.
        pxx_tot[i] = pxx_flux[i] + pxx_EM
        energy_tot[i] = energy_flux[i] + energy_EM + Γ_post/(Γ_post-1) * uₓ

        pxx_norm[i] = pxx_tot[i] / flux_px_UpS
        energy_norm[i]  = energy_tot[i]  / flux_energy_UpS

        if pxx_norm[i] > 1e-99
            pxx_norm_log = log10(pxx_norm[i])
        else
            pxx_norm_log = -99.0
        end

        if energy_norm[i] > 1e-99
            energy_norm_log = log10(energy_norm[i])
        else
            energy_norm_log = -99.0
        end

        # In a parallel shock, the z-momentum flux is irrelevant. Set it to 0
        pxz_tot[i]   = 1e-99
        pxz_norm[i]  = 1e-99
        pxz_norm_log = -99.0


        # Calculate pressure using the relativistic equations of Double+ (2004)
        # [2004ApJ...600..485D], Eqs (27) and (28) specifically. Combine the resultant
        # pressure using smooth_mom_energy_fac from input file. Note that flux_energy_UpS
        # has the rest mass-energy flux subtracted off, so add it back here
        # Note, too, that q_esc_cal_** is already in units of far upstream flux
        # TODO: per original code, "there is an unresolved question as to whether or not to
        # use the escaping fluxes in these expressions". Using the escaping fluxes sounds
        # reasonable, esp. in the nonrelativistic case. Make sure it's actually reasonable
        pₓ_numer = flux_px_UpS * (1.0 - q_esc_cal_pₓ) - γβ^2 * density_ratio * n₀*E₀_proton
        pₓ_denom = 1 + γβ^2 * Γ_pre / (Γ_pre - 1)
        pressure_pₓ  = pₓ_numer / pₓ_denom

        energy_term_1 = flux_energy_UpS * (1 - q_esc_cal_energy)
        energy_term_2 = γ₀*β₀*c * n₀*E₀_proton
        energy_term_3 = γ² * uₓ * density_ratio * n₀*E₀_proton
        pressure_energy = (energy_term_1 + energy_term_2 - energy_term_3) / (γ² * uₓ * Γ_pre/(Γ_pre-1))

        # These pressures can become negative if a sharp shock with high compression ratio
        # results in a great deal of escaping flux. Place a floor on them for plotting purposes
        pressure_pₓ = max(pressure_pₓ, 1e-99)
        pressure_energy = max(pressure_energy, 1e-99)

        # Use the tabulated pressures from the thermal crossings and the PSD to determine
        # two quantities: the total pressure and the degree of anisotropy. Note that
        # pressure_aniso will return exactly 1.0 if the pressure is isotropic, as
        # pressure_par should be half of pressure_perp
        pressure_tot_MC[i] =  pressure_psd_par[i] + pressure_psd_perp[i]
        pressure_aniso     = 2pressure_psd_par[i] / pressure_psd_perp[i]


        # Calculate the expected downstream pressure in the absence of DSA, i.e. in the test
        # particle limit. No escaping flux to worry about here, but still need to add in the
        # rest mass-energy flux
        if i == 1
            pₓ_numer = flux_px_UpS - γ₂*β₂ * γ₀*Β₀ * n₀*E₀_proton
            pₓ_denom = 1 + (γ₂*β₂)^2 * Γ₂/(Γ₂ - 1)
            pressure_pₓ_tp = pₓ_numer / pₓ_denom

            energy_numer = flux_energy_UpS + γ₀*u₀ * n₀*E₀_proton - γ₂ * c * γ₀*β₀ * n₀*E₀_proton
            energy_denom = γ₂^2 * u₂ * Γ₂/(Γ₂ - 1)
            pressure_energy_tp  = energy_numer / energy_denom
        end

        # Write it all to file
        inquire(file="./mc_grid.dat", opened=lopen)
        if !lopen
            open(newunit=mc_grid_fileunit, status="unknown", file="./mc_grid.dat")
        end

        # WARNING: these column numbers are reused in subroutine read_old_prof.
        # If they are ever modified, change that subroutine accordingly!
        write(mc_grid_fileunit,
              (
               i_iter, i,
               x_grid_rg[i],                    # 1
               x_grid_log[i],                   # 2
               x_grid_cm[i],                    # 3
               x_grid_log_cm[i],                # 4
               pxx_norm[i],                     # 5
               pxx_norm_log,                    # 6
               pxz_norm[i],                     # 7
               pxz_norm_log,                    # 8
               energy_norm[i],                  # 9
               energy_norm_log,                 # 10
               uₓ_norm,                         # 11
               log10(uₓ_norm),                  # 12
               uz_norm,                         # 13
               log10(uz_norm),                  # 14
               bmag,                            # 15
               log10(bmag),                     # 16
               θ_deg,                           # 17
               γᵤ_sf,                           # 18
               1/density_ratio,                 # 19
               density_ratio,                   # 20
               log10(pressure_pₓ),              # 21
               log10(pressure_energy),          # 22
               log10(pressure_psd_par[i]),      # 23
               log10(pressure_psd_perp[i]),     # 24
               log10(pressure_tot_MC[i]),       # 25
               pressure_aniso,                  # 26
               log10(pressure_pₓ_tp),           # 27
               log10(pressure_energy_tp),       # 28
               log10(P₀),                       # 29
               log10(1-q_esc_cal_pₓ),           # 30  Remaining fluxes for plot:
               log10(1-q_esc_cal_energy),       # 31  momentum and energy
               εB_grid[i],                      # 32
               log10(εB_grid[i])                # 33
              )
             )

    end

    print_plot_vals(mc_grid_fileunit)
    #-------------------------------------------------------------------------
    # Grid output completed


    # Return if keeping constant profile
    do_smoothing || return


    #-------------------------------------------------------------------------
    if β₀ < β_rel_fl    # Non-relativistic calculation of new velocity profile
        nonrelativistic_velocity_profile()
    else                # Relativistic calculation of new velocity profile
        # CHECKTHIS: what happens if γ*u is used as a smoothing variable instead
        # of just u? At high speeds γ*u is much more variable than just u
        relativistic_velocity_profile()
    end  # check on shock speed
    #-------------------------------------------------------------------------
    # Relativistic/nonrelativistic shock smoothing complete


    # Artificial smoothing if directed
    if x_art_start_rg < 0
        i_trans = findfirst(>(x_art_start_rg), x_grid_rg) - 1
        uₓ_scale_fac = -(uₓ_new[i_trans] - uₓ_new[n_grid]) / atan(x_grid_rg[i_trans])
        for i in i_trans:i_shock
            uₓ_new[i] = -atan(x_grid_rg[i]) * uₓ_scale_fac + uₓ_new[n_grid]
        end
    end


    # Average with previous profile
    # CHECKTHIS: what happens if γ * β is used as an averaging variable
    # instead of just u for a relativistic shock?
    for i in 1:n_grid
        uₓ_new[i] = (uₓ_new[i] + prof_weight_fac*uₓ_sk_grid[i]) / (1 + prof_weight_fac)
    end


    # Compute output arrays based on new profile
    for i in 1:n_grid
        uₓ_sk_grid[i] = uₓ_new[i]
        γ_sf_grid[i] = 1 / √(1 - (uₓ_sk_grid[i]/c)^2)
        utot_grid[i] = uₓ_new[i]
        β_ef_grid[i] = (u₀ - uₓ_sk_grid[i]) / (c - u₀*uₓ_sk_grid[i]/c)
        γ_ef_grid[i] = 1 / √(1 - β_ef_grid[i]^2)

        # Include necessary corrections for turbulence compression
        z_comp       = (γ₀ * u₀) / (γ_sf_grid[i] * uₓ_sk_grid[i])
        comp_fac     = 1 + (√(1/3 + 2/3 * z_comp^2) - 1) * bturb_comp_frac
        # Also include any additional amplification specified
        amp_fac      = 1 + (comp_fac - 1) * bfield_amp
        btot_grid[i] = bmag₀ * amp_fac

        # If a custom ε_B is in place, use that to calculate the magnetic field, not the preceding code.
        # Note that the R-H relations can be rearranged to read
        #     energy_density(x)  =  F_en₀/u(x) - F_px₀
        # assuming flux conservation everywhere.
        if use_custom_εB
            energy_density = (flux_energy_UpS + γ₀*u₀*n₀*E₀_proton) / uₓ_sk_grid[i] - flux_px_UpS
            btot_grid[i] = √(8π * εB_grid[i] * energy_density)
        end
    end

    return
end

function relativistic_velocity_profile()
    avg_DwS_uₓ_pₓ = 0.0
    avg_DwS_uₓ_energy = 0.0
    Qpₓ = q_esc_cal_pₓ * pxx_flux[1]
    Qen = q_esc_cal_energy * energy_flux[1]

    for i in 1:n_grid
        β_uₓ   = uₓ_sk_grid[i] / c
        γᵤ_sf = γ_sf_grid[i]
        γ²     = γᵤ_sf^2
        γβ     = γᵤ_sf * β_uₓ

        density_loc = γ₀ * β₀ / (γᵤ_sf*β_uₓ) * n₀

        Γ_post = Γ_grid[i,2]

        # Magnetic field components, and associated fluxes according to # Eqs. (27) & (28)
        # of Double+ (2004) [2004ApJ...600..485D], and # assume that uz = 0 (this *is* the
        # parallel smoothing subroutine) TODO: this assumes a mean field, not turbulence.
        # How do the equations change when there's turbulence?
        bmag = btot_grid[i]
        B_x = bmag * cos(θ_grid[i])
        B_z = bmag * sin(θ_grid[i])

        pxx_EM = γβ^2 / 8π * bmag^2 + γ² / 8π * (B_z^2 - B_x^2)
        energy_EM = γᵤ_sf^2 / 4π * β_uₓ * B_z^2

        # Calculate the pressure using the momentum equation only, since the energy equation
        # can give negative fluxes if fast push is used. Do not include EM flux here, since
        # pxx_flux tracked only particle contributions to F_pₓ. Also do not include escaping
        # flux, since we only care about the particles that remain
        pressure_pₓ  = (pxx_flux[i] - γβ^2 * density_loc*E₀_proton) / (1 + γβ^2 * Γ_post/(Γ_post - 1))

        # Combine flux-based pressure and PSD-based pressure as directed by user input
        pressure_loc = (1-ω)*pressure_pₓ + ω*pressure_tot_MC[i]

        # Find new velocity using just the momentum flux equation.
        # Energy flux will be used further down.
        # TODO: neither of these includes EM component of flux, which will be important
        # once turbulence is considered Newton's method, momentum.
        # Use γ*β since that retains the sign information of just uₓ, but also scales to
        # relativistic speeds
        #------------------------------------------------------------------------
        function p(γβ) # momentum
            pₓ_term = γ₀*β₀ * n₀/density_loc * γβ * (density_loc*E₀_proton + pressure_loc * Γ_post/(Γ_post - 1))
            return flux_px_UpS - Qpₓ - pxx_EM - pₓ_term - pressure_loc
        end
        γβ_found = Roots.find_zero(p, γ₀*β₀*1e-4, Roots.Newton())
        uₓ_new_pₓ[i] = γβ_found / √(1 + γβ_found^2) * c


        # Newton's method, energy
        #-------------------------------------------------------------------
        function E(γβ)
            γ = √(1 + γβ^2)
            energy_term = γβ * γ * c * (density_loc*E₀_proton + Γ_post / (Γ_post - 1) * pressure_loc)
            return flux_energy_UpS - Qen - energy_EM - energy_term
        end
        γβ_found = Roots.find_zero(E, γ₀*β₀*1e-4, Roots.Newton())
        uₓ_new_energy[i] = γβ_found / √(1 + γβ_found^2) * c

        # Find the average downstream velocity, which will be used for
        # scaling the profiles. Only average the last 10 grid positions
        if i > (n_grid-10)
            avg_DwS_uₓ_pₓ += uₓ_new_pₓ[i]
            avg_DwS_uₓ_energy += uₓ_new_energy[i]
        end

    end  # loop over grid positions


    # Smooth the velocity profile, rescale it, and average the momentum and
    # energy curves as directed by user input
    smooth_profile!(uₓ_new_pₓ, n_grid)
    smooth_profile!(uₓ_new_energy, n_grid)

    avg_DwS_uₓ_pₓ /= 10
    avg_DwS_uₓ_energy /= 10

    uₓ_scale_fac = (u₀ - u₂)/(uₓ_new_pₓ[1] - avg_DwS_uₓ_pₓ)
    for i in 1:n_grid
        uₓ_new_pₓ[i] = uₓ_scale_fac * (uₓ_new_pₓ[i] - avg_DwS_uₓ_pₓ) + u₂
        if x_grid_rg[i] ≥ 0
            uₓ_new_pₓ[i] = u₂
        end
    end

    uₓ_scale_fac = (u₀ - u₂)/(uₓ_new_energy[1] - avg_DwS_uₓ_energy)
    for i in 1:n_grid
        uₓ_new_energy[i] = uₓ_scale_fac * (uₓ_new_energy[i] - avg_DwS_uₓ_energy) + u₂
        if x_grid_rg[i] ≥ 0
            uₓ_new_energy[i] = u₂
        end
    end

    uₓ_new .= @. (1-smooth_mom_energy_fac)*uₓ_new_pₓ + smooth_mom_energy_fac*uₓ_new_energy
end

function nonrelativistic_velocity_profile()
    avg_DwS_uₓ_pₓ = 0.0
    avg_DwS_uₓ_energy = 0.0
    Qpₓ = 0.0  # By default for nonrelativistic shocks
    Qen = q_esc_cal_energy * energy_flux[1]

    for i in 1:n_grid
        uₓ    = uₓ_sk_grid[i]
        β_uₓ  = uₓ / c
        γᵤ_sf = γ_sf_grid[i]
        γ²    = γᵤ_sf^2
        γβ    = γᵤ_sf * β_uₓ

        Γ_post = Γ_grid[i,2]

        # Magnetic field components, and associated fluxes according to Eqs. (27) & (28) of
        # Double+ (2004) [2004ApJ...600..485D], and assume that uz = 0 (this *is* the
        # parallel smoothing subroutine) TODO: this assumes a mean field, not turbulence.
        # How do the equations change when there's turbulence?
        bmag   = btot_grid[i]
        B_x    = bmag * cos(θ_grid[i])
        B_z    = bmag * sin(θ_grid[i])

        pxx_EM = γβ^2 * bmag^2 / 8π + γ² * (B_z^2 - B_x^2) / 8π
        energy_EM = γᵤ_sf^2 / 4π * β_uₓ * B_z^2

        # Calculate the pressure using the momentum equation only, since the
        # energy equation can give negative fluxes if fast push is used.
        # Determining pressure relies on near cancellation of two terms, so use
        # a form for the non-relativistic equations that is expanded to β^2 to
        # allow for better joining between relativistic and non-relativistic
        # calculations. Do not include EM flux here, since pxx_flux tracked only
        # particle contributions to F_pₓ. Also do not include escaping flux,
        # since we only care about the particles that remain
        pressure_pₓ = (pxx_flux[i] - n₀*mp * u₀ * uₓ * (1+β_uₓ^2)) / (1 + β_uₓ^2 * Γ_post/(Γ_post - 1))

        # Combine flux-based pressure and PSD-based pressure as directed by user input
        pressure_loc = (1-ω)*pressure_pₓ + ω*pressure_tot_MC[i]

        # Find new velocity using newly-found pressure; here use both momentum
        # *and* energy equations. Need to include EM flux at this stage because
        # flux_**_UpS included it. Because we are keeping terms out to β² now,
        # the momentum and energy equations go from linear/quadratic to
        # cubic/quartic. Instead of solving the equations analytically, use
        # Newton's method.

        function p(u) # momentum
            β = u / c
            p_term_1 = n₀*mp * u₀ * uₓ_guess * (1 + β^2)
            p_term_2 = (1 + β^2 * Γ_post/(Γ_post - 1)) * pressure_loc
            return flux_px_UpS - Qpₓ - pxx_EM - p_term_1 - p_term_2
        end
        uₓ_new_pₓ[i] = Roots.find_zero(p, u₀ * 1e-4, Roots.Newton())

        function E(u) # energy
            β = u / c
            energy_term_1 = 1//2 * n₀*mp * u₀ * u^2 * (1 + 1.25*β^2)
            energy_term_2 = Γ_post / (Γ_post - 1) * pressure_loc * u * (1 + β^2)
            return flux_energy_UpS - Qen - energy_EM - energy_term_1 - energy_term_2
        end
        uₓ_new_energy[i] = Roots.find_zero(E, u₀ * 1e-4, Roots.Newton())


        # Find the average downstream velocity, which will be used for scaling
        # the profiles. Only average the last 10 grid positions
        if i > n_grid-10
            avg_DwS_uₓ_pₓ += uₓ_new_pₓ[i]
            avg_DwS_uₓ_energy += uₓ_new_energy[i]
        end

    end  # loop over grid positions


    # Scale the velocity profile, smooth it, and average the momentum and energy curves
    # as directed by user input
    avg_DwS_uₓ_pₓ /= 10
    avg_DwS_uₓ_energy /= 10

    uₓ_scale_fac = (u₀ - u₂)/(uₓ_new_pₓ[1] - avg_DwS_uₓ_pₓ)
    for i in 1:n_grid
        uₓ_new_pₓ[i] = uₓ_scale_fac * (uₓ_new_pₓ[i] - avg_DwS_uₓ_pₓ) + u₂
        if x_grid_rg[i] ≥ 0.0
            uₓ_new_pₓ[i] = u₂
        end
    end

    uₓ_scale_fac = (u₀ - u₂)/(uₓ_new_energy[1] - avg_DwS_uₓ_energy)
    for i in 1:n_grid
        uₓ_new_energy[i] = uₓ_scale_fac * (uₓ_new_energy[i] - avg_DwS_uₓ_energy) + u₂
        if x_grid_rg[i] ≥ 0.0
            uₓ_new_energy[i] = u₂
        end
    end

    smooth_profile!(uₓ_new_pₓ, n_grid)
    smooth_profile!(uₓ_new_energy, n_grid)

    uₓ_new = @. (1-smooth_mom_energy_fac)*uₓ_new_pₓ + smooth_mom_energy_fac*uₓ_new_energy
end

"""
    smooth_profile!(y_prof, n_grid)

Takes an input velocity profile and performs two tasks: enforces monotonicity by
smoothing out dips/bumps, and averages nearby points to smooth sharp edges

### Arguments

- `y_prof`: array holding velocity profile (modified in-place)
- `n_grid`: number of grid zones in position and velocity arrays
"""
function smooth_profile!(y_prof, n_grid)
    # Run from downstream to upstream and eliminate dips/bumps
    for i in n_grid:-1:2
        if y_prof[i-1] < y_prof[i]
            y_prof[i-1] = y_prof[i]
        end
    end

    # Now smooth sharp edges by averaging adjacent grid locations;
    # handle the edges separately for different weighting
    y_prof_dup = similar(y_prof) # only 2:n_grid-1 actually assigned
    y_prof_dup[2] = (2y_prof[1] + y_prof[2] + y_prof[3]) / 4
    for i in 3:n_grid-2
        y_prof_dup[i] = (y_prof[i-1] + y_prof[i] + y_prof[i+1]) / 3
    end
    y_prof_dup[n_grid-1] = (y_prof[n_grid-2] + y_prof[n_grid-1] + 2y_prof[n_grid]) / 4

    for i in 2:n_grid-1
        y_prof[i] = y_prof_dup[i]
    end
end
end # module
