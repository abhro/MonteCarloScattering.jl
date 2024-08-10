module smoothers

using ..constants: kB_cgs, me_cgs, mp_cgs, E₀_proton, c_cgs
using ..parameters: na_grid, na_ions, β_rel_fl
using ..io: print_plot_vals

export smooth_grid_part, smooth_profile!

const grid_smoothing_maxitrs = 10_000
const grid_smoothing_target_err = 1e-6

"""
Uses tracked fluxes of particles, and the Rankine-Hugoniot jump conditions, to determine the
smoothed shock profile for the next iteration of the code. Only valid in parallel case,
which makes subroutine, inputs and equations simpler than would be required for oblique case.

### Arguments
FIXME
- i_iter: current iteration number
- n_grid: number of grid zones
- x_grid_rg: locations[rg0] of grid zone boundaries
- uz_sk_grid: z-component of bulk fluid velocity, in shock frame. This is a parallel shock,
  so uz_sk_grid = 0 identically; just used for output
- θ_grid: angle[rad] between mean magnetic field and shock normal.
  This is a parallel shock, so θ_grid = 0 identically; just used for output
- pxx_flux: momentum flux of particles across grid zone boundaries
- energy_flux: energy flux of particles across grid zone boundaries
- γ_adiab_2: DwS adiabatic index
- flux_px_UpS: far UpS momentum flux, calculated in upstream_fluxes
- flux_energy_UpS: far UpS energy flux, calculated in upstream_fluxes

### Returns

- utot_grid: total bulk flow speed of fluid in shock frame (is a pure output b/c
  ux_sk_grid is used in the calculations instead)
- γ_ef_grid: Lorentz factor of flow relative to far UpS plasma (i.e. in the explosion frame)
- β_ef_grid: bulk flow speed associated with γ_ef_grid

### Modifies

- ux_sk_grid: x-component of bulk fluid velocity, in shock frame
- γ_sf_grid: Lorentz factor associated with utot_grid, but since this is a parallel shock
  it's the Lorentz factor associated with ux_sk_grid
- btot_grid: magnetic field strength[G] in each grid zone, including any compression of
  turbulence or additional amplification
"""
function smooth_grid_par(
        i_iter, i_shock, n_grid, x_grid_rg, x_grid_cm,
        γ_adiab_grid, uz_sk_grid, θ_grid,
        pressure_psd_par, pressure_psd_perp,
        flux_px_UpS, flux_energy_UpS, γ_adiab_2, q_esc_cal_px, q_esc_cal_energy,
        pxx_flux, energy_flux, ux_sk_grid, γ_sf_grid, btot_grid, utot_grid,
        γ_ef_grid, β_ef_grid, εB_grid,
        n_ions, aa_ion, zz_ion, tZ_ion, denZ_ion, sc_electron,
        tZ_electron, rg0, do_prof_fac_damp, prof_wt_fac, γ_Z, u_Z, β_Z,
        γ_2, β_2, u_2, do_smoothing, smooth_mom_energy_fac,
        smooth_pressure_flux_psd_fac, bturb_comp_frac, bfield_amp, bmag_Z,
        x_art_start_rg, use_custom_εB)

    pxx_norm = zeros(na_grid)
    energy_norm = zeros(na_grid)
    pxz_tot = zeros(na_grid)
    pxz_norm = zeros(na_grid)
    pressure_tot_MC = zeros(na_grid)
    ux_new_px = zeros(na_grid)
    ux_new_energy = zeros(n_grid) # zeros(na_grid)
    ux_new = zeros(n_grid) # zeros(na_grid)


    # Set constants
    # First, density in units of proton rest mass, pressure, and number density
    density_Z    = dot(denZ_ion, aa_ion)
    pressure_Z   = dot(denZ_ion, tZ_ion) * kB_cgs
    mask = (aa_ion[j] ≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])
    if ! sc_electron
        density_Z  += (me_cgs/mp_cgs)*density_electron
        pressure_Z += density_electron*kB_cgs*tZ_electron
    end


    #DEBUGLINE (for now)
    # Calculate the far upstream magnetization -- the ratio of the energy fluxes
    # in EM fields and particles. This will be used to scale the DwS decay
    # of magnetic field, linking rg0 to the ion skin depth used in PIC sims.
    ##TODO: this uses γ_Z for the KE, rather than γ_Z - 1. Okay for ultra-rel
    # shocks, but badly mistaken in transrel limit. Does this affect the results?
    #σ_Z = bmag_Z^2 / (4π * γ_Z * density_Z * E₀_proton)


    # Determine weighting factor for profile averaging
    if do_prof_fac_damp && i_iter != 1
        prof_wt_fac *=  i_iter < 6 ? 1.15 : 1.50
        prof_wt_fac = max(10.0, prof_wt_fac)
    end


    x_grid_log = zeros(na_grid)
    x_grid_log_cm = zeros(na_grid)
    pxx_tot = zeros(na_grid)
    energy_tot = zeros(na_grid)
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
            x_grid_log_cm[i] = -log10(-x_grid_rg[i] * rg0)
        elseif x_grid_rg[i] > 0
            x_grid_log_cm[i] =  log10( x_grid_rg[i] * rg0)
        else
            x_grid_log_cm[i] = 0.0
        end


        # Pull from ***_grid arrays into easier-to-use variables
        ux    = ux_sk_grid[i]
        uz    = uz_sk_grid[i]
        utot  = utot_grid[i]
        β_ux  = ux / c_cgs
        β_uz  = uz / c_cgs
        γ_usf = γ_sf_grid[i]
        γ_uef = γ_ef_grid[i]
        β_uef = β_ef_grid[i]
        θ_deg = rad2deg(θ_grid[i])
        bmag  = btot_grid[i]


        # "Pre" and "Post" refer respectively to the adiabatic index calculated before the
        # particles have propagated through the profile for this iteration, and the
        # adiabatic index calculated using the thermal crossing and PSD info of the most
        # recent iteration
        γ_adiab_pre  = γ_adiab_grid[i,1]
        γ_adiab_post = γ_adiab_grid[i,2]


        # Basic calculations using those variables
        ux_norm   = ux / ux_sk_grid[1]
        uz_norm   = 1e-99  # parallel shock; set to 0

        γ_sq = γ_usf^2
        γ_β  = γ_usf * β_ux

        density_ratio = γ_Z * β_Z / γ_β

        # Compute pressure due to electrons if they weren't included self-consistently as a
        # separate species. Equation used to do so is an extension of
        # Ellison & Moebius (1987) [1987ApJ...318..474E] Eq (2)
        if sc_electron
            pressure_electron = 0.0
        else

            # This is (γ_2 - 1)/(γ_0 - 1). Note assumption that γ_0 = 5/3.  #assumecold
            #CHECKTHIS: which adiabaic index (pre or post) gives sensible values?
            pressure_term_1 = (γ_adiab_post - 1) * 3//2
            #  Here also: #assumecold
            pressure_term_2 = 1  +  5//3 * γ_uef^2*β_uef / γ_β

            pressure_electron = (density_electron * kB_cgs * tZ_electron *
                                 pressure_term_1 * pressure_term_2)
        end


        # Magnetic field components, and associated fluxes according to
        # Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D]
        #TODO: this assumes a mean field, not turbulence. how do the
        # equations change when there's turbulence?
        B_x = bmag * cos(θ_grid[i])
        B_z = bmag * sin(θ_grid[i])

        pxx_EM = γ_β^2 / 8π * bmag^2 +  γ_sq / 8π * (B_z^2 - B_x^2) - (γ_sq - γ_usf) / 2π * (β_uz/β_ux) * B_x * B_z

        energy_EM  = γ_usf^2 / 4π * β_ux * B_z^2 - (2γ_sq - γ_usf) / 4π * β_uz * B_x * B_z


        # Total momentum/energy fluxes, including electrons (if needed) and EM.
        # Also normalized against far UpS values and in log space for plotting.
        pxx_tot[i] = pxx_flux[i] + pressure_electron + pxx_EM
        energy_tot[i] = energy_flux[i] + energy_EM +  γ_adiab_post/(γ_adiab_post-1) * pressure_electron * ux

        pxx_norm[i] = pxx_tot[i] / flux_px_UpS
        energy_norm[i]  = energy_tot[i]  / flux_energy_UpS

        if pxx_norm[i] > 1e-99
            pxx_norm_log = log10( pxx_norm[i] )
        else
            pxx_norm_log = -99.0
        end

        if energy_norm[i] > 1e-99
            energy_norm_log = log10( energy_norm[i] )
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
        # Note, too, that q_esc_cal_** is already in units of far UpS flux
        # TODO: per original code, "there is an unresolved question as to whether or not to
        # use the escaping fluxes in these expressions". Using the escaping fluxes sounds
        # reasonable, esp. in the nonrel case. Make sure it's actually reasonable
        px_term_1 = flux_px_UpS * (1.0 - q_esc_cal_px) - γ_β^2 * density_ratio * density_Z*E₀_proton
        px_term_2 = 1  +  γ_β^2 * γ_adiab_pre / (γ_adiab_pre - 1)
        pressure_px  = px_term_1 / px_term_2

        energy_term_1 = flux_energy_UpS * (1 - q_esc_cal_energy) +  γ_Z*β_Z*c_cgs * density_Z*E₀_proton -  γ_sq * ux * density_ratio * density_Z*E₀_proton
        energy_term_2 = γ_sq * ux * γ_adiab_pre / (γ_adiab_pre - 1)
        pressure_energy  = energy_term_1 / energy_term_2

        # These pressures can become negative if a sharp shock with high compression ratio
        # results in a great deal of escaping flux. Place a floor on them for plotting purposes
        pressure_px = max(pressure_px, 1e-99)
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
            density_rat_tp = γ_Z * β_Z / (γ_2*β_2)

            px_term_1 = flux_px_UpS -  (γ_2*β_2)^2 * density_rat_tp * density_Z*E₀_proton
            px_term_2 = 1  +  (γ_2*β_2)^2 * γ_adiab_2 / (γ_adiab_2 - 1)
            pressure_px_tp = px_term_1 / px_term_2

            energy_term_1 = flux_energy_UpS  +  γ_Z*u_Z * density_Z*E₀_proton  -  γ_2^2 * u_2 * density_rat_tp * density_Z*E₀_proton
            energy_term_2 = γ_2^2 * u_2 * γ_adiab_2/(γ_adiab_2 - 1)
            pressure_energy_tp  = energy_term_1 / energy_term_2
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
               ux_norm,                         # 11
               log10(ux_norm),                  # 12
               uz_norm,                         # 13
               log10(uz_norm),                  # 14
               bmag,                            # 15
               log10(bmag),                     # 16
               θ_deg,                           # 17
               γ_usf,                           # 18
               1/density_ratio,                 # 19
               density_ratio,                   # 20
               log10(pressure_px),              # 21
               log10(pressure_energy),          # 22
               log10(pressure_psd_par[i]),      # 23
               log10(pressure_psd_perp[i]),     # 24
               log10(pressure_tot_MC[i]),       # 25
               pressure_aniso,                  # 26
               log10(pressure_px_tp),           # 27
               log10(pressure_energy_tp),       # 28
               log10(pressure_Z),               # 29
               log10(1-q_esc_cal_px),           # 30  Remaining fluxes for plot:
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


    # Non-rel calculation of new velocity profile
    #-------------------------------------------------------------------------
    if β_Z < β_rel_fl

        ave_DwS_ux_px = 0.0
        ave_DwS_ux_energy = 0.0
        Qpx = 0.0  # By default for nonrel shocks
        Qen = q_esc_cal_energy * energy_flux[1]

        for i in 1:n_grid
            ux    = ux_sk_grid[i]
            β_ux  = ux / c_cgs
            γ_usf = γ_sf_grid[i]
            γ_sq  = γ_usf^2
            γ_β   = γ_usf * β_ux

            γ_adiab_post = γ_adiab_grid[i,2]

            # Magnetic field components, and associated fluxes according to
            # Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D], and
            # assume that uz = 0 (this *is* the parallel smoothing subroutine)
            # TODO: this assumes a mean field, not turbulence. how do the equations change when there's turbulence?
            bmag   = btot_grid[i]
            B_x    = bmag * cos(θ_grid[i])
            B_z    = bmag * sin(θ_grid[i])

            pxx_EM = γ_β^2 / 8π * bmag^2 +  γ_sq / 8π * (B_z^2 - B_x^2)

            energy_EM = γ_usf^2 / 4π * β_ux * B_z^2


            # Calculate the pressure using the momentum equation only, since the energy
            # equation can give negative fluxes if fast push is used. Determining pressure
            # relies on near cancellation of two terms, so use a form for the non-rel
            # equations that is expanded to β^2 to allow for better joining between rel and
            # non-rel calculations. Do not include EM flux here, since pxx_flux tracked only
            # particle contributions to F_px. Also do not include escaping flux, since we
            # only care about the particles that remain
            pressure_px = (pxx_flux[i]  -  density_Z*mp_cgs * u_Z * ux * (1+β_ux^2)) / ( 1  +  β_ux^2 * γ_adiab_post/(γ_adiab_post - 1) )

            # Combine flux-based pressure and PSD-based pressure as directed by user input
            pressure_loc = (1 - smooth_pressure_flux_psd_fac) * pressure_px + smooth_pressure_flux_psd_fac * pressure_tot_MC[i]

            # Find new velocity using newly-found pressure; here use both momentum *and*
            # energy equations. Need to include EM flux at this stage because flux_**_UpS
            # included it. Because we are keeping terms out to β^2 now, the momentum and
            # energy equations go from linear/quadratic to cubic/quartic. Instead of solving
            # the equations analytically, use Newton's method.
            # Newton's method, momentum
            #---------------------------------------------------------------------
            ux_guess     = u_Z * 1e-4
            ux_guess_p   = u_Z * 1.0001
            Δux_guess    = ux_guess_p - ux_guess

            for j in 1:grid_smoothing_maxitrs

                # Calculate flux differences associated with both ux_guess and ux_guess_p
                β_ux      = ux_guess / c_cgs
                px_term_1 = density_Z*mp_cgs * u_Z * ux_guess * (1 + β_ux^2)
                px_term_2 = (1 +  β_ux^2 * γ_adiab_post / (γ_adiab_post - 1) )  *  pressure_loc

                flux_diff = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2

                β_ux      = ux_guess_p / c_cgs
                px_term_1 = density_Z*mp_cgs * u_Z * ux_guess_p * (1 + β_ux^2)
                px_term_2 = (1 +  β_ux^2 * γ_adiab_post / (γ_adiab_post - 1) )  *  pressure_loc

                flux_diff_p = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2


                # Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
                flux_diff_prime = (flux_diff_p - flux_diff) / Δux_guess


                # Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
                ux_guess_next = ux_guess  -  flux_diff/flux_diff_prime


                # Relative change in this step
                err_curr = (ux_guess_next - ux_guess) / ux_guess

                # If the relative change is small enough, we've found our solution
                #   and can exit the loop; otherwise return for another cycle
                if abs(err_curr) < grid_smoothing_target_err
                    ux_found = ux_guess_next
                    break
                end

                ux_guess   = (ux_guess_next + ux_guess) * 0.50
                ux_guess_p = 1.0010 * ux_guess
                Δux_guess  = ux_guess_p - ux_guess

            end

            # Did we hit the maximum number of iterations without finding the flux-conserving solution?
            if j ≥ grid_smoothing_maxitrs
                @warn "smooth_grid_par: no Newton's method solution for px flux, grid zone $i"
                ux_found = ux_new_px[i-1]
            end

            ux_new_px[i] = ux_found


            # Newton's method, energy
            #------------------------------------------------------------------------
            ux_guess   = u_Z * 1e-4
            ux_guess_p = u_Z * 1.0001
            Δux_guess  = ux_guess_p - ux_guess

            for j in 1:grid_smoothing_maxitrs

                # Calculate flux differences associated with both ux_guess and ux_guess_p
                β_ux      = ux_guess   / c_cgs
                energy_term_1 = 0.5 * density_Z*mp_cgs * u_Z * ux_guess^2 * (1 + 1.25*β_ux^2)
                energy_term_2 = γ_adiab_post / (γ_adiab_post - 1) * pressure_loc * ux_guess   * (1 + β_ux^2)

                flux_diff = flux_energy_UpS - Qen - energy_EM - energy_term_1 - energy_term_2

                β_ux      = ux_guess_p / c_cgs
                energy_term_1 = 0.5 * density_Z*mp_cgs * u_Z * ux_guess_p^2 * (1 + 1.25*β_ux^2)
                energy_term_2 = γ_adiab_post / (γ_adiab_post - 1) * pressure_loc * ux_guess_p * (1 + β_ux^2)

                flux_diff_p = flux_energy_UpS - Qen - energy_EM - energy_term_1 - energy_term_2


                # Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
                flux_diff_prime = (flux_diff_p - flux_diff) / Δux_guess


                # Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
                ux_guess_next = ux_guess  -  flux_diff/flux_diff_prime


                # Relative change in this step
                err_curr = (ux_guess_next - ux_guess) / ux_guess

                # If the relative change is small enough, we've found our solution
                #   and can exit the loop; otherwise return for another cycle
                if abs(err_curr) < grid_smoothing_target_err
                    ux_found = ux_guess_next
                    break
                end

                ux_guess     = (ux_guess_next + ux_guess) * 0.5
                ux_guess_p   = 1.001 * ux_guess
                Δux_guess    = ux_guess_p - ux_guess

            end

            # Did we hit the maximum number of iterations without finding the flux-conserving solution?
            if j ≥ grid_smoothing_maxitrs
                @warn "smooth_grid_par: no Newton's method solution for en flux, grid zone $i"
                ux_found = ux_new_energy[i-1]
            end

            ux_new_energy[i] = ux_found


            # Find the average downstream velocity, which will be used for scaling
            # the profiles. Only average the last 10 grid positions
            if i > n_grid-10
                ave_DwS_ux_px += ux_new_px[i]
                ave_DwS_ux_energy += ux_new_energy[i]
            end

        end  # loop over grid positions


        # Scale the velocity profile, smooth it, and average the momentum and
        # energy curves as directed by user input
        ave_DwS_ux_px /= 10
        ave_DwS_ux_energy /= 10

        ux_scale_fac = (u_Z - u_2)/(ux_new_px[1] - ave_DwS_ux_px)
        for i in 1:n_grid
            ux_new_px[i] = ux_scale_fac * (ux_new_px[i] - ave_DwS_ux_px)  +  u_2
            if x_grid_rg[i] ≥ 0.0
                ux_new_px[i] = u_2
            end
        end

        ux_scale_fac = (u_Z - u_2)/(ux_new_energy[1] - ave_DwS_ux_energy)
        for i in 1:n_grid
            ux_new_energy[i] = ux_scale_fac * (ux_new_energy[i] - ave_DwS_ux_energy)  +  u_2
            if x_grid_rg[i] ≥ 0.0
                ux_new_energy[i] = u_2
            end
        end

        smooth_profile(n_grid, ux_new_px)
        smooth_profile(n_grid, ux_new_energy)

        for i in 1:n_grid
            ux_new[i] = (1 - smooth_mom_energy_fac) * ux_new_px[i] + smooth_mom_energy_fac * ux_new_energy[i]
        end

    end  # check on shock speed
    #-------------------------------------------------------------------------
    # Non-rel shock smoothing complete


    # Rel calculation of new velocity profile
    # CHECKTHIS: what happens if γ*u is used as a smoothing variable instead
    # of just u? At high speeds γ*u is much more variable than just u
    #-------------------------------------------------------------------------
    if β_Z ≥ β_rel_fl

        ave_DwS_ux_px = 0.0
        ave_DwS_ux_energy = 0.0
        Qpx = q_esc_cal_px * pxx_flux[1]
        Qen = q_esc_cal_energy * energy_flux[1]

        for i in 1:n_grid
            ux    = ux_sk_grid[i]
            β_ux  = ux / c_cgs
            γ_usf = γ_sf_grid[i]
            γ_sq  = γ_usf^2
            γ_β   = γ_usf * β_ux

            density_loc = γ_Z * β_Z / (γ_usf*β_ux) * density_Z

            γ_adiab_post = γ_adiab_grid(i,2)

            # Magnetic field components, and associated fluxes according to
            # Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D], and
            # assume that uz = 0 (this *is* the parallel smoothing subroutine)
            # TODO: this assumes a mean field, not turbulence. how do the equations change when there's turbulence?
            bmag   = btot_grid[i]
            B_x    = bmag * cos(θ_grid[i])
            B_z    = bmag * sin(θ_grid[i])

            pxx_EM = γ_β^2 / 8π * bmag^2 +  γ_sq / 8π * (B_z^2 - B_x^2)

            energy_EM  = γ_usf^2 / 4π * β_ux * B_z^2

            # Calculate the pressure using the momentum equation only, since
            # the energy equation can give negative fluxes if fast push is used
            # Do not include EM flux here, since pxx_flux tracked only particle contributions to F_px
            # Also do not include escaping flux, since we only care about the particles that remain
            px_term_1 = pxx_flux[i]  -  γ_β^2 * density_loc*E₀_proton
            px_term_2 = 1  +  γ_β^2 * γ_adiab_post / (γ_adiab_post - 1)
            pressure_px  = px_term_1 / px_term_2

            # Combine flux-based pressure and PSD-based pressure as directed by user input
            pressure_loc = (1.0 - smooth_pressure_flux_psd_fac) * pressure_px + smooth_pressure_flux_psd_fac * pressure_tot_MC[i]

            # Find new velocity using just the momentum flux equation.
            # Energy flux will be used further down.
            # TODO: neither of these includes EM component of flux, which will be important
            # once turbulence is considered Newton's method, momentum.
            # Use γ*β since that retains the sign information of just ux, but also scales to
            # relativistic speeds
            #------------------------------------------------------------------------
            gb_guess   = γ_Z * β_Z * 1e-4
            gb_guess_p = gb_guess * 1.0001
            Δgb_guess  = gb_guess_p - gb_guess

            for j in 1:grid_smoothing_maxitrs

                # Calculate flux differences associated with both gb_guess and gb_guess_p
                px_term_1 = γ_Z * β_Z * density_Z / density_loc  * gb_guess * ( density_loc*E₀_proton  +  pressure_loc * γ_adiab_post / (γ_adiab_post - 1) )
                px_term_2 = pressure_loc
                flux_diff = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2

                px_term_1 = γ_Z * β_Z * density_Z / density_loc  * gb_guess_p * ( density_loc*E₀_proton  +  pressure_loc * γ_adiab_post / (γ_adiab_post - 1) )
                px_term_2 = pressure_loc
                flux_diff_p = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2

                # Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
                flux_diff_prime = (flux_diff_p - flux_diff) / Δgb_guess

                # Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
                gb_guess_next = gb_guess  -  flux_diff/flux_diff_prime

                # Relative change in this step
                err_curr = (gb_guess_next - gb_guess) / gb_guess

                # If the relative change is small enough, we've found our solution
                #   and can exit the loop; otherwise return for another cycle
                if abs(err_curr) < grid_smoothing_target_err
                    gb_found = gb_guess_next
                    break
                end

                gb_guess   = (gb_guess_next + gb_guess) / 2
                gb_guess_p = 1.001 * gb_guess
                Δgb_guess  = gb_guess_p - gb_guess
            end

            # Did we hit the maximum number of iterations without finding the flux-conserving solution?
            if j ≥ grid_smoothing_maxitrs
                @warn "smooth_grid_par: no Newton's method solution for px flux, grid zone $i"
                gb_found = ux_new_px[i] / √( c_cgs^2 - ux_new_px[i-1]^2 )
            end

            ux_new_px[i] = gb_found / √( 1 + gb_found^2 ) * c_cgs


            # Newton's method, energy
            #-------------------------------------------------------------------
            gb_guess   = γ_Z * β_Z * 1e-4
            gb_guess_p = gb_guess * 1.0001
            Δgb_guess  = gb_guess_p - gb_guess

            for j in 1:grid_smoothing_maxitrs

                # Calculate flux differences associated with both gb_guess and gb_guess_p
                γ_ux    = √( 1 + gb_guess^2 )
                energy_term_1 = gb_guess   * γ_ux * c_cgs * ( density_loc*E₀_proton  +  γ_adiab_post / (γ_adiab_post - 1) * pressure_loc )
                flux_diff = flux_energy_UpS - Qen - energy_EM - energy_term_1

                γ_ux    = √( 1 + gb_guess_p^2 )
                energy_term_1 = gb_guess_p * γ_ux * c_cgs * ( density_loc*E₀_proton  +  γ_adiab_post / (γ_adiab_post - 1) * pressure_loc )
                flux_diff_p = flux_energy_UpS - Qen - energy_EM - energy_term_1


                # Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
                flux_diff_prime = (flux_diff_p - flux_diff) / Δgb_guess

                # Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
                gb_guess_next = gb_guess  -  flux_diff/flux_diff_prime


                # Relative change in this step
                err_curr = (gb_guess_next - gb_guess) / gb_guess

                # If the relative change is small enough, we've found our solution
                #  and can exit the loop; otherwise return for another cycle
                if abs(err_curr) < grid_smoothing_target_err
                    gb_found = gb_guess_next
                    break
                end

                gb_guess   = (gb_guess_next + gb_guess)/2
                gb_guess_p = 1.001 * gb_guess
                Δgb_guess  = gb_guess_p - gb_guess

            end

            # Did we hit the maximum number of iterations without finding the flux-conserving solution?
            if j ≥ grid_smoothing_maxitrs
                @warn "smooth_grid_par: no Newton's method solution for en flux, grid zone $i"
                gb_found = ux_new_energy[i] / √( c_cgs^2 - ux_new_energy[i-1]^2 )
            end

            ux_new_energy[i] = gb_found / √( 1 + gb_found^2 ) * c_cgs


            # Find the average downstream velocity, which will be used for
            # scaling the profiles. Only average the last 10 grid positions
            if i > (n_grid-10)
                ave_DwS_ux_px += ux_new_px[i]
                ave_DwS_ux_energy += ux_new_energy[i]
            end

        end  # loop over grid positions


        # Smooth the velocity profile, rescale it, and average the momentum and
        # energy curves as directed by user input
        smooth_profile(n_grid, ux_new_px)
        smooth_profile(n_grid, ux_new_energy)

        ave_DwS_ux_px /= 10
        ave_DwS_ux_energy /= 10

        ux_scale_fac = (u_Z - u_2)/(ux_new_px[1] - ave_DwS_ux_px)
        for i in 1:n_grid
            ux_new_px[i] = ux_scale_fac * (ux_new_px[i] - ave_DwS_ux_px)  +  u_2
            if x_grid_rg[i] ≥ 0
                ux_new_px[i] = u_2
            end
        end

        ux_scale_fac = (u_Z - u_2)/(ux_new_energy[1] - ave_DwS_ux_energy)
        for i in 1:n_grid
            ux_new_energy[i] = ux_scale_fac * (ux_new_energy[i] - ave_DwS_ux_energy)  +  u_2
            if x_grid_rg[i] ≥ 0
                ux_new_energy[i] = u_2
            end
        end

        for i in 1:n_grid
            ux_new[i] = (1 - smooth_mom_energy_fac) * ux_new_px[i] + smooth_mom_energy_fac * ux_new_energy[i]
        end

    end  # check on shock speed
    #-------------------------------------------------------------------------
    # Rel shock smoothing complete


    # Artificial smoothing if directed
    if x_art_start_rg < 0
        i_trans = findfirst(>(x_art_start_rg), x_grid_rg) - 1
        ux_scale_fac = -(ux_new[i_trans] - ux_new[n_grid]) / atan(x_grid_rg[i_trans])
        for i in i_trans:i_shock
            ux_new[i] = -atan(x_grid_rg[i]) * ux_scale_fac  +  ux_new[n_grid]
        end
    end


    # Average with previous profile
    # CHECKTHIS: what happens if γ * β is used as an averaging variable
    # instead of just u for a rel shock?
    for i in 1:n_grid
        ux_new[i] = (ux_new[i]  +  prof_wt_fac * ux_sk_grid[i])  /  ( 1 + prof_wt_fac )
    end


    # Compute output arrays based on new profile
    for i in 1:n_grid
        ux_sk_grid[i] = ux_new[i]
        γ_sf_grid[i] = 1 / √( 1 - (ux_sk_grid[i]/c_cgs)^2 )
        utot_grid[i] = ux_new[i]
        β_ef_grid[i] = (u_Z - ux_sk_grid[i]) / c_cgs /  (1 - u_Z*ux_sk_grid[i]/c_cgs^2 )
        γ_ef_grid[i] = 1 / √( 1 - β_ef_grid[i]^2 )

        # Include necessary corrections for turbulence compression
        z_comp       = (γ_Z * u_Z) / (γ_sf_grid[i] * ux_sk_grid[i])
        comp_fac     = √( 1/3  +  2/3 * z_comp^2 )
        comp_fac     = 1  +  (comp_fac - 1)*bturb_comp_frac
        btot_temp    = bmag_Z * comp_fac
        # Also include any additional amplification specified
        amp_fac      = 1  +  (btot_temp/bmag_Z - 1) * bfield_amp
        btot_grid[i] = bmag_Z  *  amp_fac

        # If a custom ε_B is in place, use that to calculate the magnetic field, not the preceding code.
        # Note that the R-H relations can be rearranged to read
        #     energy_density(x)  =  F_en0/u(x) - F_px0
        # assuming flux conservation everywhere.
        if use_custom_εB
            energy_density = (flux_energy_UpS + γ_Z*u_Z*density_Z*E₀_proton) / ux_sk_grid[i]  -  flux_px_UpS
            btot_grid[i] = √(8π * εB_grid[i] * energy_density)
        end
    end

    return
end

"""
Takes an input velocity profile and performs two tasks: enforces monotonicity by
smoothing out dips/bumps, and averages nearby points to smooth sharp edges

### Arguments

- y_prof: array holding velocity profile (modified in-place)
- n_grid: number of grid zones in position and velocity arrays
"""
function smooth_profile!(y_prof, n_grid)
    # Run from downstream to upstream and eliminate dips/bumps
    for i in n_grid:-2:-1
        if y_prof[i-1] < y_prof[i]
            y_prof[i-1] = y_prof[i]
        end
    end


    # Now smooth sharp edges by averaging adjacent grid locations;
    # handle the edges separately for different weighting
    y_prof_2 = similar(y_prof) # only 2:n_grid-1 actually assigned
    y_prof_2[2] = (2y_prof[1] + y_prof[2] + y_prof[3]) / 4
    for i in 3:n_grid-2
        y_prof_2[i] = (y_prof[i-1] + y_prof[i] + y_prof[i+1]) / 3
    end
    y_prof_2[n_grid-1] = ( y_prof[n_grid-2] + y_prof[n_grid-1] + 2*y_prof[n_grid] ) / 4

    for i in 2:n_grid-1
        y_prof[i] = y_prof_2[i]
    end
end
end # module
