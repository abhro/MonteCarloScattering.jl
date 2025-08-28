using Printf: @sprintf
using Unitful, UnitfulAstro
using Unitful: cm
using UnitfulAstro: pc

export print_input, print_plot_vals, tcut_print

"""
    tcut_print(...)

If tcuts were tracked during the run, print out the particle counts and spectra at each tcut

### Arguments
FIXME
- `weights_fileunit`
- `spectra_fileunit`
- `n_tcuts`: Number of times to use in particle tracking
- `tcuts`: Array to hold cutoff times for particle tracking
- `n_ions`: Number of ion species
"""
function tcut_print(
        weights_fileunit, spectra_fileunit,
        n_tcuts, tcuts, n_ions,
        num_psd_mom_bins, psd_mom_bounds,
        i_iter, weight_coupled, spectra_coupled,
    )

    # Set a floor on the array values, and normalize spectra_coupled so that each spectrum
    # has a total weight of 1 (actual weight, of course, is the entry in weight_coupled)
    for i_ion in 1:n_ions, i_cut in 1:n_tcuts
        if weight_coupled[i_cut,i_ion] < 1e-60
            weight_coupled[i_cut,i_ion] = 1e-99
        end

        if sum(spectra_coupled[:,i_cut,i_ion]) > 1e-99
            spectra_coupled[:,i_cut,i_ion] = spectra_coupled[:,i_cut,i_ion] / sum(spectra_coupled[:,i_cut,i_ion])
        end

        for i_pt in 0:num_psd_mom_bins
            if spectra_coupled[i_pt,i_cut,i_ion] < 1e-60
                spectra_coupled[i_pt,i_cut,i_ion] = 1e-99
            end
        end
    end

    # Write out the particle counts at each tcut, as a histogram
    j_plot = 0
    weights_fileunit["i_iter"] = i_iter
    weights_fileunit["tcuts_log"] = log10.(tcuts)           # 1     tcut times
    weights_fileunit["weight_coupled_log"] = log10.(weight_coupled) # 2-?   weights by species

    print_plot_vals(weights_fileunit)


    # Write out the spectra at each tcut, as a histogram. Split into groups
    # by species, and give each spectrum (i.e. tcut) its own column
    for i_ion in 1:n_ions
        # XXX how does hdf5 actually index these things? probably two levels of
        # indexing isn't even necessasry
        spectra_fileunit[i_ion]["psd_mom_bounds"] = psd_mom_bounds     # cgs units
        spectra_fileunit[i_ion]["psd_mom_bounds_nat"] = psd_mom_bounds - log10.(mp*c)  # nat units
        spectra_fileunit[i_ion]["spectra"] = log10.(spectra_coupled)

        print_plot_vals(spectra_fileunit)
    end


    # If this is the last iteration, close the files
    #if i_iter == n_itrs
    #    close(weights_fileunit)
    #    close(spectra_fileunit)
    #end
end

"""
    print_input(...)

Prints a whole mess of information to the screen and to mc_out

### Arguments

- `n_pts_inj`: target number of particles for injection distribution
- `n_pts_pcut`: target number of particle for low-energy pcuts
- `n_pts_pcut_hi`: target number of particles for hi-energy pcuts
- `n_ions`: number of ion species run
- `num_psd_mom_bins`: number of bins in momentum directions in psd
- `num_psd_θ_bins`: number of bins in θ directions in psd
- `n_xspec`: number of additional locations to track particle spectra
- `n_pcuts`: number of pcuts
- `n_grid`: number of grid zone BOUNDARIES, including `x_grid_stop`
- `r_RH`: Rankine-Hugoniot compression ratio for this shock
- `r_comp`: the compression ratio being used
- `u/β/γ ₀/₂`: total fluid velocity[cm/s, /c] and Lorentz factor for far upstream and downstream regions

### Returns
Nothing (it's all to the screen/file)
"""
function print_input(
        n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,
        num_psd_mom_bins, num_psd_θ_bins, n_xspec, n_pcuts, n_grid, r_RH,
        r_comp, u₀, β₀, γ₀, u₂, β₂, γ₂, species, bmag₀,
        bmag₂, θ_B₀, θ_B₂, θ_u₂,
        mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine, feb_upstream,
        feb_downstream, rg₀, age_max, energy_pcut_hi, do_fast_push, bturb_comp_frac,
        outfileunit)


    # Print array parameters & usage
    n_pts_max = max(n_pts_inj, n_pts_pcut, n_pts_pcut_hi)
    parameter_str = @sprintf("""

Array parameters/usage:
   na_particles = %i
        psd_max = %i
      n_pts_max = %i
         n_ions = %i
        psd_mom = %i
         θ_bins = %i

           na_c = %i
        n_xpsec = %i
        n_pcuts = %i
         n_grid = %i

""",
        na_particles, psd_max, n_pts_max, n_ions, num_psd_mom_bins, num_psd_θ_bins,
        na_c, n_xspec, n_pcuts, n_grid)

    @info(parameter_str)
    println(outfileunit, parameter_str)


    # Shock speeds, escaping fluxes, and key densities
    shock_speeds_str = @sprintf("""

  r_RH = %f
r_comp = %f

""", r_RH, r_comp)

    @info(shock_speeds_str)
    println(outfileunit, shock_speeds_str)

    fluxes_str = ("""

u₀ = $u₀                 u₂ = $u₂
β₀ = $β₀                 β₂ = $β₂
γ₀ = $γ₀                 γ₂ = $γ₂
ρ₀ = $(density(species[1])) prot/cm³        ρ₂ = $(density(species[1])*γ₀*β₀/(γ₂*β₂)) prot/cm³

""")

    @info(fluxes_str)
    println(outfileunit, fluxes_str)

    # Relevant angles and field strengths
    magfield_str = ("""

bmag₀ = $(bmag₀)              bmag₂ = $(bmag₂)
θ_B₀ = $(θ_B₀)°          θ_B₂(calc) = $(θ_B₂)°
                     θ_u₂(calc) = $(θ_u₂)°

""")
    @info(magfield_str)
    println(outfileunit, magfield_str)

    # Temperatures and Mach numbers
    temp_str = ("""

T₀(proton)   = $(temperature(species[begin]))
T₀(electron) = $(temperature(species[end]))

""")

    @info(temp_str)
    println(outfileunit, temp_str)

    mach_str = @sprintf("""

Mach(sonic)  = %.3f
Mach(Alfven) = %.3f

""", mach_sonic, mach_alfven)

    @info(mach_str)
    println(outfileunit, mach_str)


    # Divisions of gyroperiod
    gyroperiod_str = @sprintf("""

N_g(coarse) = %i
N_g(fine)   = %i

""", xn_per_coarse, xn_per_fine)
    @info(gyroperiod_str)
    println(outfileunit, gyroperiod_str)


    # FEB info and max age
    feb_str = ("""

upstream FEB   = $(feb_upstream/rg₀) rg₀ = $(uconvert(pc, feb_upstream))
downstream FEB = $(feb_downstream/rg₀) rg₀ = $(uconvert(pc, feb_downstream))

Max CR age[s] = $(age_max)

""")
    @info(feb_str)
    println(outfileunit, feb_str)


    # Test-particle index from Keshet & Waxman (2005) [2005PhRvL..94k1102K] Eq 23
    kw_idx_str = @sprintf("  Keshet & Waxman (2005) index = %f",
                          (3β₀ - 2*β₀*β₂^2 + β₂^3) / (β₀ - β₂))
    @info(kw_idx_str)
    println(outfileunit, kw_idx_str)


    # Energy to switch between low- and high-pcut particle counts
    pcut_str = @sprintf("  High pcut energy = %f keV/aa", energy_pcut_hi)
    @info(pcut_str)
    println(outfileunit, pcut_str)


    # Finally a warning about possible complications later
    if do_fast_push && bturb_comp_frac > 0
        fast_push_warn_str = ("Both fast push and amplified B-field turbulence in use. "*
                              "Check flux equations in 'init_pop' for consistency")
        @warn(fast_push_warn_str)
        println(outfileunit, " WARNING: ", fast_push_warn_str)
    end
end

"""
    print_plot_vals(...)

Prints the long list of floats at the end of each data set that are read in
by the plotting program to display information about the run.

### Arguments
- `iunit`: unit number to which the line will be written
- TODO
"""
function print_plot_vals(
        iunit,
        n_pts_inj, n_pts_pcut, do_fast_push, inp_distr,
        dont_DSA, u₀, γ₀, r_comp, r_RH, θ_B₀, θ_B₂, θ_u₂,
        bmag₀, feb_upstream, rg₀, Emax, Emax_per_aa, pmax,
        xn_per_coarse, xn_per_fine, mach_sonic, mach_alfven, x_grid_start_rg,
        x_grid_stop_rg, x_fast_stop_rg, η_mfp, x_art_start_rg, x_art_scale,
        feb_downstream, jet_rad_pc, jet_sph_frac, jet_dist_kpc, n_ions, aa_ion,
        zz_ion, n₀_ion, T₀_ion, smooth_mom_energy_fac, energy_inj,
        smooth_pressure_flux_psd_fac, energy_transfer_frac, iseed_in
    )

    x_pts_inj   = float(n_pts_inj)
    x_pts_pcut  = float(n_pts_pcut)
    xseed       = float(iseed_in)
    x_fast_push = do_fast_push ? 66.0 : 0.0
    x_in_distr  = float(inp_distr)
    x_DSA       = dont_DSA ? 66.0 : 0.0
    x_ions      = float(n_ions)

    iannt = 3333
    idum  = 333

    # WARNING: these column numbers are reused in both subroutine read_old_prof
    # and the plotting program pg_color.f90. If they are ever changed,
    # modify the other codes accordingly!
    write(iunit,
          (
           iannt, idum,
           u₀/1e5,                          # 1
           γ₀,                              # 2
           r_comp,                          # 3
           r_RH,                            # 4
           θ_B₀,                            # 5
           θ_B₂,                            # 6
           θ_u₂,                            # 7
           bmag₀,                           # 8
           feb_upstream/rg₀,                # 9
           Emax,                            # 10
           Emax_per_aa,                     # 11
           pmax/(mp*c),                     # 12
           x_pts_inj,                       # 13
           x_pts_pcut,                      # 14
           xn_per_coarse,                   # 15
           xn_per_fine,                     # 16
           mach_sonic,                      # 17
           mach_alfven,                     # 18
           x_grid_start_rg,                 # 19
           xseed,                           # 20
           x_grid_stop_rg,                  # 21
           x_fast_push,                     # 22
           x_fast_stop_rg,                  # 23
           η_mfp,                           # 24
           x_art_start_rg,                  # 25
           x_art_scale,                     # 26
           feb_downstream/rg₀,              # 27
           jet_rad_pc,                      # 28
           jet_sph_frac,                    # 29
           jet_dist_kpc,                    # 30
           smooth_mom_energy_fac,           # 31
           x_in_distr,                      # 32
           energy_inj,                      # 33
           smooth_pressure_flux_psd_fac,    # 34
           x_DSA,                           # 35
           energy_transfer_frac,            # 36

           x_ions,
           aa_ion,
           zz_ion,
           n₀_ion,
           T₀_ion,
          )
         )

end
print_plot_vals(args...) = error("Fortran holdover function")
