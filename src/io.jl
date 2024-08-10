module io
using Printf: @sprintf
using ..constants: mp_cgs, c_cgs, pc2cm
using ..parameters: na_c, na_particles, na_ions, psd_max, na_grid, na_c

export print_input, print_plot_vals, tcut_print

"""
If tcuts were tracked during the run, print out the particle counts and spectra at each tcut

### Arguments
- wts_fileunit
- spectra_fileunit
- n_tcuts - Number of times to use in particle tracking
- tcuts - Array to hold cutoff times for particle tracking
- n_ions - Number of ion species
"""
function tcut_print(
        wts_fileunit, spectra_fileunit,
        n_tcuts, tcuts, n_ions,
        num_psd_mom_bins, psd_mom_bounds,
        i_iter, wt_coupled, spectra_coupled,
    )

    # Set a floor on the array values, and normalize spectra_coupled so that each spectrum
    #  has a total weight of 1 (actual weight, of course, is the entry in wt_coupled)
    for i_ion in 1:n_ions, i_cut in 1:n_tcuts
        if wt_coupled[i_cut,i_ion] < 1e-60
            wt_coupled[i_cut,i_ion] = 1e-99
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
    wts_fileunit["i_iter"] = i_iter
    wts_fileunit["tcuts_log"] = log10.(tcuts)           # 1     tcut times
    wts_fileunit["wt_coupled_log"] = log10.(wt_coupled) # 2-?   wts by species

    print_plot_vals(wts_fileunit)


    # Write out the spectra at each tcut, as a histogram. Split into groups
    # by species, and give each spectrum (i.e. tcut) its own column
    for i_ion in 1:n_ions
        # XXX how does hdf5 actually index these things? probably two levels of
        # indexing isn't even necessasry
        spectra_fileunit[i_ion]["psd_mom_bounds"] = psd_mom_bounds     # cgs units
        spectra_fileunit[i_ion]["psd_mom_bounds_nat"] = psd_mom_bounds - log10.(mp_cgs*c_cgs)  # nat units
        spectra_fileunit[i_ion]["spectra"] = log10.(spectra_coupled)

        print_plot_vals(spectra_fileunit)
    end


    # If this is the last iteration, close the files
    #if i_iter == n_itrs
    #    close(wts_fileunit)
    #    close(spectra_fileunit)
    #end
end

"""
Prints a whole mess of information to the screen and to mc_out

### Arguments

- n_pts_inj: target # of particles for injection distribution
- n_pts_pcut: target # of particle for low-E pcuts
- n_pts_pcut_hi: target # of particles for hi-E pcuts
- n_ions: # of ion species run
- num_psd_***_bins: # of bins in momentum/θ directions in psd
- n_xspec: # of additional locations to track particle spectra
- n_pcuts: # of pcuts
- n_grid: # of grid zone BOUNDARIES, including x_grid_stop
- rRH: Rankine-Hugoniot compression ratio for this shock
- r_comp: the compression ratio being used
- u/β/γ _Z/_2: total fluid velocity[cm/s, /c] and Lorentz factor for far UpS and DwS regions

### Returns
Nothing (it's all to the screen/file)
"""
function print_input(
        n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,
        num_psd_mom_bins, num_psd_θ_bins, n_xspec, n_pcuts, n_grid, rRH,
        r_comp, u_Z, β_Z, γ_Z, u_2, β_2, γ_2, denZ_ion, bmag_Z,
        bmag_2, θ_BZ, θ_B2, θ_u2, aa_ion, tZ_ion, tZ_electron,
        mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine, feb_UpS,
        feb_DwS, rg0, age_max, energy_pcut_hi, do_fast_push, bturb_comp_frac,
        sc_electron, outfileunit)


    # Print array parameters & usage
    n_pts_max = max(n_pts_inj, n_pts_pcut, n_pts_pcut_hi)
    parameter_str = @sprintf("""

Array parameters/usage:
   na_particles = %i
        na_ions = %i
        psd_max = %i
      n_pts_max = %i
         n_ions = %i
        psd_mom = %i
         θ_bins = %i

        na_grid = %i
           na_c = %i
        n_xpsec = %i
        n_pcuts = %i
         n_grid = %i

""",
        na_particles, na_ions, psd_max, n_pts_max, n_ions, num_psd_mom_bins, num_psd_θ_bins,
        na_grid, na_c, n_xspec, n_pcuts, n_grid)

    @info(parameter_str)
    println(outfileunit, parameter_str)


    # Shock speeds, escaping fluxes, and key densities
    shock_speeds_str = @sprintf("""

   rRH = %f
r_comp = %f

""", rRH, r_comp)

    @info(shock_speeds_str)
    println(outfileunit, shock_speeds_str)

    fluxes_str = @sprintf("""

     u_Z[cm/s] = %f    u_2[cm/s]       = %f
           β_Z = %f       β_2          = %f
           γ_Z = %f        γ_2         = %f
ρ_Z[prot/cm^3] = %f     ρ_2[prot/cm^3] = %f

""", u_Z, u_2, β_Z, β_2, γ_Z, γ_2, denZ_ion[1], denZ_ion[1]*γ_Z*β_Z/(γ_2*β_2))

    @info(fluxes_str)
    println(outfileunit, fluxes_str)

    # Relevant angles and field strengths
    magfield_str = @sprintf("""

    bmag_Z[G] = %f             bmag_2[G] = %f
      θ_BZ[°] = %f         θ_B2[°](calc) = %f
θ_u2[°](calc) = %f

""", bmag_Z, bmag_2, θ_BZ, θ_B2, θ_u2)
    @info(magfield_str)
    println(outfileunit, magfield_str)

    # Temperatures and Mach numbers
    if sc_electron
        _, i_electron = findmin(aa_ion)
        temp_electron = tZ_ion[i_electron]
    else
        temp_electron = tZ_electron
    end
    temp_str = @sprintf("""

Temp_Z[K](proton)   = %f
Temp_Z[K](electron) = %f

""", tZ_ion[1], temp_electron)

    @info(temp_str)
    println(outfileunit, temp_str)

    mach_str = @sprintf("""

Mach(sonic)  = %f
Mach(Alfven) = %f

""", mach_sonic, mach_alfven)

    @info(mach_str)
    println(outfileunit, mach_str)


    # Divisions of gyroperiod
    gyroperiod_str = @sprintf("""

N_g(coarse) = %f
N_g(fine)   = %f

""", xn_per_coarse, xn_per_fine)
    @info(gyroperiod_str)
    println(outfileunit, gyroperiod_str)


    # FEB info and max age
    feb_str = @sprintf("""

UpS FEB[rg0] = %f   UpS FEB[pc] = %f
DwS FEB[rg0] = %f   DwS FEB[pc] = %f
Max CR age[sec] = %f

""", feb_UpS/rg0, feb_UpS / pc2cm, feb_DwS/rg0, feb_DwS / pc2cm, age_max)
    @info(feb_str)
    println(outfileunit, feb_str)


    # Test-particle index from Keshet & Waxman (2005) [2005PhRvL..94k1102K] Eq 23
    kw_idx_str = @sprintf("  Keshet & Waxman (2005) index = %f",
                          ( 3β_Z  -  2*β_Z*β_2^2  +  β_2^3 ) / (β_Z - β_2))
    @info(kw_idx_str)
    println(outfileunit, kw_idx_str)


    # Energy to switch between low- and high-pcut particle counts
    pcut_str = @sprintf("  High pcut energy[keV/aa] = %f", energy_pcut_hi)
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
Prints the long list of floats at the end of each data set that are read in
by the plotting program to display information about the run.

### Arguments
- iunit: unit number to which the line will be written
"""
function print_plot_vals(
        iunit,
        n_pts_inj, n_pts_pcut, do_fast_push, inp_distr,
        dont_DSA, u_Z, γ_Z, r_comp, rRH, θ_BZ, θ_B2, θ_u2,
        bmag_Z, feb_UpS, rg0, Emax_keV, Emax_keV_per_aa, pmax_cgs,
        xn_per_coarse, xn_per_fine, mach_sonic, mach_alfven, x_grid_start_rg,
        x_grid_stop_rg, x_fast_stop_rg, η_mfp, x_art_start_rg, x_art_scale,
        feb_DwS, jet_rad_pc, jet_sph_frac, jet_dist_kpc, n_ions, aa_ion,
        zz_ion, denZ_ion, tZ_ion, smooth_mom_energy_fac, energy_inj,
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
           u_Z/1e5,                         # 1
           γ_Z,                             # 2
           r_comp,                          # 3
           rRH,                             # 4
           θ_BZ,                            # 5
           θ_B2,                            # 6
           θ_u2,                            # 7
           bmag_Z,                          # 8
           feb_UpS/rg0,                     # 9
           Emax_keV,                        # 10
           Emax_keV_per_aa,                 # 11
           pmax_cgs/(mp_cgs*c_cgs),         # 12
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
           feb_DwS/rg0,                     # 27
           jet_rad_pc,                      # 28
           jet_sph_frac,                    # 29
           jet_dist_kpc,                    # 30
           smooth_mom_energy_fac,           # 31
           x_in_distr,                      # 32
           energy_inj,                      # 33
           smooth_pressure_flux_psd_fac,    # 34
           x_DSA,                           # 35
           energy_transfer_frac,                # 36

           x_ions,
           aa_ion,
           zz_ion,
           denZ_ion,
           tZ_ion,
          )
         )

end
end # module io
