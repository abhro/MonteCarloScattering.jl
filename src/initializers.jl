module initializers

import Random
using OffsetArrays
using LinearAlgebra: dot
using StaticArrays: @SVector

using ..parameters: na_grid, na_ions, β_rel_fl, psd_max, na_particles, energy_rel_pt, num_therm_bins
using ..constants: mp_cgs, me_cgs, kB_cgs, c_cgs, kB_cgs, keV2erg, E₀_proton

export init_pop, set_in_dist
export calc_DwS, calc_rRH
export set_psd_bins, set_photon_shells, setup_grid, setup_profile
export upstream_machs, upstream_fluxes

"""
Uses the Rankine-Hugoniot jump conditions to calculate the downstream conditions for a test
particle shock. Big difference between this subroutine and calc_rRH is that we already know
what the DwS speed is, courtesy of r_comp in the input.

### Arguments

- oblique: whether the shock is oblique
- θ_BZ: angle[deg] between UpS magnetic field and shock normal
- bmag_Z: far UpS magnetic field strength[Gauss]
- r_comp: compression ratio of shock
- β_Z: UpS bulk fluid speed, over c
- γ_Z: Lorentz factor associated with β_Z
- n_ions: number of different ion species
- aa_ion: array of particle species' atomic mass numbers
- zz_ion: array of particle species' charge numbers
- denZ_ion: array of far UpS densities for particle species
- tZ_ion: array of particle species' far UpS temperatures
- sc_electron: flag for whether electrons are a separate species
- tZ_electron: if electrons are not a separate species, this is their far UpS temperature

### Returns

- β: (total) bulk fluid speed DwS
- γ: Lorentz factor associated with β_2
- bmag: DwS magnetic field strength[Gauss]
- θ_B: angle[deg] between DwS magnetic field and shock normal
- θ_u: angle[deg] between DwS fluid velocity and shock normal
"""
function calc_DwS(oblique, bmag_Z, r_comp, β_Z)

    #--------------------------------------------------------------------------
    #  Possibility 1: Parallel at any shock speed
    #--------------------------------------------------------------------------
    if !oblique
        β    = β_Z / r_comp
        γ    = 1 / √( 1 - β^2 )
        bmag = bmag_Z
        θ_B  = 0.0
        θ_u  = 0.0

    #--------------------------------------------------------------------------
    #  Possibility 2: Oblique at any shock speed. Not currently supported by
    #  code, but included here in case code is extended in future.
    #--------------------------------------------------------------------------
    else
        error("ERROR in calc_DwS: not implemented for oblique shocks yet. ",
              "If this ever changes, don't forget to update the shock profile in subroutine 'setup_profile'.")
    end

    return β, γ, bmag, θ_B, θ_u
end

const rRH_maxitrs = 10000
const rRH_target_err = 1e-6
"""
Uses the Rankine-Hugoniot jump conditions to calculate the compression ratio for a shock
assuming test-particle conditions. In other words, (1) sharp shock, (2) negligible/no DSA,
and (3) no escaping flux. Additionally assumes that the inflowing plasma has
non-relativistic thermal speeds to make UpS adiabatic index exactly 5/3.

### Arguments
- β_Z, γ_Z: shock speed and Lorentz factor
- n_ions: number of different ion species
- aa_ion: array of particle species' atomic mass numbers
- zz_ion: array of particle species' charge numbers
- denZ_ion: array of far UpS densities for particle species
- tZ_ion: array of particle species' far UpS temperatures
- sc_electron: flag for whether electrons are a separate species
- tZ_electron: if electrons are not a separate species, this is their far UpS temperature
- oblique: controls whether to use parallel or oblique formulations of R-H relations

### Returns
- r_RH: Rankine-Hugoniot compression ratio
- γ_adiab_2_RH: ratio of specific heats (adiabatic index) for DwS region, assuming r_comp = r_RH
"""
function calc_rRH(
        β_Z, γ_Z, n_ions, aa_ion, zz_ion, denZ_ion,
        tZ_ion, sc_electron, tZ_electron, oblique)

    #--------------------------------------------------------------------------
    #  Four possibilities for R-H relations: nonrel/rel and parallel/oblique.
    #  Determine which of the four to use. Cutoff for nonrel/rel is set in module 'controls'
    #--------------------------------------------------------------------------
    relativistic = ( β_Z < β_rel_fl )

    # Calculate thermal pressure of far upstream gas
    pressure_Z   = dot(denZ_ion, tZ_ion) * kB_cgs
    ρ_Z          = dot(denZ_ion, aa_ion) * mp_cgs
    mask = (aa_ion .≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])

    # If electrons were not a separate species, add them in here
    if ! sc_electron
        pressure_Z +=  density_electron * kB_cgs*tZ_electron
        ρ_Z        +=  density_electron * me_cgs
    end

    #--------------------------------------------------------------------------
    #  Possibility 1: Nonrelativistic, parallel
    #  Solution comes from Ellison (1985) [1985JGR....90...29E].
    #  Uses far UpS Mach number to calculate r_RH.
    #--------------------------------------------------------------------------
    if !relativistic && !oblique

        # Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
        # to calculate the far UpS sound speed and Mach number     #assumecold
        γ_sph = 5//3
        c_s   = √( γ_sph * pressure_Z / ρ_Z )
        M_Z   = β_Z * c_cgs / c_s

        # Finally, use Equation (11) from Ellison (1985) to calculate r_RH.
        # Note that q = 0 here b/c we assume no escaping flux. This simplifies
        # the denominator quite a bit from the equation.
        r_RH = 8 / ( 2 + 6/M_Z^2 )

        # In non-rel case, downstream adiabatic index is pegged to 5/3
        γ_adiab_2_RH = 5//3


    #--------------------------------------------------------------------------
    #  Possibility 2: Relativistic, parallel
    #  Solution comes from Ellison+ (1990) [1991ApJ...378..214E].
    #  Uses relativistic Rankine-Hugoniot relations. See that paper for
    #  details of equations and associated quantities. Briefly,
    #     R-H1:         g₀ n₀ b₀  =  g₂ n₂ b₂
    #     R-H2:  g₀² w₀ b₀² + P₀  =  g₂² w₂ b₂² + P₂
    #     R-H3:  g₀² w₀ b₀        =  g₂² w₂ b₂
    #  where
    #     w    = E_rm + E_ke + P,   <--- enthalpy as total energy density + pressure
    #     E_rm = n m c²             <--- rest mass energy density
    #     E_ke = n m c² e(p)        <--- kinetic energy density, with e(p) =  √( 1 + (p/mc)² )  -  1
    #     P    = ⅓ n p v            <--- pressure
    #
    #  Assumes that downstream particle distributions are δ-functions.
    #  Solves for p2 using Newton's method, then works backwards to r_RH.
    #-------------------------------------------------------------------------
    elseif relativistic && !oblique

        # Calculate two quantities to be used during loop to find r_RH: the
        # rest mass-energy of each species, and the density relative to protons
        rm_ion          = mp_cgs .* aa_ion .* c_cgs^2 # isa Vector
        density_rel_ion = denZ_ion ./ denZ_ion[1]

        # Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
        # to calculate the far UpS enthalpy
        # #assumecold
        γ_sph = 5//3
        w_Z   = ρ_Z * c_cgs^2  +  γ_sph/(γ_sph - 1) * pressure_Z

        # Calculate the far UpS momentum flux
        UpS_mom_flux = γ_Z^2 * w_Z * β_Z^2  +  pressure_Z
        UpS_num_flux = γ_Z * denZ_ion[1] * β_Z # Protons only here; not strictly correct but appropriate for later use


        # TODO : find the equation then plug it into a root finder
        function F(p)
            # TODO define m, c, E, f_n (=ups_num_flux), f_p (=ups_mom_flux)
            p_dimless = p / (m * c)
            γ = √(1 + p_dimless^2)
            pressure = E/3 * p_dimless^2/γ
            w        = E*(γ + p_dimless^2/3γ)
            return f_n/γ * (w*γ^2 + pressure) - f_p
        end

        # Now use Newton's method to determine the downstream momentum that
        # satisfies the R-H relations.
        # Assumptions: (1) momentum distribution functions are δ-functions
        # rather than thermal, and (2) any non-proton species have p ∝ m
        #------------------------------------------------------------------------
        # We know a priori that the R-H compression ratio will be between 3 and 4.
        # Use a compression ratio of 4.5 as an upper bound, which sets a minimum
        # value for γ_2 and in turn an upper limit on the downstream momentum.
        γ_2_min = 1 / √( 1  -  (β_Z/4.5)^2 )
        w_fac_max = γ_Z * w_Z / (denZ_ion[1] * γ_2_min)

        relative_ion_energy = dot(rm_ion, density_rel_ion)

        # In the following quadratic equation, we'll use w_fac_max/rm_avg in the
        # coefficients. So calculate rm_avg here
        rm_avg = relative_ion_energy
        # Handle possibility that electrons aren't included self-consistently
        if ! sc_electron
            rm_avg += E₀_electron * density_electron/denZ_ion[1]
        end

        p2_max_A  = 16//9
        p2_max_B  = 8//3 - (w_fac_max/rm_avg)^2
        p2_max_C  = 1 - (w_fac_max/rm_avg)^2
        p2_max_sq = (-p2_max_B + √(p2_max_B^2 - 4*p2_max_A*p2_max_C)) / 2p2_max_A
        p2_max    = √p2_max_sq * mp_cgs*c_cgs

        # Initial guess for downstream proton momentum, as well as nearby
        # location to use for finite difference approximation
        p2_guess   = p2_max / 1.001
        p2_guess_p = p2_max
        Δp2_guess  = p2_guess_p - p2_guess

        i = 1
        for outer i in 1:rRH_maxitrs

            # Calculate downstream momentum fluxes associated with both p2_guess
            # and p2_guess_p
            # 1. Calculate enthalpy, which is summed over all particle species
            # 2. Calculate pressure similarly

            p2_o_mc = p2_guess / (mp_cgs * c_cgs)  # Protons only here because of how the math in P_fac & w_fac works out
            lorentz_p2 = 1/√(1 + p2_o_mc^2)

            # Pressure with proton density factored out
            P_fac = relative_ion_energy * lorentz_p2/3 * p2_o_mc^2
            # Enthalpy with proton density factored out
            w_fac = relative_ion_energy * ( 1/lorentz_p2 + lorentz_p2/3 * p2_o_mc^2 )

            # Handle possibility that electrons aren't included self-consistently
            if ! sc_electron
                P_fac += (density_electron / denZ_ion[1]) * E₀_electron * lorentz_p2/3 * p2_o_mc^2
                w_fac += (density_electron / denZ_ion[1]) * E₀_electron * ( 1/lorentz_p2 + lorentz_p2/3 * p2_o_mc^2 )
            end

            γ_2 = γ_Z * w_Z / (denZ_ion[1] * w_fac)  # Using only proton density is correct
            γβ_2 = √( γ_2^2 - 1 )

            F_p2_guess = UpS_num_flux * (w_fac * γβ_2 + P_fac / γβ_2)  -  UpS_mom_flux

            #------------------------------------------

            p2_o_mc = p2_guess_p / (mp_cgs * c_cgs)  # Protons only here because of how the math in P_fac & w_fac works out

            # Pressure with proton density factored out
            P_fac = relative_ion_energy * 1/3 * p2_o_mc^2 / √( 1  +  p2_o_mc^2 )
            # Enthalpy with proton density factored out
            w_fac = relative_ion_energy * ( √( 1 + p2_o_mc^2 ) + 1/3 * p2_o_mc^2 / √(1 + p2_o_mc^2 ) )

            # Handle possibility that electrons aren't included self-consistently
            if ! sc_electron
                P_fac += (density_electron / denZ_ion[1]) * E₀_electron * 1/3 * p2_o_mc^2 / √( 1 + p2_o_mc^2 )
                w_fac += (density_electron / denZ_ion[1]) * E₀_electron * ( √( 1 + p2_o_mc^2 ) + 1/3 * p2_o_mc^2 / √( 1 + p2_o_mc^2 ) )
            end

            γ_2 = γ_Z * w_Z / (denZ_ion[1] * w_fac)  # Using only proton density is correct
            γβ_2 = √( γ_2^2 - 1 ) # = β² / (1 - β²)

            F_p2_guess_p = UpS_num_flux * (w_fac*γβ_2 + P_fac/γβ_2) - UpS_mom_flux
            #------------------------------------------

            # Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
            Fprime_p2_guess = (F_p2_guess_p - F_p2_guess) / Δp2_guess

            # Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
            p2_guess_next = p2_guess  -  F_p2_guess/Fprime_p2_guess

            # Relative change in this step
            err_curr = (p2_guess_next - p2_guess) / p2_guess

            # If the relative change is small enough, we've found our solution and
            # can exit the loop; otherwise return for another cycle
            if abs(err_curr) < rRH_target_err
                p2_found = p2_guess_next
                break
            end

            # Make sure new value for p2_guess is less than p2_max.
            # Use a weighted average for p2_guess to converge faster.
            if p2_guess_next*1.001 ≥ p2_max
                p2_guess_p   = 0.2 * (p2_guess  +  4p2_max)
                p2_guess     = p2_guess_p / 1.001
            else
                p2_guess     = p2_guess_next
                p2_guess_p   = 1.001 * p2_guess
            end
            Δp2_guess = p2_guess_p - p2_guess

        end

        # Did we hit the maximum number of iterations without finding the flux-conserving solution?
        i ≥ rRH_maxitrs && error("ERROR in calc_rRH: Newton method did not find solution")


        # Calculate the compression ratio β_Z/β_2 associated with p2_found
        p2_o_mc = p2_found / (mp_cgs * c_cgs)  # Protons only here because of how the math in w_fac works out

        # Pressure, internal energy, and enthalpy with proton density factored out
        P_fac = relative_ion_energy * 1/3 * p2_o_mc^2 / √( 1  +  p2_o_mc^2 )
        e_fac = relative_ion_energy * ( √( 1 + p2_o_mc^2 ) - 1 )
        w_fac = relative_ion_energy * ( √( 1 + p2_o_mc^2 ) + 1/3 * p2_o_mc^2 / √(1 + p2_o_mc^2 ) )

        # Handle possibility that electrons aren't included self-consistently
        if !sc_electron
            P_fac += (density_electron / denZ_ion[1]) * E₀_electron * 1/3 * p2_o_mc^2 / √( 1 + p2_o_mc^2 )
            e_fac += (density_electron / denZ_ion[1]) * E₀_electron * ( √( 1 + p2_o_mc^2 ) - 1 )
            w_fac += (density_electron / denZ_ion[1]) * E₀_electron * ( √( 1 + p2_o_mc^2 ) + 1/3 * p2_o_mc^2 / √( 1 + p2_o_mc^2 ) )
        end

        # Calculate adiabatic index downstream
        γ_adiab_2_RH = 1 + P_fac/e_fac

        # Finally, get downstream speed and compression ratio
        β_2 = √( 1 - (denZ_ion[1] * w_fac/γ_Z * w_Z)^2 )  # Using only proton density is correct


        r_RH = β_Z/β_2
        #------------------------------------------------------------------------
        # r_RH found using Newton's method

    #--------------------------------------------------------------------------
    #  Possibility 3: Oblique at any shock speed. Not currently supported by
    #  code, but included here in case code is extended in future.
    #--------------------------------------------------------------------------
    else
        error("ERROR in calc_rRH: not implemented for oblique shocks yet. ",
              "If this ever changes, don't forget to update the shock profile in subroutine 'setup_profile'.")
    end

    return r_RH, γ_adiab_2_RH
end


"""
Sets the BOUNDARIES of the bins of the phase space distribution. The bins are numbered from
0 to num_psd_***_bins, each boundary denotes the lower edge of that # bin; the indices thus
run from 0 to num_psd_***_bins + 1.

!!! warning
    for angles, the number stored in psd_θ_bounds increases at first (increasing θ),
    then decreases because increasing the angle means decreasing the cosine.

- For total momentum: logarithmic spacing over all decades from Emin_keV to Emax_keV.
- For angle: linear spacing of cosine values for angles between psd_θ_fine and π.
  Below psd_θ_fine spacing is logarithmic in θ for some # of decades down to psd_θ_min.
- For both: values less than the minimum are equivalent to 0.0.

### Arguments

- psd_mom_min: minimum momentum[cgs] to use in PSD
- psd_mom_max: maximum momentum[cgs] to use in PSD
- psd_bins_per_dec_***: # of bins per decade to use in logarithmically-spaced regions of PSD
- psd_lin_cos_bins: # of bins to divide range [-1,cos(psd_θ_fine)] into
- psd_cos_fine: cutoff between lin/cos and log/θ spacing of PSD bins
- psd_θ_min: minimum angle for PSD

### Returns
- num_psd_***_bins: total number of bins along given dimension, not counting bin 0
- Δcos: size of each linear cosine bin
- psd_***_bounds: boundaries between bins, and upper edge of final bin
"""
function set_psd_bins(
        psd_mom_min, psd_mom_max, psd_bins_per_dec_mom, psd_bins_per_dec_θ,
        psd_lin_cos_bins, psd_cos_fine, psd_θ_min)

    num_psd_mom_bins, psd_mom_bounds = set_psd_mom_bins(
        psd_mom_min, psd_mom_max, psd_bins_per_dec_mom) # Set the momentum bins
    num_psd_θ_bins, Δcos, psd_θ_bounds = set_psd_angle_bins(
        psd_bins_per_dec_θ, psd_lin_cos_bins, psd_cos_fine, psd_θ_min) # Set angle bins

    return (num_psd_mom_bins, num_psd_θ_bins, Δcos, psd_mom_bounds, psd_θ_bounds)
end

function set_psd_mom_bins(psd_mom_min, psd_mom_max, psd_bins_per_dec_mom)
    num_psd_mom_bins = trunc(Int, log10(psd_mom_max / psd_mom_min) *  psd_bins_per_dec_mom)
    num_psd_mom_bins += 2   # Add two extra bins just to be safe
    if (num_psd_mom_bins+1) > psd_max
        error("More PSD momentum bins needed than allowed by psd_max. Bins required: ", num_psd_mom_bins)
    end

    # Fill in the array psd_mom_bounds, remembering that the array holds LOWER
    # boundaries of that bin
    psd_mom_bounds = zeros(0:psd_max)
    psd_mom_bounds[0] = -99.0
    psd_mom_bounds[1:num_psd_mom_bins+1] .= range(start  = log10(psd_mom_min),
                                                  step   = 1/psd_bins_per_dec_mom,
                                                  length = num_psd_mom_bins+1)

    return num_psd_mom_bins, psd_mom_bounds
end

function set_psd_angle_bins(psd_bins_per_dec_θ, psd_lin_cos_bins, psd_cos_fine, psd_θ_min)
    psd_θ_fine = acos(psd_cos_fine)
    ten_root_θ = exp10(1 / psd_bins_per_dec_θ)

    psd_log_θ_bins = trunc(Int, log10(psd_θ_fine/psd_θ_min) * psd_bins_per_dec_θ)
    num_psd_θ_bins = psd_log_θ_bins  +  psd_lin_cos_bins

    if (num_psd_θ_bins+1) > psd_max
        error("More PSD anglular bins needed than allowed by psd_max. Bins required: ", num_psd_θ_bins)
    end

    # Fill the logarithmic part of psd_θ_bounds using the angle (in radians), NOT its logarithm
    psd_θ_bounds = zeros(0:psd_max)
    psd_θ_bounds[0] = 1e-99
    psd_θ_bounds[1] = psd_θ_min
    for i in 2:psd_log_θ_bins
        psd_θ_bounds[i] = psd_θ_bounds[i-1] * ten_root_θ
    end
    psd_θ_bounds[psd_log_θ_bins+1:end] .= 0.0

    # Now fill in the linear part of psd_θ_bounds.
    # Note that the lower boundary of the first cell is psd_cos_fine
    Δcos = (psd_cos_fine + 1) / psd_lin_cos_bins

    for i in 1:psd_lin_cos_bins+1
        psd_θ_bounds[psd_log_θ_bins+i] = psd_cos_fine - Δcos*(i-1)
    end

    return num_psd_θ_bins, Δcos, psd_θ_bounds
end

"""
If photon calculation is desired, photons will be collected into UpS and DwS shells for
easier viewing. This subroutine sets the endpoints of the shells, as well as their midpoints.

Because the particle spectrum changes most rapidly near the shock, zones should be small
near the shock and get larger as you move further away.

First, divide the domain between the shock and the FEB (either UpS or DwS into n sections,
where n is the respective number of shells to use. HOWEVER, do this by exponent, ranging
from -1 to log10(|x_FEB|).

Keep track of the boundaries of each shells, since we will need those for calculating the
total number of particles emitting when we get to that point in photon production.
The midpoints are less useful, since photons are calculated on a zone-by-zone basis
rather than just at select points in the shock profile. Keep them anyway, since the
memory overhead is low.
"""
function set_photon_shells(
        num_UpS_shells, num_DwS_shells,
        # the Fortran subroutine gets these arguments from the controls module
        use_prp, feb_UpS, feb_DwS, rg0, x_grid_stop_rg,
    )

    x_shell = zeros(na_grid)
    x_shell_end_points = zeros(na_grid)

    # Handle UpS shells first
    x_section_width = (log10(abs(feb_UpS/rg0))+1) / num_UpS_shells
    for i in 1:num_UpS_shells
        # Calculate UpS and DwS endpoints of each region, as well as midpoint in log space
        if i == 1
            # Special case when i = 1, since our region starts at the shock.
            x_region_start = 0.0
            x_region_end   = exp10(-1 + x_section_width)
            x_region_mid   = exp10(-1 + x_section_width/2)
        else
            # In the general case, note that x_region_start should be the same as
            # the previous region's x_region_end. This can be checked with print
            # or write statements at runtime.
            x_region_start = exp10(-1 + x_section_width * (i-1) )
            x_region_end   = exp10(-1 + x_section_width *  i    )
            x_region_mid   = exp10(-1 + x_section_width * (i - 1/2) )
        end

        # Update the arrays with this information, remembering that the eventual array will
        # count downstream from the UpS FEB (so some array index juggling is necessary)
        # Also, add in the factor of -1 here, since upstream coordinates should
        # be negative in the MC code.
        x_shell[num_UpS_shells+1 - i]                = -x_region_mid
        x_shell_end_points[num_UpS_shells+1 - i]     = -x_region_end
        x_shell_end_points[num_UpS_shells+1 - i + 1] = -x_region_start
    end

    # And repeat the process for the downstream shells. The DwS limit is set
    #  differently if using a PRP or a FEB
    x_section_width = (log10(use_prp ? x_grid_stop_rg : feb_DwS/rg0)+1) / num_DwS_shells

    for i in 1:num_DwS_shells
        # Calculate UpS and DwS endpoints of each region, as well as midpoint in log space
        if i == 1
            # Special case when i = 1, since our region starts at the shock.
            x_region_start = 0.0
            x_region_end   = exp10(-1 + x_section_width)
            x_region_mid   = exp10(-1 + x_section_width/2)
        else
            # In the general case, note that x_region_start should be the same as
            # the previous region's x_region_end. This can be checked with print
            # or write statements at runtime.
            x_region_start = exp10(-1 + x_section_width * (i-1) )
            x_region_end   = exp10(-1 + x_section_width *  i    )
            x_region_mid   = exp10(-1 + x_section_width * (i - 1/2) )
        end

        # Update the arrays with the information. Less index juggling here
        x_shell[num_UpS_shells + i]                = x_region_mid
        x_shell_end_points[num_UpS_shells + i]     = x_region_start
        x_shell_end_points[num_UpS_shells + i + 1] = x_region_end
    end


    # Convert from units of rg0 to cm
    for i in 1:(num_UpS_shells + num_DwS_shells + 1)
        x_shell_end_points[i] = x_shell_end_points[i] * rg0
    end
    return (x_shell, x_shell_end_points)
end

function setup_grid(
        outfile,
        # the Fortran subroutine gets these arguments from the controls module
        x_grid_start_rg, use_prp, feb_DwS, x_grid_stop_rg, rg0)

    # Recall that rg0 is the gyroradius of a proton with speed u_Z in magnetic field bmag_Z.

    # Set the start and stop positions in units of rg0
    x_grid_start = x_grid_start_rg * rg0
    if !use_prp
        x_grid_stop_rg = feb_DwS / rg0
        println(outfile,
                "DwS FEB set at x = $x_grid_stop_rg rg0. Overwriting entered value for 'XGDDW'.")
        println(outfile)
    end
    x_grid_stop = x_grid_stop_rg * rg0

    # Logarithmically-spaced grid zones run from x_grid_start_rg to -10 rg0. Set them here.
    n_log_UpS = 27
    Δlog = (log10(-x_grid_start_rg) - 1) / n_log_UpS-1

    x_grid_rg = zeros(0:na_grid)
    for i in 1:n_log_UpS
        x_grid_rg[i] = -exp10( log10(-x_grid_start_rg) - (i-1)*Δlog )
    end
    # Many grid zones are set manually; zones can easily be added/removed, but
    # make sure to change the number in the log-spaced regions UpS or DwS
    first_zone = @SVector [-9.0, -8.0, -7.0, -6.0, -5.0, -4.5, -4.0, -3.5, -3.0,
                           -2.5, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0,
                           -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2,
                           -0.15, -0.1,
                           -0.07, -0.05, -0.04, -0.03, -0.02, -0.015, -0.01,
                           -3e-3, -1e-3]
    i_grid_ct = n_log_UpS + 1
    n = length(first_zone)
    x_grid_rg[i_grid_ct:(i_grid_ct+n-1)] .= first_zone
    i_grid_ct += n

    # Extremely fine spacing right around the shock
    extremely_fine_spacing = @SVector [-1e-4, -1e-7, 0.0, 1e-7, 1e-4]
    n = length(extremely_fine_spacing)
    i_grid_ct += 1
    x_grid_rg[i_grid_ct:i_grid_ct+n-1] .= extremely_fine_spacing
    i_grid_ct += n

    # Downstream from the shock, spacing doesn't need to be quite so fine
    # because velocity gradients aren't as extreme, if they exist at all
    downstream_spacing = @SVector [1e-3, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 0.10,
                                   0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    i_grid_ct += 1
    n = length(downstream_spacing)
    x_grid_rg[i_grid_ct:i_grid_ct+n-1] .= downstream_spacing
    i_grid_ct += n

    # As seen above, the manually-set grid zones end at x = +1 rg0.
    # DwS from there, more log-spaced zones.
    n_log_DwS = 16
    x_end_man = x_grid_rg[i_grid_ct]
    Δlog      = (log10(x_grid_stop_rg) - log10(x_end_man)) / n_log_DwS

    for i in 1:n_log_DwS
        x_grid_rg[i_grid_ct+i] = exp10(log10(x_end_man) + i*Δlog)
    end

    i_grid_ct += n_log_DwS


    # Set n_grid, and the extreme boundaries of the grid
    n_grid = i_grid_ct
    x_grid_rg[0] = -1e30
    x_grid_rg[n_grid+1] = 1e30

    # Convert everything from rg0 units to cgs units
    x_grid_cm = zeros(0:na_grid)
    x_grid_cm[0:n_grid+1] .= x_grid_rg[0:n_grid+1] .* rg0

    x_grid_cm[n_grid+1] = 1e30 * rg0

    return (n_grid, x_grid_start, x_grid_stop, x_grid_rg, x_grid_cm)
end

"""
Calculates the far upstream fluxes for the shock.

Two different cases considered:

1. Non-relativistic oblique shock. Uses equations of Ellison+ (1996) [1996ApJ...473.1029E]
2. Relativistic shock, any obliquity. Uses equations of Double+ (2004) [2004ApJ...600..485D]

Only oblique equations used because they reduce trivially to parallel cases when θ_BZ = 0.

HOWEVER, assumes that z-component of far UpS velocity is 0 in all cases; in practice oblique
shocks would induce some z-velocity in the shock profile even though particles initially
arrive with no bulk z component. Also assumes isotropic initial pressure, so no off-diagonal
components in pressure tensor.

### Arguments

FIXME

No inputs; pulls everything from module 'controls'

### Returns
- flux_px_UpS: far UpS momentum flux, x component
- flux_pz_UpS: far UpS momentum flux, z component
- flux_energy_UpS: far UpS energy flux
"""
function upstream_fluxes(
        # the Fortran subroutine gets these arguments from the controls module
        oblique, n_ions, denZ_ion, tZ_ion, aa_ion, zz_ion, sc_electron, tZ_electron,
        bmag_Z, θ_BZ, γ_Z, β_Z, u_Z
    )

    oblique && error("Cannot handle oblique shocks yet")

    # UpS internal energy density and pressure, assuming isotropic particle
    # distribution. Note that this INCLUDES the mass-energy density, which
    # is typically omitted in non-rel calculations
    pressure_Z   = dot(denZ_ion, tZ_ion) * kB_cgs
    ρ_Z          = dot(denZ_ion, aa_ion) * mp_cgs
    mask = (aa_ion .≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])

    # If electrons were not a separate species, add them in here
    if ! sc_electron
        pressure_Z += density_electron * kB_cgs*tZ_electron
        ρ_Z        += density_electron * me_cgs
    end

    # Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
    # to calculate the far UpS internal energy
    # #assumecold
    γ_sph = 5//3
    e_Z   = ρ_Z * c_cgs^2  +  1/(γ_sph - 1) * pressure_Z

    # Quantities related to the UpS magnetic field. Note that B_z is the
    # z-component of the magnetic field, not B_0
    B_x = bmag_Z * cosd(θ_BZ)
    B_z = bmag_Z * sind(θ_BZ)


    # Momentum, x-component
    F_px_fl = (γ_Z*β_Z)^2 * (e_Z + pressure_Z)  +  pressure_Z   # Fluid part (Double+ Eq 23)
    F_px_EM = (γ_Z^2*β_Z^2*bmag_Z^2 + γ_Z^2*(B_z^2-B_x^2)) / 8π # EM part (Double+ Eq 25)
    flux_px_UpS = F_px_fl  +  F_px_EM                           # Total

    # Momentum, z-component (Fluid Part 0, from Double+ Eq 24)
    # Total = EM part (Double+ Eq 26)
    flux_pz_UpS = -γ_Z/4π * B_x * B_z

    # Energy
    F_energy_fl = γ_Z^2 * β_Z * (e_Z + pressure_Z) # Fluid part (Double+ Eq 20)
    F_energy_EM = γ_Z^2 * β_Z * B_z^2/4π           # EM part (Double+ Eq 21)
    # Total -- convert to cgs units!
    flux_energy_UpS = c_cgs * (F_energy_fl + F_energy_EM)

    # And subtract off the mass-energy flux to bring it in line with non-rel
    # calculations and what the MC code actually tracks
    flux_energy_UpS -= γ_Z * u_Z * ρ_Z*c_cgs^2

    # Non-relativistic version. Note that it's missing the ρc² flux present in the
    # relativistic forms above. It is also expanded to second order in β_Z (only in the
    # hydro terms, for now) to allow for more precise matching with the relativistic version
    if β_Z < β_rel_fl
        flux_px_UpS = ρ_Z * u_Z^2 * (1 + β_Z^2) +  pressure_Z * (1 + γ_sph/(γ_sph-1)*β_Z^2) +  B_z^2/8π
        flux_pz_UpS = - B_x * B_z / 4π
        flux_energy_UpS = ρ_Z * u_Z^3 * (1 + 1.25*β_Z^2)/2 + pressure_Z * u_Z * γ_sph/(γ_sph-1) * (1+β_Z^2) + u_Z*B_z^2/4π
    end

    return (flux_px_UpS, flux_pz_UpS, flux_energy_UpS)
end

"""
Calculates the sonic and Alfven mach numbers for the shock.

For speed of sound, uses Equation (13) of Fujimura & Kennel (1979) [1979A%26A....79..299F]
For Alfven wave speed, uses Equation (46) of Gedalin (1993) [1993PhRvE..47.4354G]

### Arguments

FIXME

No input arguments; these come from module "controls"

### Returns
- mach_sonic: sonic Mach number of UpS flow
- mach_alfven: Alfvenic Mach number of UpS flow
"""
function upstream_machs(
        β_Z, n_ions, denZ_ion, tZ_ion, aa_ion, zz_ion, sc_electron, tZ_electron, bmag_Z)

    # Assume cold UpS plasma, so that the adiabatic index is 5/3 identically
    # #assumecold
    γ_adiab = 5//3


    # Find FK1979's R factor, the ratio of pressure to rest mass energy density
    pressure_Z = dot(denZ_ion, tZ_ion) * kB_cgs
    ρ_Z        = dot(denZ_ion, aa_ion) * mp_cgs
    mask = (aa_ion .≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])

    # If electrons were not a separate species, add them in here
    if ! sc_electron
        pressure_Z +=  density_electron * kB_cgs*tZ_electron
        ρ_Z        +=  density_electron * me_cgs
    end

    R_fac = pressure_Z / (ρ_Z * c_cgs^2)


    # Calculate the sound speeds differently based on whether the shock is non-rel or rel
    if β_Z < β_rel_fl # non-relativistic
        # Use standard forms for both speeds
        c_s = √( γ_adiab * pressure_Z / ρ_Z )
        v_A = bmag_Z / √( 4π * ρ_Z )

    else

        # Plug everything into FK1979's Equation (13)
        # cₛ²/c² = ΓP/w̄ = ΓR(aR + 1)
        a_fac = γ_adiab / (γ_adiab - 1)
        c_s   = c_cgs * √( γ_adiab * R_fac / (a_fac*R_fac + 1) )

        # And into Gedalin (1993)'s Equation (46); note assumption that
        # equation of state is    e = ρc² + P/(Γ-1)
        #     v_A² = (B²/4π) / (ε + p + B²/4π)..................Gedalin Eq. 46
        enthalpy = a_fac * pressure_Z  +  ρ_Z * c_cgs^2
        v_A = c_cgs  /  √( 1  +  4π * enthalpy / bmag_Z^2 )

    end

    # Now calculate the Mach numbers using the derived wave speeds
    mach_sonic  = β_Z * c_cgs / c_s
    mach_alfven = β_Z * c_cgs / v_A

    return (mach_sonic, mach_alfven)
end

"""
Sets the initial values of the shock profile

### Arguments

TODO

### Returns
- ux_sk_grid: bulk fluid velocity along x axis (i.e., perpendicular to shock face) in shock frame
- uz_sk_grid: bulk fluid velocity along z axis (i.e., parallel to shock face) in shock frame
- utot_grid: total bulk fluid velocity in shock frame
- γ_sf_grid: bulk flow Lorentz factor in shock frame
- β_ef_grid: relative x-axis speed between plasma and explosion frames
- γ_ef_grid: Lorentz factor associated with β_ef_grid
- btot_grid: total magnetic field
- θ_grid: angle of magnetic field[radians] relative to shock normal (i.e., to x axis)
- εB_grid: user-defined function for fraction of energy density in magnetic field.
  Sets value of btot_grid
- bmag_2: field strength in DwS region. Initially set in calc_DwS, it may be reset here
  depending on values of bturb_comp_frac & bfield_amp
"""
function setup_profile(
        u_Z, β_Z, γ_Z, bmag_Z, θ_BZ, r_comp,
        bturb_comp_frac, bfield_amp, use_custom_εB, n_ions, aa_ion,
        denZ_ion, sc_electron, zz_ion, flux_px_UpS, flux_energy_UpS,
        n_grid, x_grid_cm, x_grid_rg,
    )

    ux_sk_grid = zeros(0:n_grid)
    uz_sk_grid = zeros(0:n_grid)
    utot_grid = zeros(0:n_grid)
    γ_sf_grid = zeros(0:n_grid)
    β_ef_grid = zeros(0:n_grid)
    γ_ef_grid =  ones(0:n_grid)
    btot_grid = zeros(0:n_grid)
    θ_grid    = zeros(0:n_grid)

    comp_fac = 0.0
    for i in 0:n_grid
        if x_grid_cm[i] < 0
            ux_sk_grid[i] = u_Z
            γ_sf_grid[i] = γ_Z
            #β_ef_grid[i] = 0.0
            #γ_ef_grid[i] = 1.0
            btot_grid[i] = bmag_Z
        else
            ux_sk_grid[i] = u_Z / r_comp
            γ_sf_grid[i] = 1 / √( 1 - (utot_grid[i]/c_cgs)^2 )
            β_ef_grid[i] = (β_Z - ux_sk_grid[i]/c_cgs) /  ( 1 - β_Z*ux_sk_grid[i]/c_cgs )
            γ_ef_grid[i] = 1 / √( 1 - β_ef_grid[i]^2 )

            # When initializing magnetic field, include necessary corrections for turbulence compression
            z_comp         = (γ_Z * u_Z) / (γ_sf_grid[i] * ux_sk_grid[i])
            local comp_fac = 1  +  (√(( 1  +  2*z_comp^2 )/3) - 1)*bturb_comp_frac
            # Also include any additional amplification specified
            amp_fac        = 1  +  (comp_fac - 1) * bfield_amp
            btot_grid[i]   = bmag_Z  *  amp_fac

        end
        uz_sk_grid[i]  = 0.0
        utot_grid[i]  = ux_sk_grid[i] # uz_sk_grid[i] is 0, so don't bother adding
        θ_grid[i] = deg2rad(θ_BZ)
    end

    εB_grid   = zeros(0:n_grid)
    # If directed in data_input, use a custom-defined ε_B to set btot_grid.
    #-------------------------------------------------------------------------------
    if use_custom_εB

        # Calculate ε_B0, which depends on far UpS magnetic field and mass density.
        # If electrons aren't a separate species, they don't contribute enough mass
        # to be important
        density_Z = dot(denZ_ion, aa_ion)
        εB0 = bmag_Z^2 / (8π * density_Z * E₀_proton)


        # The Monte Carlo length is rg0 = γZ * βZ * E₀_proton / (e * BmagZ).
        # The plasma skin depth is λ_SD = γZ * E₀_proton / (4π * e^2 * denZ),
        # where denZ refers to the upstream number density of electrons
        # With the definition
        #     σ = 2εB0/γZ = BmagZ^2 / (4π * γZ * denZ * E₀_proton),
        # where denZ here refers to the number density of *protons*, one can
        # show that in the shock frame (where grid exists),
        #     λ_SD = βZ/√(σ*density_p/density_e) * rg0.
        density_Z_electron = 0.0
        if sc_electron
            density_Z_electron = denZ_ion[n_ions]
        else
            density_Z_electron = dot(denZ_ion, zz_ion)
        end
        σ = 2εB0 / γ_Z
        rg2sd = β_Z / √(σ*density_Z/density_Z_electron)


        # Also need the final value of ε_B downstream, in case our DwS
        # region is long enough that the magnetic field can decay to this value
        # Note that the R-H relations can be rearranged to read
        #     energy_density(x)  =  F_en0/u(x) - F_px0
        # assuming flux conservation everywhere.
        energy_density = (flux_energy_UpS + γ_Z*u_Z*density_Z*E₀_proton) / ux_sk_grid[n_grid]  -  flux_px_UpS
        εB2 = (bmag_Z*comp_fac)^2 / (8π * energy_density)
        # Use this value to compute the distance downstream at which the field
        # will have decayed to it
        # Per the Blandford-McKee solution, energy ∝ 1/χ ∝ 1/distance DwS.
        # Since we do not actually modify our pressures and densities according
        # to the BM solution, instead modify ε_B
        end_decay_rg = (5e-3 / εB2) / rg2sd


        # Now we can calculate εB_grid...
        # Per the Blandford-McKee solution, energy ∝ 1/χ ∝ 1/distance DwS.
        # Since we do not actually modify our pressures and densities according
        # to the BM solution, instead modify ε_B
        for i in 1:n_grid
            if x_grid_rg[i]*rg2sd < -50
                εB_grid[i] = max( 10.4e-4 / abs(x_grid_rg[i]*rg2sd)^0.6, εB0 )
            elseif x_grid_rg[i]*rg2sd < 50
                εB_grid[i] = 1e-4
            elseif x_grid_rg[i] < end_decay_rg
                εB_grid[i] = 5e-3 / (x_grid_rg[i]*rg2sd)
            else
                εB_grid[i] = εB2
            end
        end


        # ...and use it to calculate btot_grid
        for i in 1:n_grid
            energy_density = (flux_energy_UpS + γ_Z*u_Z*density_Z*E₀_proton) / ux_sk_grid[i]  -  flux_px_UpS
            btot_grid[i] = √(8π * εB_grid[i] * energy_density)
        end

    else
        # No custom ε_B requested
        fill!(εB_grid, 1e-99)

    end
    #-------------------------------------------------------------------------
    # Custom εB_grid defined if needed

    bmag_2 = btot_grid[n_grid]

    return (ux_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
            β_ef_grid, γ_ef_grid, btot_grid, θ_grid, εB_grid, bmag_2)
end

"""
Initializes the particle populations that will propagate through the shock structure.
Handles fast push and associated flux-tracking & changes to the population

### Arguments

TODO

- inp_distr
- i_ion
- do_fast_push
- aa

### Returns

TODO

- n_pts_use
- i_grid_in
- xwt_in
- pt_pf_in
- pb_pf_in
- x_PT_cm_in
- pxx_flux
- pxz_flux
- energy_flux
"""
function init_pop(
        do_fast_push, inp_distr, i_ion, aa,
        # from controls module
        tZ_ion, energy_inj, inj_wt, n_pts_inj, denZ_ion, x_grid_start, rg0, η_mfp,
        x_fast_stop_rg, β_Z, γ_Z, u_Z, n_ions, aa_ion, zz_ion, sc_electron, tZ_electron, oblique,
        # from grid_vars module
        n_grid, x_grid_rg, ux_sk_grid, γ_sf_grid,
        # from iteration_vars module
        ptot_inj, wt_inj, n_pts_MB,
    )

    # If not using fast push, this procedure is quite quick.
    # Fill the arrays and return to the main loop
    if ! do_fast_push
        if inp_distr == 1
            T_or_E = tZ_ion[i_ion]
        elseif inp_distr == 2
            T_or_E = energy_inj
        else
            error("not set to handle inp_distr > 2")
        end

        ptot_inj[:,i_ion], wt_inj[:,i_ion], n_pts_MB[i_ion] = set_in_dist(
            inj_wt, n_pts_inj, inp_distr, T_or_E, aa, denZ_ion[i_ion])

        n_pts_use = n_pts_MB[i_ion]
        xwt_in = wt_inj[1:n_pts_use, i_ion]
        pt_pf_in = ptot_inj[1:n_pts_use, i_ion]
        pb_pf_in = pt_pf_in[1:n_pts_use] * 2*(rand(n_pts_use) .- 0.5)
        x_PT_cm_in = fill(x_grid_start - 10*rg0*η_mfp, n_pts_use)
        pxx_flux = zeros(na_grid)
        pxz_flux = zeros(na_grid)
        energy_flux = zeros(na_grid)

        i_grid_in = zeros(Int, n_pts_use)
        return n_pts_use, i_grid_in, xwt_in, pt_pf_in, pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux
    end


    ##
    ## Everything from here on assumes fast push
    ##


    # Fast push won't work with any distribution besides thermal (yet?), so
    # stop if this occurs
    inp_distr > 1 && error("fast push will only work with thermal input distr.")


    # Compute density, pressure, and temperature at termination of fast push;
    # check to ensure that #assumecold still holds
    # Start counting at zero so that same loop applies even if fast push *isn't* enabled
    i_stop = findfirst(>(x_fast_stop_rg), x_grid_rg) - 1
    @debug "Found i_stop" i_stop

    relativistic = ( β_Z ≥ β_rel_fl )
    density_ratio = u_Z / ux_sk_grid[i_stop]
    if relativistic
        density_ratio *= γ_Z / γ_sf_grid[i_stop]
    end

    # Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
    # to calculate the far UpS internal energy   #assumecold
    γ_sph = 5//3

    temp_ratio = density_ratio^γ_sph / density_ratio

    if (kB_cgs * tZ_ion[i_ion] * temp_ratio) > (4 * aa*mp_cgs*c_cgs^2 * energy_rel_pt)
        error("Fast push cannot work b/c highest energy thermal particles become mildly relativistic. ",
              "Move fast push location UpS or disable entirely.")
    end
    pxx_flux = zeros(na_grid)
    pxz_flux = zeros(na_grid)
    energy_flux = zeros(na_grid)


    # Only run through the flux updates for the first particle species (i.e.
    # protons, not that it matters here); skip thereafter
    if i_ion == 1
        flux_update!(
            pxx_flux, pxz_flux, energy_flux,
            aa_ion, zz_ion, denZ_ion, tZ_ion, sc_electron, oblique, relativistic,
            i_stop,
            γ_Z, u_Z, γ_sf_grid, ux_sk_grid,
        )
    end  # check on i_ion


    # With fast push fluxes taken care of, create the particle distribution
    # that will be injected at x_fast_stop
    T_or_E = tZ_ion[i_ion] * temp_ratio

    ptot_inj[:,i_ion], wt_inj[:,i_ion], n_pts_MB[i_ion] = set_in_dist(
        inj_wt, n_pts_inj, inp_distr, T_or_E, aa, denZ_ion[i_ion])
    n_pts_use = n_pts_MB[i_ion]

    xwt_in  = wt_inj[1:n_pts_use, i_ion]
    pt_pf_in = ptot_inj[1:n_pts_use, i_ion]
    x_PT_cm_in = fill(x_fast_stop_rg * rg0, n_pts_use)
    i_grid_in  = fill(i_stop, n_pts_use)

    pb_pf_in = zeros(n_pts_use)
    for i_prt in 1:n_pts_use

        # Per Vladimirov+ (2009) [PhD], particle velocities should not be isotropic
        # in plasma frame. Must be weighted according to shock-frame pitch angle.
        #------------------------------------------------------------------------
        rand = Random.rand()

        if relativistic
            γ_pt_pf  = hypot(1, pt_pf_in[i_prt] / (aa*mp_cgs*c_cgs))
            vt_pf    = pt_pf_in[i_prt] / (γ_pt_pf * aa*mp_cgs)
            vmin_sf = (ux_sk_grid[i_stop] - vt_pf) / ( 1  -  ux_sk_grid[i_stop]*vt_pf/c_cgs^2 )
            vmax_sf = (ux_sk_grid[i_stop] + vt_pf) / ( 1  +  ux_sk_grid[i_stop]*vt_pf/c_cgs^2 )

            vxsf    = √( (vmax_sf^2 - vmin_sf^2)*rand  +  vmin_sf^2 )
            vx_pf    = (vxsf - ux_sk_grid[i_stop]) / ( 1  -  vxsf*ux_sk_grid[i_stop]/c_cgs^2 )
        else
            γ_pt_pf  = 1.0
            vt_pf    = pt_pf_in[i_prt] / (aa*mp_cgs)
            vmin_sf = ux_sk_grid[i_stop] - vt_pf
            vmax_sf = ux_sk_grid[i_stop] + vt_pf

            vxsf    = √( (vmax_sf^2 - vmin_sf^2)*rand  +  vmin_sf^2 )
            vx_pf    = vxsf - ux_sk_grid[i_stop]
        end

        pb_pf_in[i_prt] = γ_pt_pf * aa*mp_cgs * vx_pf
        #------------------------------------------------------------------------
        # Velocity-weighted pitch angles finished
    end

    return n_pts_use, i_grid_in, xwt_in, pt_pf_in, pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux
end


function flux_update!(
        pxx_flux, pxz_flux, energy_flux,
        aa_ion, zz_ion, denZ_ion, tZ_ion, sc_electron, oblique, relativistic,
        i_stop,
        γ_Z, u_Z, γ_sf_grid, ux_sk_grid
    )
    # Update the flux arrays as if the particles had actually crossed them
    #-----------------------------------------------------------------------
    # Obtain UpS pressure and density
    pressure_Z = dot(denZ_ion, tZ_ion) * kB_cgs
    ρ_Z        = dot(denZ_ion, aa_ion) * mp_cgs
    mask = (aa_ion .≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])

    # If electrons were not a separate species, add them in here
    if ! sc_electron
        pressure_Z += density_electron * kB_cgs*tZ_electron
        ρ_Z        += density_electron * me_cgs
    end

    # Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
    # to calculate the far UpS pressure and internal energy
    # #assumecold
    γ_sph = 5//3


    # Calculate fluxes and update arrays; note that if fast push isn't
    # enabled i_stop = 0, and this loop never executes
    for i in 1:i_stop

        density_ratio = (γ_Z * u_Z) / (γ_sf_grid[i] * ux_sk_grid[i])
        ρ_curr        = ρ_Z * density_ratio

        # Note assumption that γ_sph doesn't change from zone to zone: #assumecold
        pressure_curr = pressure_Z * density_ratio^γ_sph

        β_curr   = ux_sk_grid[i] / c_cgs
        γ_β_curr = γ_sf_grid[i] * ux_sk_grid[i] / c_cgs

        # Determine fluxes while handling different possible orientations and shock speeds.
        # For non-rel fluxes, expand out to β^2 to allow for better matching
        # with relativistic versions.
        # WARNING: these fluxes do not include contributions from a strong
        # magnetic field. This is incorporated during the smoothing process.
        #----------------------------------------------------------------------
        if !oblique
            flux_pz = 0.0
            if !relativistic
                flux_px = ρ_curr * ux_sk_grid[i]^2 * ( 1 + β_curr^2 ) + pressure_curr * ( 1 + γ_sph/(γ_sph-1) * β_curr^2 )
                flux_energy = (ρ_curr/2 * ux_sk_grid[i]^3 * ( 1 + 1.25*β_curr^2 )
                               + pressure_curr * ux_sk_grid[i] * γ_sph/(γ_sph-1) * ( 1 + β_curr^2 ))
            else

                flux_px = pressure_curr + γ_β_curr^2 * ( ρ_curr * c_cgs^2 + γ_sph/(γ_sph-1)*pressure_curr)
                flux_energy = γ_β_curr^2 *c_cgs / (ux_sk_grid[i]/c_cgs) * (ρ_curr*c_cgs^2 + γ_sph/(γ_sph-1)*pressure_curr )

                # Subtract mass-energy flux from flux_energy to bring it in line with non-rel calculations
                flux_energy -= γ_β_curr*c_cgs * ρ_curr * c_cgs^2
            end
        else
            error("Fast push cannot handle oblique shocks yet")
        end
        #--------------------------------------------------------------------
        # Fluxes calculated

        pxx_flux[i] = flux_px
        pxz_flux[i] = flux_pz
        energy_flux[i]  = flux_energy

    end  # loop over grid location
    #------------------------------------------------------------------------
    # Arrays updated through i_fast_stop
end

"""
Sets the injected particle distributions for all species. Initially, particles are placed in
a Maxwell-Boltzmann (thermal) distribution based on supplied temperature and particle mass.
This is corrected at the end if a δ-function distribution was requested (not worried about
the wasted computation because this subroutine runs only rarely).

### Arguments
FIXME
- inj_wt: whether each particle or each bin has equal weights (T for equal weight particles,
  F for equal weight bins)
- n_pts_inj: target # of particles for distribution
- inp_distr: thermal, δ-function, or some other distribution
- T_or_E: if using thermal distribution, this is temperature[K]; if δ function, it's injection energy[keV]
- aa: mass number for this particle species
- density_Z: far UpS number density for this species

### Returns
- ptot_out: array holding plasma frame total momenta for all particles in the distribution
- wt_out: array holding particle weights
- n_pts_use: number of particles in the distribution; will almost surely be different from
  n_pts_inj if using thermal dist

CHECKTHIS: that output distribution matches M-B, just to make sure I haven't made a typo
"""
function set_in_dist(inj_wt::Bool, n_pts_inj, inp_distr, T_or_E, aa, density_Z)

    # Error prevention
    0 < inp_distr < 3 || throw(DomainError(inp_distr, "Code can only do inp_distr = 1 or 2."))

    ptot_out = zeros(na_particles)
    wt_out = zeros(na_particles)
    # Administrative constants, e.g. total number of particles to distribute
    #--------------------------------------------------------------------------
    if inj_wt
        # Can't say anything about total number of particles in this case, b/c
        # haven't split them into M-B distribution yet
    else
        n_per_bin = n_pts_inj ÷ num_therm_bins   # Integer math loses excess particles, but that's intended behavior
        n_pts_tot = n_per_bin * num_therm_bins

        if n_per_bin < 5
            throw(ArgumentError("too few particles per bin ($n_per_bin; need at least 5). Increase n_pts_inj."))
        end

    end
    #------------------------------------------------------------------------
    # End administrative section


    # Set a mess of constants to be referred to routinely
    #------------------------------------------------------------------------
    particle_mass = aa*mp_cgs
    rm_energy = particle_mass * c_cgs^2

    kT      = kB_cgs * T_or_E  # Working under assumption of thermal dist now
    kT_o_rm = kT / rm_energy
    # Minimum, maximum extent of Maxwell-Boltzmann distribution
    kT_min  = 2e-3 * kT
    kT_max  = 10   * kT

    kT_rel_div = energy_rel_pt

    # Find min and max momenta of M-B curve
    if kT_o_rm < kT_rel_div
        # In non-rel case, kinetic energy ≈ thermal energy
        p_min = √( 2particle_mass * kT_min )
        p_max = √( 2particle_mass * kT_max )
    else
        # Once particles are relatvistic, rest-mass energy becomes important:
        #     E² = p²c² + m²c⁴  ≈  (kT + mc²)²
        p_min = √( (kT_min + rm_energy)^2 - rm_energy^2 )  /  c_cgs
        p_max = √( (kT_max + rm_energy)^2 - rm_energy^2 )  /  c_cgs
    end

    Δp = (p_max - p_min) / num_therm_bins
    #------------------------------------------------------------------------
    # End of constants section


    # Generate the Maxwell-Boltzmann distribution. The actual calculation of
    # f1 and f2 does not need modification for relativistic momenta.
    # Such a modification would only affect the normalization of the curve,
    # not the dependence on momentum for a particular value of kT. Since we
    # only care about the relative fractional area of each bin, overall
    # normalization doesn't matter.
    #------------------------------------------------------------------------
    # Find total area under M-B curve
    area_tot = 0.0
    for i in 1:num_therm_bins
        p1 = p_min + (i-1)*Δp
        p2 = p_min +  i   *Δp

        if kT_o_rm < kT_rel_div
            energy_o_kT1 = p1^2 / (2particle_mass * kT)
            energy_o_kT2 = p2^2 / (2particle_mass * kT)
        else
            energy_o_kT1 = hypot(p1*c_cgs, rm_energy) / kT
            energy_o_kT2 = hypot(p2*c_cgs, rm_energy) / kT
        end

        # Start working in log space because of potentially huge exponents
        f1 = exp(2log(p1) - energy_o_kT1)
        f2 = exp(2log(p2) - energy_o_kT2)

        # Integrate using the trapezoid rule
        area_tot += Δp * 0.5 * (f1 + f2)
    end


    # Fill the bins with particles. Needs to be done in two separate loops,
    # despite the similarities, because particle weights are handled
    # differently based on value of inj_wt
    if inj_wt # First, particles have equal weight

        area_per_pt = area_tot / n_pts_inj # Area each particle gets if inj_wt = T
        n_pts_tot   = 1                    # Total number of particles in M-B distribution

        for i in 1:num_therm_bins
            p1 = p_min + (i-1)*Δp
            p2 = p_min +  i   *Δp

            if kT_o_rm < kT_rel_div
                energy_o_kT1 = p1^2 / (2particle_mass * kT)
                energy_o_kT2 = p2^2 / (2particle_mass * kT)
            else
                energy_o_kT1 = hypot(p1*c_cgs, rm_energy) / kT
                energy_o_kT2 = hypot(p2*c_cgs, rm_energy) / kT
            end

            # Work in log space because of potentially huge exponents
            f1 = exp(2log(p1) - energy_o_kT1)
            f2 = exp(2log(p2) - energy_o_kT2)

            # Calculate bin area using trapezoid rule
            bin_area = Δp * (f1 + f2)/2

            area_frac = bin_area / area_per_pt

            # Rounded particle count for this bin
            n_pts_this_bin = round(Int, area_frac)

            # Geometric center of bin; particles in this bin will receive this momentum
            # particles in this bin will have momentum = the geometric center of the bin
            ptot_out[n_pts_tot:n_pts_tot+n_pts_this_bin] .= √(p1*p2)
            n_pts_tot += n_pts_this_bin
        end

        # If each particle has equal weight, then the total weight of the
        # distribution should be proportional to the density of the species
        wt_out[1:n_pts_tot] .= density_Z / n_pts_tot

    else # Bins have equal weight
        for i in 1:num_therm_bins
            p1 = p_min + (i-1)*Δp
            p2 = p_min +  i   *Δp

            if kT_o_rm < kT_rel_div
                energy_o_kT1 = p1^2 / (2particle_mass * kT)
                energy_o_kT2 = p2^2 / (2particle_mass * kT)
            else
                energy_o_kT1 = hypot(p1*c_cgs, rm_energy) / kT
                energy_o_kT2 = hypot(p2*c_cgs, rm_energy) / kT
            end

            # Work in log space because of potentially huge exponents
            f1 = exp(2log(p1) - energy_o_kT1)
            f2 = exp(2log(p2) - energy_o_kT2)

            # Calculate bin area using trapezoid rule
            bin_area = Δp * (f1 + f2)/2

            area_frac = bin_area / area_tot

            # particles in this bin will have momentum = the geometric center of the bin
            jstart = (i-1)*n_per_bin + 1
            jend   = jstart + (n_per_bin-1)
            ptot_out[jstart:jend] .= √(p1*p2)
            wt_out[jstart:jend] .= area_frac / n_per_bin * density_Z

        end  # loop over i

    end  # test of inj_wt

    n_pts_use = n_pts_tot
    #------------------------------------------------------------------------
    # End of Maxwell-Boltzmann section


    # If the particles are to use a δ-function distribution, set it here
    #--------------------------------------------------------------------------
    if inp_distr == 2
        n_pts_use = n_pts_inj

        rm_energy      = aa * mp_cgs * c_cgs^2
        energy_inj_cgs = T_or_E * keV2erg

        if energy_inj_cgs/rm_energy < energy_rel_pt
            p1 = √( 2 * aa * mp_cgs * energy_inj_cgs )
        else
            p1 = √( energy_inj_cgs^2 - rm_energy^2 )  /  c_cgs
        end

        ptot_out[1:n_pts_inj] .= p1
        wt_out[1:n_pts_inj] .= density_Z / n_pts_tot

    end  # check of inp_distr
    #--------------------------------------------------------------------------
    # δ-function handled


    return (ptot_out, wt_out, n_pts_use)
end
end # module
