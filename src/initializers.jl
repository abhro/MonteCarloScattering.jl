module initializers
export calc_DwS, calc_rRH, set_psd_mom_bins, init_pop, upstream_fluxes, upstream_machs, setup_profile, setup_grid, set_psd_angle_bins
using Random: rand
using OffsetArrays: OffsetVector, Origin
using LinearAlgebra: dot
using StaticArrays: SVector
using Roots: find_zero, Newton
using Unitful, UnitfulAstro
using Unitful: g, cm, s, dyn, erg, keV
using Unitful: mp, c, k as kB
using Distributions: Uniform
using ..constants: E₀_proton
using ..parameters
import ..density, ..temperature, ..mass, ..number_density
using ..CGSTypes

"""
    calc_DwS(...)

Uses the Rankine-Hugoniot jump conditions to calculate the downstream conditions for a test
particle shock. Big difference between this subroutine and calc_rRH is that we already know
what the DwS speed is, courtesy of r_comp in the input.

Shock must be parallel (not oblique).

### Arguments

- `B₀`: far UpS magnetic field strength[Gauss]
- `r_comp`: compression ratio of shock
- `β₀`: UpS bulk fluid speed, over c

### Returns

- `β`: (total) bulk fluid speed DwS
- `γ`: Lorentz factor associated with β₂
- `B`: DwS magnetic field strength[Gauss]
- `θ_B`: angle[deg] between DwS magnetic field and shock normal
- `θᵤ`: angle[deg] between DwS fluid velocity and shock normal
"""
function calc_DwS(B₀, r_comp, β₀)
    β   = β₀ / r_comp
    γ   = 1 / √(1 - β^2)
    B   = B₀
    θ_B = 0.0
    θᵤ  = 0.0
    return β, γ, B, θ_B, θᵤ
end

"""
    calc_rRH(...)

Uses the Rankine-Hugoniot jump conditions to calculate the compression ratio for a shock
assuming test-particle conditions. In other words, (1) sharp shock, (2) negligible/no DSA,
and (3) no escaping flux. Additionally assumes that the inflowing plasma has
non-relativistic thermal speeds to make UpS adiabatic index exactly 5/3.

### Arguments
- `β₀`, `γ₀`: shock speed and Lorentz factor
- `n_ions`: number of different ion species
- `m_ion`: array of particle species' mass
- `n₀_ion`: array of far UpS densities for particle species
- `T₀_ion`: array of particle species' far UpS temperatures

Shock must be parallel (not oblique)

Unused:
- `oblique`: controls whether to use parallel or oblique formulations of R-H relations

### Returns
- `r_RH`: Rankine-Hugoniot compression ratio
- `Γ₂_RH`: ratio of specific heats (adiabatic index) for DwS region, assuming r_comp = r_RH
"""
function calc_rRH(u₀, β₀, γ₀, species)

    #--------------------------------------------------------------------------
    #  Two possibilities for R-H relations: nonrelativistic/relativistic
    #  Determine which of the four to use. Cutoff for nonrelativistic/relativistic is set in module 'parameters'
    #--------------------------------------------------------------------------
    relativistic = (β₀ < β_rel_fl)

    # Calculate thermal pressure of far upstream gas
    P₀ = dot(density.(species), temperature.(species)) * Unitful.k
    ρ₀ = dot(density.(species), mass.(species))

    #--------------------------------------------------------------------------
    #  Possibility 1: Nonrelativistic, parallel
    #  Solution comes from Ellison (1985) [1985JGR....90...29E].
    #  Uses far UpS Mach number to calculate r_RH.
    #--------------------------------------------------------------------------
    if !relativistic
        r_RH, Γ₂_RH = calc_rRH_nonrelativistic(P₀, ρ₀, β₀)

    #--------------------------------------------------------------------------
    #  Possibility 2: Relativistic, parallel
    #  Solution comes from Ellison & Reynolds (1990) [1991ApJ...378..214E].
    #  Uses relativistic Rankine-Hugoniot relations. See that paper for
    #  details of equations and associated quantities. Briefly,
    #      R-H1:  γ₁  n₁ b₁        =  γ₂  n₂ b₂
    #      R-H2:  γ₁² w₁ b₁² + P₁  =  γ₂² w₂ b₂² + P₂
    #      R-H3:  γ₁² w₁ b₁        =  γ₂² w₂ b₂
    #  where
    #      w    = E_rm + E_ke + P   ← enthalpy as total energy density + pressure
    #      E_rm = n m c²            ← rest mass energy density
    #      E_ke = n m c² (γ - 1)    ← kinetic energy density, with γ = √(1 + (p/mc)²)
    #      P    = ⅓ n p v           ← pressure
    #
    #  Assumes that downstream particle distributions are δ-functions.
    #  Solves for p2 using Newton's method, then works backwards to r_RH.
    #-------------------------------------------------------------------------
    else
        r_RH, Γ₂_RH = calc_rRH_relativistic(species, ρ₀, P₀, β₀, n₀_ion)
    end

    return r_RH, Γ₂_RH
end

"""
    calc_rRH_nonrelativistic(P₀, ρ₀, β₀)

TODO
"""
function calc_rRH_nonrelativistic(P₀, ρ₀, β₀)

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal
    # gas, to calculate the far UpS sound speed and Mach number   #assumecold
    Γ_sph = 5//3
    cₛ    = √(Γ_sph * P₀ / ρ₀)
    M_Z   = β₀ * Unitful.c / cₛ |> NoUnits

    # Finally, use Equation (11) from Ellison (1985) to calculate r_RH.
    # Note that q = 0 here b/c we assume no escaping flux. This simplifies
    # the denominator quite a bit from the equation.
    r_RH = 8 / (2 + 6/M_Z^2)

    # In non-relativistic case, downstream adiabatic index is pegged to 5/3
    Γ₂_RH = 5//3

    return r_RH, Γ₂_RH
end
"""
    calc_rRH_relativistic(species, ρ₀, P₀, β₀, n₀_ion)

TODO
"""
function calc_rRH_relativistic(species, ρ₀, P₀, β₀, n₀_ion)

    # FIXME the comment refers to old version of variables
    # Calculate two quantities to be used during loop to find r_RH: the
    # rest mass-energy of each species, and the (number) density relative to protons
    relative_ion_energy = (
        dot(mass.(species),     # rest energy of each species (to be multiplied by c²)
            density.(species))  # density relative to protons (to be divided by n₀_proton)
        * c^2                   # turn mass into rest energy of all species
        / first(density.(species))) # turn densities into density relative to protons

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal gas,
    # to calculate the far UpS enthalpy      #assumecold
    Γ_sph = 5//3
    w₀ = ρ₀ * c^2 + Γ_sph/(Γ_sph-1) * P₀

    # Calculate the far UpS momentum flux
    UpS_mom_flux = γ₀^2 * w₀ * β₀^2 + P₀
    UpS_num_flux = γ₀ * n₀_ion[1] * β₀ # Protons only here; not strictly correct but appropriate for later use

    function F(p)
        γβ = p / (mp*c) # since p is relativistic, p/mc = γmv/mc = γ⋅β
        γ = √(1 + γβ^2)
        P = relative_ion_energy/3 * γβ^2/γ # pressure
        w = relative_ion_energy*(γ + γβ^2/3γ)
        return UpS_num_flux/γ * (w*γ^2 + P) - UpS_mom_flux
    end

    # Now use Newton's method to determine the downstream momentum that satisfies the
    # R-H relations. Assumptions: (1) momentum distribution functions are δ-functions
    # rather than thermal, and (2) any non-proton species have p ∝ m
    p₂_found = find_zero(F, 0, Newton())

    # Calculate the compression ratio β₀/β₂ associated with p₂_found
    γβ = p₂_found / (mp*c)  # Protons only here because of how the math in w_fac works out

    # Pressure, internal energy, and enthalpy with proton density factored out
    P_fac = relative_ion_energy/3 * γβ^2 / √(1 + γβ^2)
    e_fac = relative_ion_energy * (√(1 + γβ^2) - 1)
    w_fac = relative_ion_energy * (√(1 + γβ^2) + 1/3 * γβ^2 / √(1 + γβ^2))

    # Calculate adiabatic index downstream
    Γ₂_RH = 1 + P_fac/e_fac

    # Finally, get downstream speed and compression ratio
    β₂ = √(1 - (n₀_ion[1] * w_fac/γ₀ * w₀)^2)  # Using only proton density is correct


    r_RH = β₀/β₂
    #------------------------------------------------------------------------
    # r_RH found using Newton's method

    return r_RH, Γ₂_RH
end

"""
    set_psd_mom_bins(...)

Sets the BOUNDARIES of the bins of the phase space distribution. The bins are numbered from
0 to num_psd_mom_bins, each boundary denotes the lower edge of that # bin; the indices thus
run from 0 to num_psd_mom_bins + 1.

Logarithmic spacing over all decades from Emin_keV to Emax_keV is used.
Values less than the minimum are equivalent to 0.0.

### Arguments
- `psd_mom_min`: minimum momentum to use in PSD
- `psd_mom_max`: maximum momentum to use in PSD
- `psd_bins_per_dec_mom`: # of bins per decade to use in logarithmically-spaced regions of PSD

### Returns
- `num_psd_mom_bins`: total number of bins along given dimension, not counting bin 0
- `psd_mom_bounds`: boundaries between bins, and upper edge of final bin
"""
function set_psd_mom_bins(psd_mom_min, psd_mom_max, psd_bins_per_dec_mom)
    num_psd_mom_bins = trunc(Int, log10(psd_mom_max / psd_mom_min) *  psd_bins_per_dec_mom)
    num_psd_mom_bins += 2   # Add two extra bins just to be safe

    # Fill in the array psd_mom_bounds, remembering that the array holds LOWER boundaries
    # of that bin
    psd_mom_bounds = Origin(0)(Float64[-99.0])
    append!(psd_mom_bounds,
            range(start  = log10(psd_mom_min/(mp*c)),
                  step   = 1/psd_bins_per_dec_mom,
                  length = num_psd_mom_bins+1))

    length(psd_mom_bounds) == num_psd_mom_bins+2 || error() # sanity check

    return num_psd_mom_bins, psd_mom_bounds
end

"""
    set_psd_angle_bins(...)

Sets the BOUNDARIES of the bins of the phase space distribution. The bins are numbered from
0 to num_psd_θ_bins, each boundary denotes the lower edge of that # bin; the indices thus
run from 0 to num_psd_θ_bins + 1.

!!! warning
    the number stored in psd_θ_bounds increases at first (increasing θ),
    then decreases because increasing the angle means decreasing the cosine.

Linear spacing of cosine values for angles between psd_θ_fine and π. Below psd_θ_fine
spacing is logarithmic in θ for some number of decades down to psd_θ_min. Values less than
the minimum are equivalent to 0.0.

### Arguments
- `psd_bins_per_dec_θ`: number of bins per decade to use in logarithmically-spaced regions of PSD
- `psd_lin_cos_bins`: number of bins to divide range [-1,cos(psd_θ_fine)] into
- `psd_cos_fine`: cutoff between lin/cos and log/θ spacing of PSD bins
- `psd_θ_min`: minimum angle for PSD

### Returns
- `num_psd_θ_bins`: total number of bins along given dimension, not counting bin 0
- `Δcos`: size of each linear cosine bin
- `psd_θ_bounds`: boundaries between bins, and upper edge of final bin
"""
function set_psd_angle_bins(psd_bins_per_dec_θ, psd_lin_cos_bins, psd_cos_fine, psd_θ_min)
    psd_θ_fine = acos(psd_cos_fine)
    ten_root_θ = exp10(1 / psd_bins_per_dec_θ)

    psd_log_θ_bins = trunc(Int, log10(psd_θ_fine/psd_θ_min) * psd_bins_per_dec_θ)
    #num_psd_θ_bins = psd_log_θ_bins + psd_lin_cos_bins

    # Fill the logarithmic part of psd_θ_bounds using the angle (in radians), NOT its logarithm
    psd_θ_bounds = Origin(0)(Float64[1e-99])
    append!(psd_θ_bounds, psd_θ_min * ten_root_θ.^range(0, length=psd_log_θ_bins))

    # Now fill in the linear part of psd_θ_bounds.
    # Note that the lower boundary of the first cell is psd_cos_fine
    Δcos = (psd_cos_fine + 1) / psd_lin_cos_bins
    append!(psd_θ_bounds, range(start=psd_cos_fine, step=-Δcos, length=psd_lin_cos_bins+1))

    #return num_psd_θ_bins, Δcos, psd_θ_bounds
    return Δcos, psd_θ_bounds
end

"""
    set_photon_shells(...)

If photon calculation is desired, photons will be collected into UpS and DwS shells for
easier viewing. This subroutine sets the endpoints of the shells, as well as their midpoints.

Because the particle spectrum changes most rapidly near the shock, zones should be small
near the shock and get larger as you move further away.

First, divide the domain between the shock and the FEB (either UpS or DwS into n sections,
where n is the respective number of shells to use. HOWEVER, do this by exponent, ranging
from -1 to log10(|x_FEB|).

Keep track of the boundaries of each shells, since we will need those for calculating the
total number of particles emitting when we get to that point in photon production. The
midpoints are less useful, since photons are calculated on a zone-by-zone basis rather than
just at select points in the shock profile. Keep them anyway, since the memory overhead is low.
"""
function set_photon_shells(
        num_UpS_shells, num_DwS_shells,
        use_prp, feb_UpS, feb_DwS, rg₀, x_grid_stop_rg,
    )

    total_shells = num_UpS_shells+num_DwS_shells
    x_shell_midpoints = zeros(total_shells)
    x_shell_endpoints = zeros(total_shells+1) # fencepost problem

    # Handle UpS shells first
    set_UpS_photon_shells!(x_shell_midpoints, x_shell_endpoints, num_UpS_shells, feb_UpS, rg₀)

    # And repeat the process for the downstream shells.
    set_DwS_photon_shells!(x_shell_midpoints, x_shell_endpoints, num_UpS_shells, num_DwS_shells,
                           use_prp, feb_DwS, rg₀, x_grid_stop_rg)

    # Convert from units of rg₀ to cm
    @. x_shell_endpoints *= rg₀
    return (x_shell_midpoints, x_shell_endpoints)
end

"""
    set_UpS_photon_shells!(...)

TODO
"""
function set_UpS_photon_shells!(
        x_shell_midpoints, x_shell_endpoints,
        num_UpS_shells, feb_UpS, rg₀,
    )
    x_section_width = (log10(abs(feb_UpS/rg₀))+1) / num_UpS_shells
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
        # Also, add in the factor of -1 here, since upstream coordinates should be negative
        # in the MC code.
        x_shell_midpoints[num_UpS_shells+1 - i]     = -x_region_mid
        x_shell_endpoints[num_UpS_shells+1 - i]     = -x_region_end
        x_shell_endpoints[num_UpS_shells+1 - i + 1] = -x_region_start
    end
end
"""
    set_DwS_photon_shells(...)

TODO
"""
function set_DwS_photon_shells(
        x_shell_midpoints, x_shell_endpoints,
        num_UpS_shells, num_DwS_shells, use_prp, feb_DwS, rg₀, x_grid_stop_rg,
    )
    # The DwS limit is set differently if using a PRP or a FEB
    limitDwS = use_prp ? x_grid_stop_rg : feb_DwS/rg₀
    x_section_width = (log10(limitDwS)+1) / num_DwS_shells

    for i in 1:num_DwS_shells
        # Calculate UpS and DwS endpoints of each region, as well as midpoint in log space
        if i == 1
            # Special case when i = 1, since our region starts at the shock.
            x_region_start = 0.0
            x_region_end   = exp10(-1 + x_section_width)
            x_region_mid   = exp10(-1 + x_section_width/2)
        else
            # In the general case, note that x_region_start should be the same as the previous
            # region's x_region_end. This can be checked with print statements at runtime.
            x_region_start = exp10(-1 + x_section_width * (i-1)  )
            x_region_end   = exp10(-1 + x_section_width *  i     )
            x_region_mid   = exp10(-1 + x_section_width * (i-1/2))
        end

        # Update the arrays with the information. Less index juggling here
        x_shell_endpoints[num_UpS_shells + i]     = x_region_start
        x_shell_midpoints[num_UpS_shells + i]     = x_region_mid
        x_shell_endpoints[num_UpS_shells + i + 1] = x_region_end
    end
end

# Many grid zones are set manually; zones can easily be added/removed, but
# make sure to change the number in the log-spaced regions UpS or DwS
const FIRST_ZONE = SVector(-9.0, -8.0, -7.0, -6.0, -5.0, -4.5, -4.0, -3.5, -3.0,
                           -2.5, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0,
                           -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2,
                           -0.15, -0.1,
                           -0.07, -0.05, -0.04, -0.03, -0.02, -0.015, -0.01,
                           -3e-3, -1e-3)
# Extremely fine spacing right around the shock
const EXTREMELY_FINE_SPACING = SVector(-1e-4, -1e-7, 0.0, 1e-7, 1e-4)

# Downstream from the shock, spacing doesn't need to be quite so fine
# because velocity gradients aren't as extreme, if they exist at all
const DOWNSTREAM_SPACING = SVector(1e-3, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 0.10,
                                   0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

"""
    setup_grid(...)

TODO
"""
function setup_grid(x_grid_start_rg, x_grid_stop_rg, use_prp, feb_DwS, rg₀)

    # Recall that rg₀ is the gyroradius of a proton with speed u₀ in magnetic field bmag₀.

    # Set the start and stop positions in units of rg₀
    x_grid_start = x_grid_start_rg * rg₀
    if !use_prp
        x_grid_stop_rg = feb_DwS / rg₀
        @info("DwS FEB set at x = $x_grid_stop_rg rg₀. Overwriting entered value for 'XGDDW'.")
    end
    x_grid_stop = x_grid_stop_rg * rg₀

    # Logarithmically-spaced grid zones run from x_grid_start_rg to -10rg₀. Set them here.
    n_log_UpS = 27
    Δlog = (log10(-x_grid_start_rg) - 1) / n_log_UpS-1

    # we build this chunk-by-chunk, which is inefficient
    x_grid_rg = Origin(0)(Float64[])

    push!(x_grid_rg, -1e30) # set left boundary of grid

    append!(x_grid_rg, -exp10.(range(start=log10(-x_grid_start_rg), step=-Δlog, length=n_log_UpS)))
    append!(x_grid_rg, FIRST_ZONE)
    append!(x_grid_rg, EXTREMELY_FINE_SPACING)
    append!(x_grid_rg, DOWNSTREAM_SPACING)

    # As seen above, the manually-set grid zones end at x = +1rg₀.
    # DwS from there, more log-spaced zones.
    n_log_DwS = 16
    x_end_man = x_grid_rg[end]
    Δlog      = (log10(x_grid_stop_rg) - log10(x_end_man)) / n_log_DwS

    append!(x_grid_rg, exp10.(range(start=log10(x_end_man), step=Δlog, length=n_log_DwS)))

    push!(x_grid_rg, 1e30)  # set right boundary of the grid

    return x_grid_rg, x_grid_start, x_grid_stop
end

"""
    upstream_fluxes(...)

Calculates the far upstream fluxes for the shock.

Two different cases considered:

1. Non-relativistic oblique shock. Uses equations of Ellison+ (1996) [1996ApJ...473.1029E]
2. Relativistic shock, any obliquity. Uses equations of Double+ (2004) [2004ApJ...600..485D]

Only oblique equations used because they reduce trivially to parallel cases when θ_B₀ = 0.

HOWEVER, assumes that z-component of far UpS velocity is 0 in all cases; in practice oblique
shocks would induce some z-velocity in the shock profile even though particles initially
arrive with no bulk z component. Also assumes isotropic initial pressure, so no off-diagonal
components in pressure tensor.

Shock must be parallel (not oblique). Set `oblique` to `false`.

### Arguments

FIXME

No inputs; pulls everything from module 'controls'

### Returns
- flux_px_UpS: far UpS momentum flux, x component
- flux_pz_UpS: far UpS momentum flux, z component
- flux_energy_UpS: far UpS energy flux
"""
function upstream_fluxes(n₀_ion, T₀_ion, m_ion, B₀, θ_B₀, u₀, β₀, γ₀)

    # UpS internal energy density and pressure, assuming isotropic particle distribution.
    # Note that this INCLUDES the mass-energy density, which is typically omitted in
    # nonrelativistic calculations
    P₀ = dot(n₀_ion, T₀_ion) * kB |> dyn/cm^2 # pressure
    ρ₀ = dot(n₀_ion, m_ion)       |> g/cm^3   # mass density
    @debug "calculated params" P₀ ρ₀

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal gas,
    # to calculate the far UpS internal energy             #assumecold
    Γ_sph = 5//3
    # internal energy density
    e₀ = ρ₀*c^2 + 1/(Γ_sph - 1) * P₀ |> erg/cm^3

    # Quantities related to the UpS magnetic field. Note that B_z is the
    # z-component of the magnetic field, not B₀
    B_x = B₀ * cosd(θ_B₀)
    B_z = B₀ * sind(θ_B₀)

    relativistic = β₀ ≥ β_rel_fl

    if relativistic
        flux_px_UpS, flux_pz_UpS = upstream_momentum_flux_relativistic(β₀, γ₀, e₀, P₀, B₀, B_x, B_z)
        flux_energy_UpS = upstream_energy_flux_relativistic(u₀, β₀, γ₀, e₀, ρ₀, P₀, B_z)
    else
        # Non-relativistic version. Note that it's missing the ρc² flux present in the
        # relativistic forms above. It is also expanded to second order in β₀ (only in the
        # hydro terms, for now) to allow for more precise matching with the relativistic version
        flux_px_UpS, flux_pz_UpS = upstream_momentum_flux_nonrelativistic(u₀, β₀, γ₀, e₀, ρ₀, P₀, B_x, B_z)
        flux_energy_UpS = upstream_energy_flux_nonrelativistic(u₀, β₀, γ₀, e₀, ρ₀, P₀)
    end

    return (flux_px_UpS, flux_pz_UpS, flux_energy_UpS)
end

"""
    upstream_momentum_flux_relativistic(...)

TODO
"""
function upstream_momentum_flux_relativistic(β₀, γ₀, e₀, P₀, B₀, B_x, B_z)

    # Momentum flux, x-component
    # Fluid part (Double+ Eq 23)
    F_pₓ_fl = (γ₀*β₀)^2 * (e₀ + P₀) + P₀ |> g/(cm*s^2)
    # EM part (Double+ Eq 25)
    F_pₓ_EM = γ₀^2 * ((β₀*B₀)^2 + B_z^2 - B_x^2) / 8π |> g/(cm*s^2)
    @debug "Found partial fluxes" F_pₓ_fl F_pₓ_EM
    flux_px_UpS = F_pₓ_fl + F_pₓ_EM                         # Total

    # Momentum flux, z-component (Fluid Part = 0, from Double+ Eq 24)
    # Total = EM part (Double+ Eq 26)
    flux_pz_UpS = -γ₀/4π * B_x * B_z

    return flux_px_UpS, flux_pz_UpS
end

"""
    upstream_momentum_flux_nonrelativistic(...)

TODO
"""
function upstream_momentum_flux_nonrelativistic(u₀, β₀, γ₀, e₀, ρ₀, P₀, B_x, B_z)
    flux_px_UpS = ρ₀ * u₀^2 * (1 + β₀^2) +
                  P₀ * (1 + Γ_sph/(Γ_sph-1)*β₀^2) +
                  B_z^2/8π
    flux_pz_UpS = - B_x * B_z / 4π
    return flux_px_UpS, flux_pz_UpS
end

"""
    upstream_energy_flux_nonrelativistic(...)

TODO
"""
function upstream_energy_flux_nonrelativistic(u₀, β₀, γ₀, e₀, ρ₀, P₀)
    return ρ₀ * u₀^3 * (1 + 1.25*β₀^2)/2 + P₀ * u₀ * Γ_sph/(Γ_sph-1) * (1+β₀^2) + u₀*B_z^2/4π
end

"""
    upstream_energy_flux_relativistic(...)

TODO
"""
function upstream_energy_flux_relativistic(u₀, β₀, γ₀, e₀, ρ₀, P₀, B_z)
    F_energy_fl = γ₀^2 * β₀ * (e₀ + P₀) # Fluid part (Double+ Eq 20)
    F_energy_EM = γ₀^2 * β₀ * B_z^2/4π  # EM part (Double+ Eq 21)
    # Total -- convert to cgs units!
    flux_energy_UpS = c * (F_energy_fl + F_energy_EM)

    # And subtract off the mass-energy flux to bring it in line with nonrelativistic
    # calculations and what the MC code actually tracks
    flux_energy_UpS -= γ₀ * u₀ * ρ₀*c^2

    return flux_energy_UpS
end

"""
    upstream_machs(β₀, species, B₀)

Calculates the sonic and Alfvén mach numbers for the shock.

- For speed of sound, uses Equation (13) of Fujimura & Kennel (1979) [1979A%26A....79..299F]
- For Alfvén wave speed, uses Equation (46) of Gedalin (1993) [1993PhRvE..47.4354G]

### Arguments

TODO

### Returns
- mach_sonic: sonic Mach number of UpS flow
- mach_alfven: Alfvénic Mach number of UpS flow
"""
function upstream_machs(β₀, species, B₀)

    # Assume cold UpS plasma, so that the adiabatic index is 5/3 identically   #assumecold
    Γ = 5//3

    P₀ = dot(number_density.(species), temperature.(species)) * kB  # pressure
    ρ₀ = dot(number_density.(species), mass.(species))              # mass density
    @debug "Found parameters" P₀ ρ₀

    relativistic = (β₀ ≥ β_rel_fl)

    return (mach_sonic_func(β₀*c, P₀, ρ₀, Γ, relativistic),
            mach_alfven_func(β₀*c, P₀, ρ₀, Γ, B₀, relativistic))
end

# TODO name these functions better. the `_func` suffix is because there are
# globals with the same name
"""
    mach_sonic_func(u, P, ρ, Γ, relativistic)

Return sonic mach number (ratio of speed to local speed of sound)
"""
function mach_sonic_func(u, P, ρ, Γ, relativistic)
    cₛ = relativistic ?
            sound_speed_relativistic(P, ρ, Γ) :
            sound_speed_nonrelativistic(P, ρ, Γ)
    return u / cₛ |> NoUnits
end

"""
    mach_alfven_func(u, P, ρ, Γ, B, relativistic)

Return Alfvénic mach number (ratio of speed to Alfvén wave group velocity)
"""
function mach_alfven_func(u, P, ρ, Γ, B, relativistic)
    v_A = relativistic ?
            alfven_speed_relativistic(P, ρ, Γ, B) :
            alfven_speed_nonrelativistic(ρ, B)
    return u / v_A |> NoUnits
end

function sound_speed_relativistic(P, ρ, Γ)
    # Find FK1979's R factor, the ratio of pressure to rest mass energy density
    R = P / (ρ * c^2) |> NoUnits # dimensionless auxiliary variable

    # Compute the speed of sound using Fujimura & Kennel
    #     cₛ²/c² = ΓR/(aR + 1)                              FK1979 Eq. 13
    a = Γ / (Γ - 1)     # defined near FK1979 Equation (6)
    @debug "Called with" P ρ Γ
    @debug "Calculated params" R a
    return c * √(Γ * R / (a*R + 1))
end

sound_speed_nonrelativistic(P, ρ, Γ) = √(Γ * P / ρ) # cₛ = √(K/ρ), where K = ΓP is the bulk modulus

function alfven_speed_relativistic(P, ρ, Γ, B)
    # And into Gedalin (1993)'s Equation (46); note assumption that
    # equation of state is    e = ρc² + P/(Γ-1)
    #     v_A² = (B²/4π) / (ε + p + B²/4π)                  Gedalin Eq. 46
    enthalpy = Γ/(Γ-1) * P + ρ * c^2
    # In principle 1 erg/cm³ == 1 G², but uconvert breaks when doing it
    v_A = c / √(1 + 4π * enthalpy / B^2)
    # more on equation of state
    #    e = ρc² + P/(Γ-1)
    # where
    # - e is the (total internal) energy density
    # - ρc² is the rest energy density
    # - P/(Γ-1) is the thermal component (internal kinetic energy)
    return v_A
end
alfven_speed_nonrelativistic(ρ, B) = B / √(4π * ρ) # Alfvén wave group velocity

"""
    setup_profile(...)

Sets the initial values of the shock profile

### Arguments

- `u₀`, `β₀`, `γ₀`
- `bmag₀`
- `θ_B₀`
- `r_comp`
- `bturb_comp_frac`
- `bfield_amp`
- `use_custom_εB,`
- `n_ions`
- `species`
- `flux_px_UpS`
- `flux_energy_UpS`
- `grid_axis`
- `x_grid_cm`
- `x_grid_rg`

### Returns
- uₓ_sk_grid: bulk fluid velocity along x axis (i.e., perpendicular to shock face) in shock frame
- uz_sk_grid: bulk fluid velocity along z axis (i.e., parallel to shock face) in shock frame
- utot_grid: total bulk fluid velocity in shock frame
- γ_sf_grid: bulk flow Lorentz factor in shock frame
- β_ef_grid: relative x-axis speed between plasma and explosion frames
- γ_ef_grid: Lorentz factor associated with β_ef_grid
- btot_grid: total magnetic field
- θ_grid: angle of magnetic field[radians] relative to shock normal (i.e., to x axis)
- εB_grid: user-defined function for fraction of energy density in magnetic field.
  Sets value of btot_grid
- bmag₂: field strength in DwS region. Initially set in calc_DwS, it may be reset here
  depending on values of bturb_comp_frac & bfield_amp
"""
function setup_profile(
        u₀, β₀, γ₀, bmag₀, θ_B₀,
        r_comp, bturb_comp_frac, bfield_amp, use_custom_εB,
        n_ions, species, flux_px_UpS, flux_energy_UpS,
        grid_axis, x_grid_cm, x_grid_rg,
    )

    uₓ_sk_grid = OffsetVector{typeof(u₀)}(undef, grid_axis)
    uz_sk_grid = zeros(typeof(u₀), grid_axis)
    γ_sf_grid  = OffsetVector{Float64}(undef, grid_axis)
    β_ef_grid  = OffsetVector{Float64}(undef, grid_axis)
    γ_ef_grid  = OffsetVector{Float64}(undef, grid_axis)
    btot_grid  = OffsetVector{typeof(bmag₀)}(undef, grid_axis)
    θ_grid     = fill(deg2rad(θ_B₀), grid_axis)

    comp_fac = 0.0
    for i in grid_axis
        if x_grid_cm[i] < 0cm
            uₓ_sk_grid[i] = u₀
            γ_sf_grid[i] = γ₀
            β_ef_grid[i] = 0.0
            γ_ef_grid[i] = 1.0
            btot_grid[i] = bmag₀
        else
            u = u₀ / r_comp
            β = u / c |> NoUnits
            uₓ_sk_grid[i] = u
            γ_sf_grid[i] = 1 / √(1 - β^2)
            β_ef_grid[i] = (β₀ - β) /  (1 - β₀*β)
            γ_ef_grid[i] = 1 / √(1 - β_ef_grid[i]^2)

            # When initializing magnetic field, include necessary corrections for turbulence compression
            z_comp         = (γ₀ * u₀) / (γ_sf_grid[i] * u)
            local comp_fac = 1 + (√((1 + 2*z_comp^2)/3) - 1) * bturb_comp_frac
            # Also include any additional amplification specified
            amp_fac        = 1 + (comp_fac - 1) * bfield_amp
            btot_grid[i]   = bmag₀ * amp_fac
        end
    end

    utot_grid = copy(uₓ_sk_grid) # uz_sk_grid is 0, so don't bother adding
    #utot_grid = uₓ_sk_grid + uz_sk_grid
    #utot_grid = hypot.(uₓ_sk_grid, uz_sk_grid)

    εB_grid = OffsetVector{Float64}(undef, grid_axis)
    # If directed in data_input, use a custom-defined ε_B to set btot_grid.
    #-------------------------------------------------------------------------------
    if use_custom_εB
        set_custom_εB!(εB_grid, btot_grid, grid_axis,
                       n_ions, species, bmag₀,
                       flux_px_UpS, flux_energy_UpS, uₓ_sk_grid, x_grid_rg,
                       comp_fac,
                       γ₀, β₀, u₀)
    else
        fill!(εB_grid, 1e-99)
    end
    #-------------------------------------------------------------------------
    # Custom εB_grid defined if needed

    bmag₂ = btot_grid[end]

    return (uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
            β_ef_grid, γ_ef_grid, btot_grid, θ_grid, εB_grid, bmag₂)
end

"""
    set_custom_εB!(...)

TODO
"""
function set_custom_εB!(
        εB_grid, btot_grid,
        grid_axis,
        n_ions, species, bmag₀,
        flux_px_UpS, flux_energy_UpS, uₓ_sk_grid, x_grid_rg,
        comp_fac,
        γ₀, β₀, u₀)

    @debug("Input parameters", εB_grid, btot_grid, grid_axis, n_ions, species, bmag₀,
           flux_px_UpS, flux_energy_UpS, uₓ_sk_grid, x_grid_rg, comp_fac, γ₀, β₀, u₀)

    # Calculate ε_B₀ which depends on far UpS magnetic field and mass density. If electrons
    # aren't a separate species, they don't contribute enough mass to be important
    n₀ = dot(density.(species), mass.(species))/mp # total number density
    εB₀ = bmag₀^2 / (8π * n₀ * E₀_proton)

    # The Monte Carlo length is rg₀ = γ₀ ⋅ β₀ ⋅ E₀_proton / (e ⋅ Bmag₀). The plasma skin
    # depth is λ_SD = γ₀ ⋅ E₀_proton / (4π ⋅ e² ⋅ den₀), where den₀ refers to the upstream
    # number density of electrons. With the definition
    #     σ = 2εB₀/γ₀ = Bmag₀² / (4π ⋅ γ₀ ⋅ n₀ ⋅ E₀_proton),
    # where n₀ here refers to the number density of *protons*, one can show that in
    # the shock frame (where grid exists),
    #     λ_SD = β₀ / √(σ ⋅ density_p/density_e) ⋅ rg₀.
    n₀_electron = density(species[end]) # electron number density
    σ = 2εB₀ / γ₀
    rg2sd = β₀ / √(σ*n₀/n₀_electron)

    # Also need the final value of ε_B downstream, in case our DwS region is long enough
    # that the magnetic field can decay to this value. Note that the R-H relations can be
    # rearranged to read
    #     energy_density(x) = F_en₀/u(x) - F_px₀
    # assuming flux conservation everywhere.
    energy_density₂ = (flux_energy_UpS + γ₀*u₀*n₀*E₀_proton) / uₓ_sk_grid[end] - flux_px_UpS
    εB₂ = (bmag₀*comp_fac)^2 / (8π * energy_density₂)
    # Use this value to compute the distance downstream at which the field will have decayed to it.
    # Per the Blandford-McKee solution, energy ∝ 1/χ ∝ 1/distance DwS. Since we do not actually
    # modify our pressures and densities according to the BM solution, instead modify εB
    end_decay_rg = (5e-3 / εB₂) / rg2sd

    @debug("Setting custom εB", n₀, εB₀, n₀_electron, σ, rg2sd, energy_density₂, εB₂, end_decay_rg)

    # Now we can calculate εB_grid and use it to calculate btot_grid. Per the Blandford-McKee
    # solution, energy ∝ 1/χ ∝ 1/distance DwS. Since we do not actually modify our pressures
    # and densities according to the BM solution, instead modify ε_B
    for i in grid_axis
        x_grid_sd = x_grid_rg[i]*rg2sd
        if x_grid_sd < -50
            εB_grid[i] = max(1.04e-5 / abs(x_grid_sd)^0.6, εB₀)
        elseif x_grid_sd < 50
            εB_grid[i] = 1e-4
        elseif x_grid_rg[i] < end_decay_rg
            εB_grid[i] = 5e-3 / x_grid_sd
        else
            εB_grid[i] = εB₂
        end
        energy_density = (flux_energy_UpS + γ₀*u₀*n₀*E₀_proton) / uₓ_sk_grid[i] - flux_px_UpS
        #@debug("Setting εB_grid array elements", i, εB_grid[i], energy_density)
        # FIXME this tries to be a square root of a negative number sometimes
        #btot_grid[i] = √(8π * εB_grid[i] * energy_density)
        btot_grid[i] = √abs(8π * εB_grid[i] * energy_density)
    end

    return εB_grid, btot_grid
end

"""
    init_pop(...)

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
- weight_in
- ptot_pf_in
- pb_pf_in
- x_PT_cm_in
- pxx_flux
- pxz_flux
- energy_flux
"""
function init_pop(
        do_fast_push, inp_distr, i_ion, m,
        # from controls module
        T₀_ion, energy_inj, inj_weight, n_pts_inj, n₀_ion, x_grid_start, rg₀, η_mfp,
        x_fast_stop_rg, β₀, γ₀, u₀, n_ions, m_ion,
        # from grid_vars module
        n_grid, x_grid_rg, uₓ_sk_grid, γ_sf_grid,
        # from iteration_vars module
        ptot_inj, weight_inj, n_pts_MB,
    )

    # If not using fast push, this procedure is quite quick.
    # Fill the arrays and return to the main loop
    if !do_fast_push
        if inp_distr == 1
            T_or_E = T₀_ion[i_ion]
        elseif inp_distr == 2
            T_or_E = energy_inj
        else
            error("not set to handle inp_distr > 2")
        end

        ptot_inj[:,i_ion], weight_inj[:,i_ion], n_pts_MB[i_ion] = set_inj_dist(
            inj_weight, n_pts_inj, inp_distr, T_or_E, m, n₀_ion[i_ion])

        n_pts_use = n_pts_MB[i_ion]
        weight_in = weight_inj[1:n_pts_use, i_ion]
        ptot_pf_in = ptot_inj[1:n_pts_use, i_ion]
        pb_pf_in = ptot_pf_in[1:n_pts_use] * 2*(rand(n_pts_use) .- 0.5)
        x_PT_cm_in = fill(x_grid_start - 10*rg₀*η_mfp, n_pts_use)
        pxx_flux = zeros(MomentumDensityFluxCGS, n_grid)
        pxz_flux = zeros(MomentumDensityFluxCGS, n_grid)
        energy_flux = zeros(EnergyDensityFluxCGS, n_grid)

        i_grid_in = zeros(Int, n_pts_use)
        return (n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in,
                x_PT_cm_in, pxx_flux, pxz_flux, energy_flux)
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
    @debug("Found i_stop", i_stop)

    relativistic = (β₀ ≥ β_rel_fl)
    density_ratio = u₀ / uₓ_sk_grid[i_stop]
    if relativistic
        density_ratio *= γ₀ / γ_sf_grid[i_stop]
    end

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal
    # gas, to calculate the far UpS internal energy   #assumecold
    Γ_sph = 5//3

    temp_ratio = density_ratio^Γ_sph / density_ratio

    if (kB * T₀_ion[i_ion] * temp_ratio) > (4 * m*c^2 * energy_rel_pt)
        error("Fast push cannot work because highest energy thermal particles become mildly relativistic. ",
              "Move fast push location UpS or disable entirely.")
    end
    pxx_flux = zeros(MomentumDensityFluxCGS, n_grid)
    pxz_flux = zeros(MomentumDensityFluxCGS, n_grid)
    energy_flux = zeros(EnergyDensityFluxCGS, n_grid)


    # Only run through the flux updates for the first particle species (i.e.
    # protons, not that it matters here); skip thereafter
    if i_ion == 1
        flux_update!(
            pxx_flux, pxz_flux, energy_flux,
            m_ion, n₀_ion, T₀_ion, relativistic,
            i_stop,
            γ₀, u₀, γ_sf_grid, uₓ_sk_grid,
        )
    end  # check on i_ion


    # With fast push fluxes taken care of, create the particle distribution
    # that will be injected at x_fast_stop
    T_or_E = T₀_ion[i_ion] * temp_ratio

    ptot_inj[:,i_ion], weight_inj[:,i_ion], n_pts_MB[i_ion] = set_inj_dist(
        inj_weight, n_pts_inj, inp_distr, T_or_E, m, n₀_ion[i_ion])
    n_pts_use = n_pts_MB[i_ion]

    weight_in  = weight_inj[1:n_pts_use, i_ion]
    ptot_pf_in = ptot_inj[1:n_pts_use, i_ion]
    x_PT_cm_in = fill(x_fast_stop_rg * rg₀, n_pts_use)
    i_grid_in  = fill(i_stop, n_pts_use)

    pb_pf_in = zeros(MomentumCGS, n_pts_use)
    for i_prt in 1:n_pts_use

        # Per Vladimirov+ (2009) [PhD], particle velocities should not be isotropic
        # in plasma frame. Must be weighted according to shock-frame pitch angle.
        #------------------------------------------------------------------------

        if relativistic
            γₚ_pf = hypot(1, ptot_pf_in[i_prt] / (m*c))
            vt_pf = ptot_pf_in[i_prt] / (γₚ_pf * m)

            vmin² = ((uₓ_sk_grid[i_stop] - vt_pf) / (1 - uₓ_sk_grid[i_stop]*vt_pf/c^2))^2
            vmax² = ((uₓ_sk_grid[i_stop] + vt_pf) / (1 + uₓ_sk_grid[i_stop]*vt_pf/c^2))^2
            # Unitful distributions not yet supported :(
            dist_v_sf = Uniform(ustrip(cm^2/s^2, vmin²), ustrip(cm^2/s^2, vmax²))

            vx_sf = cm/s * √(rand(dist_v_sf))
            vx_pf = (vx_sf - uₓ_sk_grid[i_stop]) / (1 - vx_sf*uₓ_sk_grid[i_stop]/c^2)
        else
            γₚ_pf = 1.0
            vt_pf = ptot_pf_in[i_prt] / m

            vmin² = (uₓ_sk_grid[i_stop] - vt_pf)^2
            vmax² = (uₓ_sk_grid[i_stop] + vt_pf)^2
            # Unitful distributions not yet supported :(
            dist_v_sf = Uniform(ustrip(cm^2/s^2, vmin²), ustrip(cm^2/s^2, vmax²))

            vx_sf = cm/s * √(rand(dist_v_sf))
            vx_pf = vx_sf - uₓ_sk_grid[i_stop]
        end

        pb_pf_in[i_prt] = γₚ_pf * m * vx_pf
        #------------------------------------------------------------------------
        # Velocity-weighted pitch angles finished
    end

    return n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux
end


"""
    flux_update!(...)

TODO
"""
function flux_update!(
        pxx_flux, pxz_flux, energy_flux,
        m_ion, n₀_ion, T₀_ion, relativistic,
        i_stop,
        γ₀, u₀, γ_sf_grid, uₓ_sk_grid
    )
    # Update the flux arrays as if the particles had actually crossed them
    #-----------------------------------------------------------------------
    # Obtain UpS pressure and density
    P₀ = dot(n₀_ion, T₀_ion) * kB
    ρ₀ = dot(n₀_ion, m_ion)

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal gas,
    # to calculate the far UpS pressure and internal energy   #assumecold
    Γ_sph = 5//3

    # Calculate fluxes and update arrays; note that if fast push isn't enabled
    # i_stop = 0, and this loop never executes
    for i in 1:i_stop

        density_ratio = (γ₀ * u₀) / (γ_sf_grid[i] * uₓ_sk_grid[i])
        ρ_curr        = ρ₀ * density_ratio

        # Note assumption that Γ_sph doesn't change from zone to zone: #assumecold
        pressure_curr = P₀ * density_ratio^Γ_sph

        β_curr   = uₓ_sk_grid[i] / c
        γβ_curr = γ_sf_grid[i] * uₓ_sk_grid[i] / c

        # Determine fluxes while handling different possible orientations and shock speeds.
        # For non-relativistic fluxes, expand out to β^2 to allow for better
        # matching with relativistic versions.
        # WARNING: these fluxes do not include contributions from a strong
        # magnetic field. This is incorporated during the smoothing process.
        #----------------------------------------------------------------------
        flux_pz = 0.0erg/cm^3
        if !relativistic
            flux_pₓ = (ρ_curr * uₓ_sk_grid[i]^2 * (1 + β_curr^2)
                       + pressure_curr * (1 + Γ_sph/(Γ_sph-1) * β_curr^2))
            flux_energy = (ρ_curr/2 * uₓ_sk_grid[i]^3 * (1 + 1.25*β_curr^2)
                           + pressure_curr * uₓ_sk_grid[i] * Γ_sph/(Γ_sph-1) * (1 + β_curr^2))
        else
            e_curr = ρ_curr*c^2 # energy density

            flux_pₓ = pressure_curr + γβ_curr^2 * (e_curr + Γ_sph/(Γ_sph-1)*pressure_curr)
            flux_energy = (γβ_curr^2 * c / (uₓ_sk_grid[i]/c) *
                           (e_curr + Γ_sph/(Γ_sph-1)*pressure_curr)
                           # Subtract mass-energy flux from flux_energy to bring
                           # it in line with non-relativistic calculations
                           - γβ_curr*c * e_curr)
        end
        #--------------------------------------------------------------------
        # Fluxes calculated

        @debug "Updating fluxes" flux_pₓ flux_pz flux_energy pxx_flux pxz_flux energy_flux
        pxx_flux[i] = flux_pₓ
        pxz_flux[i] = flux_pz
        energy_flux[i] = flux_energy

    end  # loop over grid location
    #------------------------------------------------------------------------
    # Arrays updated through i_fast_stop
end

"""
    set_inj_dist(inj_weight, n_pts_inj, inp_distr, T_or_E, m, n₀)

Sets the injected particle distributions for all species. Initially, particles are placed in
a Maxwell-Boltzmann (thermal) distribution based on supplied temperature and particle mass.
This is corrected at the end if a δ-function distribution was requested (not worried about
the wasted computation because this subroutine runs only rarely).

### Arguments
FIXME
- `inj_weight`: whether each particle or each bin has equal weights (`true` for equal weight particles,
  `false` for equal weight bins)
- `n_pts_inj`: target # of particles for distribution
- `inp_distr`: thermal, δ-function, or some other distribution
- `T_or_E`: if using thermal distribution, this is temperature[K]; if δ function, it's injection energy[keV]
- `m`: mass for this particle species
- `n₀`: far UpS number density for this species

### Returns
- `ptot_out`: array holding plasma frame total momenta for all particles in the distribution
- `weight_out`: array holding particle weights
- `n_pts_use`: number of particles in the distribution; will almost surely be different from
  `n_pts_inj` if using thermal distribution

CHECKTHIS: that output distribution matches M-B, just to make sure I haven't made a typo
"""
function set_inj_dist(inj_weight::Bool, n_pts_inj, inp_distr, T_or_E, m, n₀)

    # Error prevention
    0 < inp_distr < 3 || throw(DomainError(inp_distr, "Code can only do inp_distr = 1 or 2."))

    # Administrative constants, e.g. total number of particles to distribute
    #--------------------------------------------------------------------------
    if inj_weight
        # Can't say anything about total number of particles in this case, b/c
        # haven't split them into M-B distribution yet
    else
        n_per_bin = n_pts_inj ÷ num_therm_bins  # Integer math loses excess particles, but that's intended behavior
        n_pts_tot = n_per_bin * num_therm_bins

        if n_per_bin < 5
            throw(ArgumentError("too few particles per bin ($n_per_bin; need at least 5). Increase n_pts_inj."))
        end
    end
    #------------------------------------------------------------------------
    # End administrative section


    # Set a mess of constants to be referred to routinely
    #------------------------------------------------------------------------
    rm_energy = m * c^2

    kT      = kB * T_or_E  # Working under assumption of thermal dist now
    # Minimum, maximum extent of Maxwell-Boltzmann distribution
    kT_min  = 2e-3 * kT
    kT_max  = 10   * kT

    kT_rel_div = energy_rel_pt
    relativistic = (kT/rm_energy ≥ kT_rel_div)

    # Find min and max momenta of M-B curve
    if !relativistic
        # In non-relativistic case, kinetic energy ≈ thermal energy
        p_min = √(2m * kT_min)
        p_max = √(2m * kT_max)
    else
        # Once particles are relativistic, rest-mass energy becomes important:
        #     E² = p²c² + m²c⁴  ≈  (kT + mc²)²
        p_min = √((kT_min + rm_energy)^2 - rm_energy^2) / c
        p_max = √((kT_max + rm_energy)^2 - rm_energy^2) / c
    end

    Δp = (p_max - p_min) / num_therm_bins
    p_range = range(start = p_min, step = Δp, length = num_therm_bins+1)
    # define energy over kT
    if !relativistic
        E_range = p_range.^2 / (2m * kT)
    else
        E_range = hypot.(p_range*c, rm_energy) / kT
    end
    #------------------------------------------------------------------------
    # End of constants section


    # Generate the Maxwell-Boltzmann distribution. The actual calculation of f1 and f2 does
    # not need modification for relativistic momenta. Such a modification would only affect
    # the normalization of the curve, not the dependence on momentum for a particular value
    # of kT. Since we only care about the relative fractional area of each bin, overall
    # normalization doesn't matter.
    #------------------------------------------------------------------------
    # Find total area under M-B curve
    area_tot = 0.0g*cm/s
    for (i, p1) in enumerate(p_range[begin:end-1])
        p2 = p_range[i+1]

        energy_o_kT1 = E_range[i]
        energy_o_kT2 = E_range[i+1]

        # Start working in log space because of potentially huge exponents
        f1 = exp(2log(p1/(g*cm/s)) - energy_o_kT1)
        f2 = exp(2log(p2/(g*cm/s)) - energy_o_kT2)
        #@debug "In log space" i p1 p2 f1 f2

        area_tot += Δp * 0.5 * (f1 + f2) # Integrate using the trapezoid rule
    end

    ptot_out = zeros(MomentumCGS, na_particles)
    weight_out = zeros(na_particles)

    # Fill the bins with particles. Needs to be done in two separate loops, despite the
    # similarities, since particle weights are handled differently based on value of inj_weight
    if inj_weight # First, particles have equal weight
        n_pts_tot = set_inj_dist_particle_equal_weight!(
            ptot_out, weight_out, p_range, E_range, area_tot, n_pts_inj, Δp, n₀)
    else # Bins have equal weight
        set_inj_dist_bin_equal_weight!(
            ptot_out, weight_out, p_range, E_range, area_tot, Δp, n₀)
    end  # test of inj_weight

    n_pts_use = n_pts_tot
    #------------------------------------------------------------------------
    # End of Maxwell-Boltzmann section


    # If the particles are to use a δ-function distribution, set it here
    #--------------------------------------------------------------------------
    if inp_distr == 2
        n_pts_use = n_pts_inj

        rm_energy      = m * c^2
        energy_inj_cgs = ustrip(erg, T_or_E*keV)

        if energy_inj_cgs/rm_energy < energy_rel_pt
            p1 = √(2m * energy_inj_cgs)
        else
            p1 = √(energy_inj_cgs^2 - rm_energy^2) / c
        end

        ptot_out[1:n_pts_inj] .= p1
        weight_out[1:n_pts_inj] .= ustrip(cm^-3, n₀) / n_pts_tot

    end  # check of inp_distr
    #--------------------------------------------------------------------------
    # δ-function handled


    return (ptot_out, weight_out, n_pts_use)
end

function set_inj_dist_particle_equal_weight!(
        ptot_out, weight_out,
        p_range, E_range, area_tot, n_pts_inj, Δp, n₀)

    area_per_pt = area_tot / n_pts_inj # Area each particle gets if inj_weight = T
    n_pts_tot   = 1                    # Total number of particles in M-B distribution

    for (i, p1) in enumerate(@view p_range[begin:end-1])
        p2 = p_range[i+1]

        energy_o_kT1 = E_range[i]
        energy_o_kT2 = E_range[i+1]

        # Work in log space because of potentially huge exponents
        f1 = exp(2log(p1/(g*cm/s)) - energy_o_kT1)
        f2 = exp(2log(p2/(g*cm/s)) - energy_o_kT2)

        bin_area = Δp * (f1 + f2)/2 # Calculate bin area using trapezoid rule

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
    weight_out[1:n_pts_tot] .= ustrip(cm^-3, n₀) / n_pts_tot

    return n_pts_tot
end

"""
    set_inj_dist_bin_equal_weight!(...)

TODO
"""
function set_inj_dist_bin_equal_weight!(
        ptot_out, weight_out,
        p_range, E_range, area_tot, Δp, n₀)
    for (i, p1) in enumerate(p_range[begin:end-1])
        p2 = p_range[i+1]

        energy_o_kT1 = E_range[i]
        energy_o_kT2 = E_range[i+1]

        # Work in log space because of potentially huge exponents
        f1 = exp(2log(p1) - energy_o_kT1)
        f2 = exp(2log(p2) - energy_o_kT2)

        bin_area = Δp * (f1 + f2)/2 # Calculate bin area using trapezoid rule

        area_frac = bin_area / area_tot

        # particles in this bin will have momentum = the geometric center of the bin
        jstart = (i-1)*n_per_bin + 1
        jend   = jstart + (n_per_bin-1)
        ptot_out[jstart:jend] .= √(p1*p2)
        weight_out[jstart:jend] .= area_frac / n_per_bin * n₀

    end  # loop over i
end
end # module
