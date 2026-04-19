module initializers
export calc_downstream, calc_rRH
export upstream_fluxes, upstream_machs
export set_psd_mom_bins, init_pop, setup_profile, setup_grid, set_psd_angle_bins

using Random: rand
using OffsetArrays: OffsetVector, Origin
using LinearAlgebra: dot
using StaticArrays: SVector
using Roots: find_zero, Newton
using Unitful, UnitfulAstro
using Unitful: g, cm, s, dyn, erg, keV
using Unitful: mp, c, k as kB
using Distributions: TriangularDist
using ..parameters: num_therm_bins, na_particles, E_rel_pt, β_rel_fl
import ..density, ..temperature, ..mass, ..number_density
using ..CGSTypes: MomentumCGS, BFieldCGS, MomentumDensityFluxCGS, EnergyDensityFluxCGS

"""
    calc_downstream(...)

Uses the Rankine-Hugoniot jump conditions to calculate the downstream conditions for a test
particle shock. Big difference between this subroutine and `calc_rRH` is that we already know
what the downstream speed is, courtesy of `r_comp` in the input.

Shock must be parallel (not oblique).

### Arguments

- `B₀`: far upstream magnetic field strength [Gauss]
- `r_comp`: compression ratio of shock
- `β₀`: upstream bulk fluid speed, over c

### Returns

- `β`: (total) bulk fluid speed downstream
- `γ`: Lorentz factor associated with `β₂`
- `B`: downstream magnetic field strength [Gauss]
- `θ_B`: angle[deg] between downstream magnetic field and shock normal
- `θᵤ`: angle[deg] between downstream fluid velocity and shock normal
"""
function calc_downstream(B₀, r_comp, β₀)
    β = β₀ / r_comp
    γ = 1 / √(1 - β^2)
    B = B₀
    θ_B = 0.0
    θᵤ = 0.0
    return β, γ, B, θ_B, θᵤ
end

"""
    calc_rRH(...)

Uses the Rankine-Hugoniot jump conditions to calculate the compression ratio for a shock
assuming test-particle conditions. In other words, (1) sharp shock, (2) negligible/no DSA,
and (3) no escaping flux. Additionally assumes that the inflowing plasma has
non-relativistic thermal speeds to make upstream adiabatic index exactly 5/3.

### Arguments
- `β₀`, `γ₀`: shock speed and Lorentz factor
- `n_ions`: number of different ion species
- `m_ion`: array of particle species' mass
- `n₀_ion`: array of far upstream densities for particle species
- `T₀_ion`: array of particle species' far upstream temperatures

Shock must be parallel (not oblique)

### Returns
- `r_RH`: Rankine-Hugoniot compression ratio
- `Γ₂_RH`: ratio of specific heats (adiabatic index) for downstream region, assuming `r_comp` = `r_RH`
"""
function calc_rRH(β₀::Real, γ₀::Real, species)

    # Two possibilities for R-H relations: nonrelativistic/relativistic
    # Cutoff for (non-)relativistic is set in module 'parameters'
    relativistic = (β₀ < β_rel_fl)

    # Calculate thermal pressure of far upstream gas
    P₀ = dot(density.(species), temperature.(species)) * kB
    ρ₀ = dot(density.(species), mass.(species))

    if !relativistic    #  Possibility 1: Nonrelativistic, parallel
        r_RH, Γ₂_RH = calc_rRH_nonrelativistic(P₀, ρ₀, β₀)
    else                #  Possibility 2: Relativistic, parallel
        r_RH, Γ₂_RH = calc_rRH_relativistic(species, ρ₀, P₀, (β₀, γ₀))
    end

    return r_RH, Γ₂_RH
end

"""
    calc_rRH_nonrelativistic(P₀, ρ₀, β₀)

Nonrelativistic, parallel.

Solution comes from Ellison (1985) [1985JGR....90...29E].
Uses far upstream Mach number to calculate `r_RH`.
"""
function calc_rRH_nonrelativistic(P₀, ρ₀, β₀::Real)

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal
    # gas, to calculate the far upstream sound speed and Mach number   #assumecold
    Γ_sph = 5//3
    cₛ = √(Γ_sph * P₀ / ρ₀)
    M_Z = β₀ * Unitful.c / cₛ |> NoUnits

    # Finally, use Equation (11) from Ellison (1985) to calculate r_RH.
    # Note that q = 0 here because we assume no escaping flux. This
    # simplifies the denominator quite a bit from the equation.
    r_RH = 8 / (2 + 6 / M_Z^2)

    # In non-relativistic case, downstream adiabatic index is pegged to 5/3
    Γ₂_RH = 5//3

    return r_RH, Γ₂_RH
end

"""
    calc_rRH_relativistic(species, ρ₀, P₀, β₀, n₀_ion)

Relativistic, parallel.

Solution comes from [Ellison & Reynolds (1990) 1991ApJ...378..214E](@cite ellison_determination_1991).
Uses relativistic Rankine-Hugoniot relations. See that paper for
details of equations and associated quantities. Briefly,
```math
\\begin{align*}
    γ_1   n_1 b_1          &=  γ_2  n_2 b_2           \\tag{R-H1} \\\\
    γ_1^2 w_1 b_1^2 + P_1  &=  γ_2^2 w_2 b_2^2 + P_2  \\tag{R-H2} \\\\
    γ_1^2 w_1 b_1          &=  γ_2^2 w_2 b_2          \\tag{R-H3}
\\end{align*}
```
where
- ``w   = e_0 + e_K + P``   is the enthalpy as total energy density + pressure
- ``e_0 = n m c^2``         is the rest mass energy density
- ``e_K = n m c^2 (γ - 1)`` is the kinetic energy density, with ``γ = \\sqrt{1 + (p/mc)^2}``
- ``P   = \\frac{1}{3} n p v``  is the pressure

Assumes that downstream particle distributions are δ-functions.
Solves for `p₂` using Newton's method, then works backwards to `r_RH`.
"""
function calc_rRH_relativistic(species, ρ₀, P₀, β₀::Real, γ₀::Real)

    # Calculate two quantities to be used during loop to find r_RH: the
    # rest energy of each species, and the (number) density relative to protons
    n₀_ion = density.(species)          # number density of each species
    m_ion = mass.(species)              # rest mass of each species
    e₀_ion = dot(m_ion, n₀_ion) * c^2   # rest energy density of each species
    relative_ion_energy = e₀_ion / first(n₀_ion)

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal gas,
    # to calculate the far upstream enthalpy      #assumecold
    Γ_sph = 5//3
    w₀ = ρ₀ * c^2 + Γ_sph/(Γ_sph-1) * P₀

    # Calculate the far upstream momentum flux
    upstream_mom_flux = γ₀^2 * w₀ * β₀^2 + P₀
    upstream_num_flux = γ₀ * first(n₀_ion) * β₀ # Protons only here; not strictly correct but appropriate for later use

    function F(p)
        γβ = p / (mp * c) # since p is relativistic, p/mc = γmv/mc = γ⋅β
        γ = √(1 + γβ^2)
        P = relative_ion_energy / 3 * γβ^2 / γ # pressure
        w = relative_ion_energy * (γ + γβ^2 / 3γ)
        return upstream_num_flux / γ * (w * γ^2 + P) - upstream_mom_flux
    end

    # Now use Newton's method to determine the downstream momentum that satisfies the
    # R-H relations. Assumptions: (1) momentum distribution functions are δ-functions
    # rather than thermal, and (2) any non-proton species have p ∝ m
    p₂_found = find_zero(F, 0, Newton())

    # Calculate the compression ratio β₀/β₂ associated with p₂_found
    γβ = p₂_found / (mp * c)  # Protons only here because of how the math in w_fac works out

    # Pressure, internal energy, and enthalpy with proton density factored out
    P_fac = relative_ion_energy / 3 * γβ^2 / √(1 + γβ^2)
    e_fac = relative_ion_energy * (√(1 + γβ^2) - 1)
    w_fac = relative_ion_energy * (√(1 + γβ^2) + 1 / 3 * γβ^2 / √(1 + γβ^2))

    # Calculate adiabatic index downstream
    Γ₂_RH = 1 + P_fac / e_fac

    # Finally, get downstream speed and compression ratio
    β₂ = √(1 - (n₀_ion[1] * w_fac / γ₀ * w₀)^2)  # Using only proton density is correct


    r_RH = β₀ / β₂
    #------------------------------------------------------------------------
    # r_RH found using Newton's method

    return r_RH, Γ₂_RH
end

"""
    set_psd_mom_bins(...)

Sets the BOUNDARIES of the bins of the phase space distribution. The bins are numbered from
0 to `num_psd_mom_bins`, each boundary denotes the lower edge of that # bin; the indices thus
run from 0 to `num_psd_mom_bins + 1`.

Logarithmic spacing over all decades from `E_min` to `E_max` is used.
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
    num_psd_mom_bins = trunc(Int, log10(psd_mom_max / psd_mom_min) * psd_bins_per_dec_mom)
    num_psd_mom_bins += 2   # Add two extra bins just to be safe

    # Fill in the array psd_mom_bounds, remembering that the array holds LOWER boundaries
    # of that bin
    psd_mom_bounds = Origin(0)(Float64[-99.0])
    append!(
        psd_mom_bounds,
        range(
            start = log10(psd_mom_min / (mp * c)),
            step = 1 / psd_bins_per_dec_mom,
            length = num_psd_mom_bins + 1
        )
    )

    length(psd_mom_bounds) == num_psd_mom_bins + 2 || error() # sanity check

    return num_psd_mom_bins, psd_mom_bounds
end

"""
    set_psd_angle_bins(...)

Sets the BOUNDARIES of the bins of the phase space distribution. The bins are numbered from
0 to `num_psd_θ_bins`, each boundary denotes the lower edge of that # bin; the indices thus
run from 0 to `num_psd_θ_bins + 1`.

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

    psd_log_θ_bins = trunc(Int, log10(psd_θ_fine / psd_θ_min) * psd_bins_per_dec_θ)
    #num_psd_θ_bins = psd_log_θ_bins + psd_lin_cos_bins

    # Fill the logarithmic part of psd_θ_bounds using the angle (in radians), NOT its logarithm
    psd_θ_bounds = Origin(0)(Float64[1.0e-99])
    append!(psd_θ_bounds, psd_θ_min * ten_root_θ .^ range(0, length = psd_log_θ_bins))

    # Now fill in the linear part of psd_θ_bounds.
    # Note that the lower boundary of the first cell is psd_cos_fine
    Δcos = (psd_cos_fine + 1) / psd_lin_cos_bins
    append!(psd_θ_bounds, range(start = psd_cos_fine, step = -Δcos, length = psd_lin_cos_bins + 1))

    sort!(psd_θ_bounds)

    #return num_psd_θ_bins, Δcos, psd_θ_bounds
    return Δcos, psd_θ_bounds
end

"""
    set_photon_shells(...)

If photon calculation is desired, photons will be collected into upstream and downstream shells for
easier viewing. This subroutine sets the endpoints of the shells, as well as their midpoints.

Because the particle spectrum changes most rapidly near the shock, zones should be small
near the shock and get larger as you move further away.

First, divide the domain between the shock and the FEB (either upstream or downstream into
n sections, where n is the respective number of shells to use. HOWEVER, do this by exponent,
ranging from -1 to log10(|x_FEB|).

Keep track of the boundaries of each shells, since we will need those for calculating the
total number of particles emitting when we get to that point in photon production. The
midpoints are less useful, since photons are calculated on a zone-by-zone basis rather than
just at select points in the shock profile. Keep them anyway, since the memory overhead is low.
"""
function set_photon_shells(
        num_upstream_shells, num_downstream_shells,
        use_prp, feb_upstream, feb_downstream, rg₀, x_grid_stop_rg,
    )

    total_shells = num_upstream_shells + num_downstream_shells
    x_shell_midpoints = zeros(total_shells)
    x_shell_endpoints = zeros(total_shells + 1) # fencepost problem

    # Handle upstream shells first
    set_upstream_photon_shells!(x_shell_midpoints, x_shell_endpoints, num_upstream_shells, feb_upstream, rg₀)

    # And repeat the process for the downstream shells.
    set_downstream_photon_shells!(
        x_shell_midpoints, x_shell_endpoints, num_upstream_shells, num_downstream_shells,
        use_prp, feb_downstream, rg₀, x_grid_stop_rg
    )

    # Convert from units of rg₀ to cm
    @. x_shell_endpoints *= rg₀
    return (x_shell_midpoints, x_shell_endpoints)
end

"""
    set_upstream_photon_shells!(...)

TODO
"""
function set_upstream_photon_shells!(
        x_shell_midpoints, x_shell_endpoints,
        num_upstream_shells, feb_upstream, rg₀,
    )
    x_section_width = (log10(abs(feb_upstream / rg₀)) + 1) / num_upstream_shells
    for i in 1:num_upstream_shells
        # Calculate upstream and downstream endpoints of each region,
        # as well as midpoint in log space
        if i == 1
            # Special case when i = 1, since our region starts at the shock.
            x_region_start = 0.0
            x_region_end = exp10(-1 + x_section_width)
            x_region_mid = exp10(-1 + x_section_width / 2)
        else
            # In the general case, note that x_region_start should be the same as
            # the previous region's x_region_end. This can be checked with print
            # or write statements at runtime.
            x_region_start = exp10(-1 + x_section_width * (i -   1) )
            x_region_end   = exp10(-1 + x_section_width *  i        )
            x_region_mid   = exp10(-1 + x_section_width * (i - 1/2) )
        end

        # Update the arrays with this information, remembering that the eventual array will
        # count downstream from the upstream FEB (so some array index juggling is necessary)
        # Also, add in the factor of -1 here, since upstream coordinates should be negative
        # in the MC code.
        N = num_upstream_shells + 1 - i
        x_shell_midpoints[N] = -x_region_mid
        x_shell_endpoints[N] = -x_region_end
        x_shell_endpoints[N + 1] = -x_region_start
    end
    return
end
"""
    set_downstream_photon_shells!(...)

TODO
"""
function set_downstream_photon_shells!(
        x_shell_midpoints, x_shell_endpoints,
        num_upstream_shells, num_downstream_shells, use_prp, feb_downstream, rg₀, x_grid_stop_rg,
    )
    # The downstream limit is set differently if using a PRP or a FEB
    limitdownstream = use_prp ? x_grid_stop_rg : feb_downstream / rg₀
    x_section_width = (log10(limitdownstream) + 1) / num_downstream_shells

    for i in 1:num_downstream_shells
        # Calculate upstream and downstream endpoints of each region, as well as midpoint in log space
        if i == 1
            # Special case when i = 1, since our region starts at the shock.
            x_region_start = 0.0
        else
            # In the general case, note that x_region_start should be the same as the previous
            # region's x_region_end. This can be checked with print statements at runtime.
            x_region_start = exp10(-1 + x_section_width * (i - 1))
        end
        x_region_mid = exp10(-1 + x_section_width * (i - 1 // 2))
        x_region_end = exp10(-1 + x_section_width * i)

        # Update the arrays with the information. Less index juggling here
        x_shell_endpoints[num_upstream_shells + i] = x_region_start
        x_shell_midpoints[num_upstream_shells + i] = x_region_mid
        x_shell_endpoints[num_upstream_shells + i + 1] = x_region_end
    end
    return
end

# Many grid zones are set manually; zones can easily be added/removed, but
# make sure to change the number in the log-spaced regions upstream or downstream
const FIRST_ZONE = SVector(
    -9.0, -8.0, -7.0, -6.0, -5.0, -4.5, -4.0, -3.5, -3.0,
    -2.5, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0,
    -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2,
    -0.15, -0.1,
    -0.07, -0.05, -0.04, -0.03, -0.02, -0.015, -0.01,
    -3.0e-3, -1.0e-3
)
# Extremely fine spacing right around the shock
const EXTREMELY_FINE_SPACING = SVector(-1.0e-4, -1.0e-7, 0.0, 1.0e-7, 1.0e-4)

# Downstream from the shock, spacing doesn't need to be quite so fine
# because velocity gradients aren't as extreme, if they exist at all
const DOWNSTREAM_SPACING = SVector(
    1.0e-3, 1.0e-2, 2.0e-2, 3.0e-2, 5.0e-2, 7.0e-2, 0.1,
    0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0
)

"""
    setup_grid(...)

### Arguments
- `x_grid_start_rg`
- `x_grid_stop_rg`
- `use_prp`
- `feb_downstream`
- `rg₀`

### Returns
- `x_grid_rg`
- `x_grid_start`
- `x_grid_stop`
"""
function setup_grid(x_grid_start_rg::Float64, x_grid_stop_rg::Float64, use_prp::Bool, feb_downstream::Float64, rg₀::Float64)

    # Recall that rg₀ is the gyroradius of a proton with speed u₀ in magnetic field B₀.

    # Set the start and stop positions in units of rg₀
    x_grid_start = x_grid_start_rg * rg₀
    if !use_prp
        x_grid_stop = feb_downstream
        @info("downstream FEB set at x = $x_grid_stop_rg rg₀. Overwriting entered value for 'XGDDW'.")
    else
        x_grid_stop = x_grid_stop_rg * rg₀
    end

    # Logarithmically-spaced grid zones run from x_grid_start_rg to -10rg₀. Set them here.
    n_log_upstream = 27
    Δlogx = (log10(-x_grid_start_rg) - 1) / n_log_upstream - 1

    # we build this chunk-by-chunk, which is inefficient
    x_grid_rg = Origin(0)(Float64[])

    push!(x_grid_rg, -1.0e30) # set left boundary of grid

    log_x_grid_upstream = range(start = log10(-x_grid_start_rg), step = -Δlogx, length = n_log_upstream)
    append!(x_grid_rg, -exp10.(log_x_grid_upstream))
    append!(x_grid_rg, FIRST_ZONE)
    append!(x_grid_rg, EXTREMELY_FINE_SPACING)
    append!(x_grid_rg, DOWNSTREAM_SPACING)

    # As seen above, the manually-set grid zones end at x = +1rg₀.
    # Downstream from there, more log-spaced zones.
    n_log_downstream = 16
    x_end_man = x_grid_rg[end]
    Δlogx = (log10(x_grid_stop*rg₀) - log10(x_end_man)) / n_log_downstream

    log_x_grid_downstream = range(start = log10(x_end_man), step = Δlogx, length = n_log_downstream)
    append!(x_grid_rg, exp10.(log_x_grid_downstream))

    push!(x_grid_rg, 1.0e30)  # set right boundary of the grid

    return x_grid_rg, x_grid_start, x_grid_stop
end

"""
    upstream_fluxes(...)

Calculates the far upstream fluxes for the shock.

Two different cases considered:

1. Non-relativistic oblique shock. Uses equations of Ellison+ (1996) [1996ApJ...473.1029E]
2. Relativistic shock, any obliquity. Uses equations of Double+ (2004) [2004ApJ...600..485D]

Only oblique equations used because they reduce trivially to parallel cases when θ_B₀ = 0.

HOWEVER, assumes that z-component of far upstream velocity is 0 in all cases; in practice
oblique shocks would induce some z-velocity in the shock profile even though particles
initially arrive with no bulk z component. Also assumes isotropic initial pressure, so
no off-diagonal components in pressure tensor.

Shock must be parallel (not oblique). Set `oblique` to `false`.

### Arguments

- `n₀_ion`
- `T₀_ion`
- `m_ion`
- `B₀`
- `θ_B₀`
- `u₀`, `β₀`, `γ₀`

No inputs; pulls everything from module 'controls'

### Returns
- `F_px_upstream`: far upstream momentum flux, x component
- `F_pz_upstream`: far upstream momentum flux, z component
- `F_energy_upstream`: far upstream energy flux
"""
function upstream_fluxes(n₀_ion, T₀_ion, m_ion, B₀, θ_B₀, u₀, β₀, γ₀)

    # Upstream internal energy density and pressure, assuming isotropic particle distribution.
    # Note that this INCLUDES the mass-energy density, which is typically omitted in
    # nonrelativistic calculations
    P₀ = dot(n₀_ion, T₀_ion) * kB |> dyn / cm^2 # pressure
    ρ₀ = dot(n₀_ion, m_ion)       |> g / cm^3   # mass density
    @debug("Calculated params", P₀, ρ₀)

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal gas,
    # to calculate the far upstream internal energy             #assumecold
    Γ_sph = 5//3
    # internal energy density
    e₀ = ρ₀ * c^2 + 1 / (Γ_sph - 1) * P₀ |> erg / cm^3

    # Quantities related to the upstream magnetic field. Note that B_z is the
    # z-component of the magnetic field, not B₀
    B_x = B₀ * cosd(θ_B₀)
    B_z = B₀ * sind(θ_B₀)

    relativistic = β₀ ≥ β_rel_fl

    if relativistic
        regime = Val(:relativistic)
        F_px_upstream, F_pz_upstream = upstream_momentum_flux(regime, β₀, γ₀, e₀, P₀, B₀, B_x, B_z)
        F_energy_upstream = upstream_energy_flux(regime, u₀, β₀, γ₀, e₀, ρ₀, P₀, B_z)
    else
        regime = Val(:classical)
        # Non-relativistic version. Note that it's missing the ρc² flux present in the
        # relativistic forms above. It is also expanded to second order in β₀ (only in the
        # hydro terms, for now) to allow for more precise matching with the relativistic version
        F_px_upstream, F_pz_upstream = upstream_momentum_flux(regime, u₀, β₀, ρ₀, P₀, B_x, B_z, Γ_sph)
        F_energy_upstream = upstream_energy_flux(regime, u₀, β₀, ρ₀, P₀, B_z, Γ_sph)
    end

    return (F_px_upstream, F_pz_upstream, F_energy_upstream)
end

"""
    upstream_momentum_flux(regime=Val(:classical), u₀, β₀, ρ₀, P₀, B_x, B_z, Γ)
    upstream_momentum_flux(regime=Val(:relativistic), β₀, γ₀, e₀, P₀, B₀, B_x, B_z)

### Arguments
- `regime`: Whether to use Newtonian or relativistic calculations
- `u₀`, `β₀`, `γ₀`
- `e₀`
- `ρ₀`
- `P₀`
- `B₀`, `B_x`, `B_z`
- `Γ`: Adiabatic index
"""
function upstream_momentum_flux end
function upstream_momentum_flux(::Val{:classical}, u₀, β₀, ρ₀, P₀, B_x, B_z, Γ)
    Ξ_sph = Γ / (Γ - 1)
    u_B = B_z^2 / 8π                # magnetic field energy density
    F_px_upstream = ρ₀ * u₀^2 * (1 + β₀^2) + P₀ * (1 + Ξ_sph * β₀^2) + u_B
    F_pz_upstream = - B_x * B_z / 4π
    return F_px_upstream, F_pz_upstream
end
function upstream_momentum_flux(::Val{:relativistic}, β₀, γ₀, e₀, P₀, B₀, B_x, B_z)

    # Momentum flux, x-component
    # Fluid part (Double+ Eq 23)
    F_pₓ_fl = (γ₀ * β₀)^2 * (e₀ + P₀) + P₀ |> g / (cm * s^2)
    # EM part (Double+ Eq 25)
    F_pₓ_EM = γ₀^2 * ((β₀ * B₀)^2 + B_z^2 - B_x^2) / 8π |> g / (cm * s^2)
    @debug("Found partial fluxes", F_pₓ_fl, F_pₓ_EM)
    F_px_upstream = F_pₓ_fl + F_pₓ_EM                         # Total

    # Momentum flux, z-component (Fluid Part = 0, from Double+ Eq 24)
    # Total = EM part (Double+ Eq 26)
    F_pz_upstream = -γ₀ * B_x * B_z / 4π

    return F_px_upstream, F_pz_upstream
end

"""
    upstream_energy_flux(regime=Val(:classical), u₀, β₀, ρ₀, P₀, Γ, B_z)
    upstream_energy_flux(regime=Val(:relativistic), u₀, β₀, γ₀, e₀, ρ₀, P₀, B_z)

### Arguments
- `regime`: Whether to use Newtonian or relativistic calculations
- `u₀`, `β₀`, `γ₀`
- `e₀`
- `ρ₀`
- `P₀`
- `B₀`, `B_x`, `B_z`
- `Γ`: Adiabatic index
"""
function upstream_energy_flux end
function upstream_energy_flux(::Val{:classical}, u₀, β₀, ρ₀, P₀, B_z, Γ)
    Ξ = Γ / (Γ - 1)
    return (
        ρ₀ * u₀^3 * (1 + 1.25 * β₀^2) / 2
        + P₀ * u₀ * Ξ * (1 + β₀^2)
        + u₀ * B_z^2 / 4π
    )
end
function upstream_energy_flux(::Val{:relativistic}, u₀, β₀, γ₀, e₀, ρ₀, P₀, B_z)
    F_energy_fl = γ₀^2 * β₀ * (e₀ + P₀)     # Fluid part (Double+ Eq 20)
    F_energy_EM = γ₀^2 * β₀ * B_z^2 / 4π    # EM part (Double+ Eq 21)
    # Total -- convert to cgs units!
    F_energy_upstream = c * (F_energy_fl + F_energy_EM)

    # And subtract off the mass-energy flux to bring it in line with nonrelativistic
    # calculations and what the MC code actually tracks
    F_energy_upstream -= γ₀ * u₀ * ρ₀ * c^2

    return F_energy_upstream
end

"""
    upstream_machs(β₀, species, B₀)

Calculates the sonic and Alfvén mach numbers for the shock.

- For speed of sound, uses Equation (13) of Fujimura & Kennel (1979) [1979A%26A....79..299F]
- For Alfvén wave speed, uses Equation (46) of Gedalin (1993) [1993PhRvE..47.4354G]

### Arguments

- `β₀`: Far upstream plasma speed
- `species`: Vector of species in the plasma (of type `Species`)
- `B₀`: Far upstream magnetic field

### Returns
- `mach_sonic`: sonic Mach number of upstream flow
- `mach_alfven`: Alfvénic Mach number of upstream flow
"""
function upstream_machs(β₀, species, B₀)

    # Assume cold upstream plasma, so that the adiabatic index is 5/3 identically   #assumecold
    Γ = 5//3
    n₀ = number_density.(species)
    P₀ = dot(n₀, temperature.(species)) * kB  # pressure
    ρ₀ = dot(n₀, mass.(species))              # mass density

    relativistic = (β₀ ≥ β_rel_fl)

    return (
        mach_sonic_func(β₀ * c, P₀, ρ₀, Γ, relativistic),
        mach_alfven_func(β₀ * c, P₀, ρ₀, Γ, B₀, relativistic),
    )
end

# TODO name these functions better. the `_func` suffix is because there are
# globals with the same name
"""
    mach_sonic_func(u, P, ρ, Γ, relativistic)

Return sonic mach number (ratio of speed to local speed of sound)
"""
function mach_sonic_func(u, P, ρ, Γ, relativistic::Bool)
    method = relativistic ? Val(:relativistic) : Val(:classical)
    cₛ = sound_speed(method, P, ρ, Γ)
    return u / cₛ |> NoUnits
end

"""
    mach_alfven_func(u, P, ρ, Γ, B, relativistic)

Return Alfvénic mach number (ratio of speed to Alfvén wave group velocity)
"""
function mach_alfven_func(u, P, ρ, Γ, B, relativistic::Bool)
    v_A = relativistic ?
        alfven_speed(Val(:relativistic), ρ, B, P, Γ) :
        alfven_speed(Val(:classical), ρ, B)
    return u / v_A |> NoUnits
end

"""
    sound_speed(regime, P, ρ, Γ)

Return the speed of sound in a plasma.

### Arguments
- `regime`: Must be one of `Val(:classical)` or `Val(:relativistic)`
- `P`: plasma pressure
- `ρ`: plasma mass density
- `Γ`: plasma adiabatic index
"""
function sound_speed end
function sound_speed(::Val{:relativistic}, P, ρ, Γ)
    # Find FK1979's R factor, the ratio of pressure to rest energy density
    R = P / (ρ * c^2) |> NoUnits # dimensionless auxiliary variable

    # Compute the speed of sound using Fujimura & Kennel
    #     cₛ²/c² = ΓR/(aR + 1)                              FK1979 Eq. 13
    a = Γ / (Γ - 1)     # defined near FK1979 Equation (6)
    return c * √(Γ * R / (a * R + 1))
end
sound_speed(::Val{:classical}, P, ρ, Γ) = √(Γ * P / ρ) # cₛ = √(K/ρ), where K = ΓP is the bulk modulus

"""
    alfven_speed(Val(:relativistic), ρ, B, P, Γ)

Calculate the Alfvén speed for a relativistic plasma with pressure `P`,
density `ρ`, Lorentz factor `Γ`, and ambient magnetic field `B`.

It uses the Equation (46) given by Gedalin (1993), where we assume an
equation of state ``e = ρc^2 + P/(Γ-1)``.
```math
\\begin{align*}
    v_A &= c \\sqrt{\\frac{u_B}{ε + p + u_B}} \\\\
        &= \\frac{c}{\\sqrt{1 + \\frac{w}{u_B}}}
\\end{align*}
```
where ``u_B = B^2/4π`` is the magnetic field energy density,
``w = \\frac{Γ}{Γ-1} P + ρ c^2`` is the enthalpy density, ``ε`` is ??,
``p`` is ??, ``e`` is the (total internal) energy density,
``ρ c^2`` is the rest energy density, and
``P/(Γ-1)`` is the thermal component (internal kinetic energy).
"""
function alfven_speed(::Val{:relativistic}, ρ, B, P, Γ)
    enthalpy = Γ/(Γ-1) * P + ρ * c^2
    v_A = c / √(1 + 4π * enthalpy / B^2)
    return v_A
end
"""
    alfven_speed(Val(:classical), ρ, B)

Alfvén wave group velocity.
"""
alfven_speed(::Val{:classical}, ρ, B) = B / √(4π * ρ)

"""
    setup_profile(...)

Sets the initial values of the shock profile

### Arguments

- `u₀`, `β₀`, `γ₀`
- `B₀`
- `θ_B₀`
- `r_comp`
- `bturb_comp_frac`
- `bfield_amp`
- `use_custom_εB,`
- `n_ions`
- `species`
- `F_px_upstream`
- `F_energy_upstream`
- `grid_axis`
- `x_grid_cm`
- `x_grid_rg`

### Returns
- `uₓ_sk_grid`: bulk fluid velocity along x axis (i.e., perpendicular to shock face) in shock frame
- `uz_sk_grid`: bulk fluid velocity along z axis (i.e., parallel to shock face) in shock frame
- `utot_grid`: total bulk fluid velocity in shock frame
- `γ_sf_grid`: bulk flow Lorentz factor in shock frame
- `β_ef_grid`: relative x-axis speed between plasma and explosion frames
- `γ_ef_grid`: Lorentz factor associated with β_ef_grid
- `btot_grid`: total magnetic field
- `θ_grid`: angle of magnetic field[radians] relative to shock normal (i.e., to x axis)
- `εB_grid`: user-defined function for fraction of energy density in magnetic field.
  Sets value of btot_grid
- `B₂`: field strength in downstream region. Initially set in calc_downstream, it may be reset here
  depending on values of bturb_comp_frac & bfield_amp
"""
function setup_profile(
        u₀, β₀, γ₀, B₀, θ_B₀,
        r_comp, bturb_comp_frac, bfield_amp, use_custom_εB,
        n_ions, species, F_px_upstream, F_energy_upstream,
        grid_axis, x_grid_cm, x_grid_rg,
    )

    uₓ_sk_grid = OffsetVector{typeof(u₀)}(undef, grid_axis)
    uz_sk_grid = zeros(typeof(u₀), grid_axis)
    γ_sf_grid = OffsetVector{Float64}(undef, grid_axis)
    β_ef_grid = OffsetVector{Float64}(undef, grid_axis)
    γ_ef_grid = OffsetVector{Float64}(undef, grid_axis)
    btot_grid = OffsetVector{BFieldCGS}(undef, grid_axis)
    θ_grid = fill(deg2rad(θ_B₀), grid_axis)

    comp_fac = 0.0
    for i in grid_axis
        if x_grid_cm[i] < 0cm
            uₓ_sk_grid[i] = u₀
            γ_sf_grid[i] = γ₀
            β_ef_grid[i] = 0.0
            γ_ef_grid[i] = 1.0
            btot_grid[i] = B₀
        else
            u = u₀ / r_comp
            β = u / c |> NoUnits
            uₓ_sk_grid[i] = u
            γ_sf_grid[i] = 1 / √(1 - β^2)
            β_ef_grid[i] = (β₀ - β) / (1 - β₀ * β)
            γ_ef_grid[i] = 1 / √(1 - β_ef_grid[i]^2)

            # When initializing magnetic field, include necessary corrections for turbulence compression
            z_comp = (γ₀ * u₀) / (γ_sf_grid[i] * u)
            aux_fac = √((1 + 2 * z_comp^2) / 3)
            local comp_fac = 1 + (aux_fac - 1) * bturb_comp_frac
            # Also include any additional amplification specified
            amp_fac = 1 + (comp_fac - 1) * bfield_amp
            btot_grid[i] = B₀ * amp_fac
        end
    end

    utot_grid = copy(uₓ_sk_grid) # uz_sk_grid is 0, so don't bother adding
    #utot_grid = uₓ_sk_grid + uz_sk_grid
    #utot_grid = hypot.(uₓ_sk_grid, uz_sk_grid)

    εB_grid = OffsetVector{Float64}(undef, grid_axis)
    # If directed in data_input, use a custom-defined ε_B to set btot_grid.
    #-------------------------------------------------------------------------------
    if use_custom_εB
        set_custom_εB!(
            εB_grid, grid_axis,
            n_ions, species, B₀,
            F_px_upstream, F_energy_upstream, uₓ_sk_grid, x_grid_rg,
            comp_fac,
            γ₀, β₀, u₀
        )
        n₀ = dot(density.(species), mass.(species)) / mp # total number density
        e₀ = n₀ * mp * c^2  # rest energy density
        for i in grid_axis
            energy_density = (F_energy_upstream + γ₀ * u₀ * e₀) / uₓ_sk_grid[i] - F_px_upstream
            # FIXME this tries to be a square root of a negative number sometimes
            #btot_grid[i] = √(8π * εB_grid[i] * energy_density)
            btot_grid[i] = √abs(8π * εB_grid[i] * energy_density)
        end
    else
        fill!(εB_grid, 1.0e-99)
    end
    #-------------------------------------------------------------------------
    # Custom εB_grid defined if needed

    B₂ = btot_grid[end]

    return (
        uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
        β_ef_grid, γ_ef_grid, btot_grid, θ_grid, εB_grid, B₂,
    )
end

"""
    set_custom_εB!(...)

### Arguments
- `εB_grid`
- `grid_axis`
- `n_ions`
- `species`
- `B₀`
- `F_px_upstream`
- `F_energy_upstream`
- `uₓ_sk_grid`
- `x_grid_rg`
- `comp_fac`
- `γ₀`, `β₀`, `u₀`
"""
function set_custom_εB!(
        εB_grid,
        grid_axis,
        n_ions, species, B₀,
        F_px_upstream, F_energy_upstream, uₓ_sk_grid, x_grid_rg,
        comp_fac,
        γ₀, β₀, u₀
    )

    #@debug(
    #    "Input parameters", εB_grid, grid_axis, n_ions, species, B₀,
    #    F_px_upstream, F_energy_upstream, uₓ_sk_grid, x_grid_rg, comp_fac, u₀, γ₀, β₀
    #)

    # Calculate ε_B₀ which depends on far upstream magnetic field and mass density.
    # If electrons aren't a separate species, they don't contribute enough mass to be important.
    n₀ = dot(density.(species), mass.(species)) / mp # total number density
    e₀ = n₀ * mp * c^2   # upstream rest energy density
    εB₀ = B₀^2 / (8π * e₀) |> NoUnits

    # The Monte Carlo length is rg₀ = γ₀ ⋅ β₀ ⋅ mₚc² / (q ⋅ B₀). The plasma skin
    # depth is λ_SD = γ₀ ⋅ mₚc² / (4π ⋅ q² ⋅ den₀), where den₀ refers to the upstream
    # number density of electrons. With the definition
    #     σ = 2εB₀/γ₀ = B₀² / (4π ⋅ γ₀ ⋅ n₀ ⋅ mₚc²),
    # where n₀ here refers to the number density of *protons*, one can show that in
    # the shock frame (where grid exists),
    #     λ_SD = β₀ / √(σ ⋅ density_p/density_e) ⋅ rg₀.
    n₀_electron = density(species[end]) # electron number density
    σ = 2εB₀ / γ₀
    rg2sd = β₀ / √(σ * n₀ / n₀_electron) |> NoUnits

    # Also need the final value of ε_B downstream, in case our downstream region is long enough
    # that the magnetic field can decay to this value. Note that the R-H relations can be
    # rearranged to read
    #     energy_density(x) = F_en₀/u(x) - F_px₀
    # assuming flux conservation everywhere.
    energy_density₂ = (F_energy_upstream + γ₀ * u₀ * e₀) / uₓ_sk_grid[end] - F_px_upstream
    εB₂ = (B₀ * comp_fac)^2 / (8π * energy_density₂) |> NoUnits
    # Use this value to compute the distance downstream at which the field will have decayed to it.
    # Per the Blandford-McKee solution, energy ∝ 1/χ ∝ 1/distance downstream. Since we do not actually
    # modify our pressures and densities according to the BM solution, instead modify εB
    end_decay_rg = (5.0e-3 / εB₂) / rg2sd |> NoUnits

    @debug("Setting custom εB", n₀, εB₀, n₀_electron, σ, rg2sd, energy_density₂, εB₂, end_decay_rg)

    # Now we can calculate εB_grid. Per the Blandford-McKee solution,
    # energy ∝ 1/χ ∝ 1/distance downstream. Since we do not actually modify our
    # pressures and densities according to the BM solution, instead modify ε_B
    for i in grid_axis
        x_grid_sd = x_grid_rg[i] * rg2sd
        if x_grid_sd < -50
            εB_grid[i] = max(1.04e-5 / abs(x_grid_sd)^0.6, εB₀)
        elseif x_grid_sd < 50
            εB_grid[i] = 1.0e-4
        elseif x_grid_rg[i] < end_decay_rg
            εB_grid[i] = 5.0e-3 / x_grid_sd
        else
            εB_grid[i] = εB₂
        end
    end

    return εB_grid
end

"""
    init_pop(...)

Initializes the particle populations that will propagate through the shock structure.
Handles fast push and associated flux-tracking & changes to the population

### Arguments

- `do_fast_push`
- `inp_distr`
- `i_ion`
- `m`
- `rng`
- `T₀_ion`
- `energy_inj`
- `inj_weight`
- `n_pts_inj`
- `n₀_ion`
- `x_grid_start`
- `rg₀`
- `η_mfp`
- `x_fast_stop_rg`
- `β₀`, `γ₀`, `u₀`
- `n_ions`
- `m_ion`
- `n_grid`
- `x_grid_rg`
- `uₓ_sk_grid`
- `γ_sf_grid`
- `ptot_inj`
- `weight_inj`
- `n_pts_MB`

### Returns

- `n_pts_use`
- `i_grid_in`
- `weight_in`
- `ptot_pf_in`
- `pb_pf_in`
- `x_PT_cm_in`
- `pxx_flux`
- `pxz_flux`
- `energy_flux`
"""
function init_pop(
        do_fast_push, inp_distr, i_ion, m, rng,
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

        ptot_inj[:, i_ion], weight_inj[:, i_ion], n_pts_MB[i_ion] = set_inj_dist(
            inj_weight, n_pts_inj, inp_distr, T_or_E, m, n₀_ion[i_ion]
        )

        n_pts_use = n_pts_MB[i_ion]
        weight_in = weight_inj[1:n_pts_use, i_ion]
        ptot_pf_in = ptot_inj[1:n_pts_use, i_ion]
        pb_pf_in = ptot_pf_in[1:n_pts_use] .* 2 * (rand(rng, n_pts_use) .- 0.5)
        x_PT_cm_in = fill(x_grid_start - 10 * rg₀ * η_mfp, n_pts_use)
        pxx_flux = zeros(MomentumDensityFluxCGS, n_grid)
        pxz_flux = zeros(MomentumDensityFluxCGS, n_grid)
        energy_flux = zeros(EnergyDensityFluxCGS, n_grid)

        i_grid_in = zeros(Int, n_pts_use)
        return (
            n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in,
            x_PT_cm_in, pxx_flux, pxz_flux, energy_flux,
        )
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
    # gas, to calculate the far upstream internal energy   #assumecold
    Γ_sph = 5//3

    temp_ratio = density_ratio^Γ_sph / density_ratio

    if (kB * T₀_ion[i_ion] * temp_ratio) > (4 * m * c^2 * E_rel_pt)
        error(
            "Fast push cannot work because highest energy thermal particles become mildly relativistic. ",
            "Move fast push location upstream or disable entirely."
        )
    end
    pxx_flux = zeros(MomentumDensityFluxCGS, n_grid)
    pxz_flux = zeros(MomentumDensityFluxCGS, n_grid)
    energy_flux = zeros(EnergyDensityFluxCGS, n_grid)


    # Only run through the flux updates for the first particle species (i.e.
    # protons, not that it matters here); skip thereafter
    if i_ion == 1
        F_update!(
            pxx_flux, pxz_flux, energy_flux,
            m_ion, n₀_ion, T₀_ion, relativistic,
            i_stop,
            u₀, γ₀, γ_sf_grid, uₓ_sk_grid,
        )
    end  # check on i_ion


    # With fast push fluxes taken care of, create the particle distribution
    # that will be injected at x_fast_stop
    T_or_E = T₀_ion[i_ion] * temp_ratio

    ptot_inj[:, i_ion], weight_inj[:, i_ion], n_pts_MB[i_ion] = set_inj_dist(
        inj_weight, n_pts_inj, inp_distr, T_or_E, m, n₀_ion[i_ion]
    )
    n_pts_use = n_pts_MB[i_ion]

    weight_in = weight_inj[1:n_pts_use, i_ion]
    ptot_pf_in = ptot_inj[1:n_pts_use, i_ion]
    x_PT_cm_in = fill(x_fast_stop_rg * rg₀, n_pts_use)
    i_grid_in = fill(i_stop, n_pts_use)

    u = uₓ_sk_grid[i_stop]
    βᵤ = u / c |> NoUnits

    pb_pf_in = zeros(MomentumCGS, n_pts_use)
    for i_prt in 1:n_pts_use

        # Per Vladimirov+ (2009) [PhD], particle velocities should not be isotropic
        # in plasma frame. Must be weighted according to shock-frame pitch angle.
        #------------------------------------------------------------------------

        if relativistic
            γ_pf = hypot(1, ptot_pf_in[i_prt] / (m * c))
            β_pf = sqrt(1 - 1 / γ_pf^2)

            # We want v² to have a uniform distribution, so v must have a
            # triangular/linear distribution (with parameter b=c, i.e., triangle
            # peaks at right vertex)
            # Unitful distributions aren't yet supported, so do calculations as
            # β. Also saves having to divide/multiply a lot by c.
            βmin = abs((βᵤ - β_pf) / (1 - βᵤ * β_pf))
            βmax = abs((βᵤ + β_pf) / (1 + βᵤ * β_pf))
            dist_β_sf = TriangularDist(βmin, βmax, βmax)

            βx_sf = rand(rng, dist_β_sf)

            vx_pf = (βx_sf - βᵤ) / (1 - βx_sf * βᵤ) * c
        else
            γ_pf = 1.0
            vt_pf = ptot_pf_in[i_prt] / m

            # We want v² to have a uniform distribution, so v must have a
            # triangular/linear distribution (with parameter b=c, i.e., triangle
            # peaks at right vertex)
            # Unitful distributions not yet supported :(
            vmin = ustrip(cm / s, abs(u - vt_pf))
            vmax = ustrip(cm / s, abs(u + vt_pf))
            dist_v_sf = TriangularDist(vmin, vmax, vmax)

            vx_sf = cm / s * rand(rng, dist_v_sf)

            vx_pf = vx_sf - u
        end

        pb_pf_in[i_prt] = γ_pf * m * vx_pf
        #------------------------------------------------------------------------
        # Velocity-weighted pitch angles finished
    end

    return n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux
end


"""
    F_update!(...)

Update the flux arrays as if the particles had actually crossed them.

### Arguments
- `pxx_flux`: Modified in-place.
- `pxz_flux`: Modified in-place.
- `energy_flux`: Modified in-place.
- `m_ion`
- `n₀_ion`
- `T₀_ion`
- `relativistic`
- `i_stop`
- `u₀`
- `γ₀`
- `γ_sf_grid`
- `uₓ_sk_grid`
"""
function F_update!(
        pxx_flux, pxz_flux, energy_flux,
        m_ion, n₀_ion, T₀_ion, relativistic::Bool,
        i_stop,
        u₀, γ₀, γ_sf_grid, uₓ_sk_grid
    )
    P₀ = dot(n₀_ion, T₀_ion) * kB       # Upstream thermal pressure
    ρ₀ = dot(n₀_ion, m_ion)             # Upstream mass density

    # Assume an adiabatic index of 5/3, appropriate for non-relativistic ideal gas,
    # to calculate the far upstream pressure and internal energy   #assumecold
    Γ_sph = 5//3
    Ξ_sph = Γ_sph / (Γ_sph - 1)

    # Calculate fluxes and update arrays; note that if fast push isn't enabled
    # i_stop = 0, and this loop never executes
    for i in 1:i_stop

        u_curr = uₓ_sk_grid[i]              # current speed
        β_curr = u_curr / c                 # current speed (in units of c)
        γ_curr = γ_sf_grid[i]
        γβ_curr = γ_curr * β_curr

        density_ratio = (γ₀ * u₀) / (γ_curr * u_curr)
        ρ_curr = ρ₀ * density_ratio         # current mass density
        # Note assumption that Γ_sph doesn't change from zone to zone: #assumecold
        P_curr = P₀ * density_ratio^Γ_sph   # current pressure

        # Determine fluxes while handling different possible orientations and shock speeds.
        # For non-relativistic fluxes, expand out to β^2 to allow for better
        # matching with relativistic versions.
        # WARNING: these fluxes do not include contributions from a strong
        # magnetic field. This is incorporated during the smoothing process.
        #----------------------------------------------------------------------
        F_pz = 0.0erg / cm^3
        if !relativistic
            F_pₓ = (
                ρ_curr * u_curr^2 * (1 + β_curr^2)
                    + P_curr * (1 + Ξ_sph * β_curr^2)
            )
            F_energy = (
                ρ_curr / 2 * u_curr^3 * (1 + 1.25 * β_curr^2)
                    + P_curr * u_curr * Ξ_sph * (1 + β_curr^2)
            )
        else
            e_curr = ρ_curr * c^2 # energy density

            F_pₓ = P_curr + γβ_curr^2 * (e_curr + Ξ_sph * P_curr)
            F_energy = (
                γβ_curr * γ_curr * c * (e_curr + Ξ_sph * P_curr)
                    # Subtract mass-energy flux from F_energy to bring
                    # it in line with non-relativistic calculations
                    - γβ_curr * c * e_curr
            )
        end
        #--------------------------------------------------------------------
        # Fluxes calculated

        pxx_flux[i] = F_pₓ
        pxz_flux[i] = F_pz
        energy_flux[i] = F_energy

    end  # loop over grid location
    #------------------------------------------------------------------------
    # Arrays updated through i_fast_stop
    return
end

"""
    set_inj_dist(inj_weight, n_pts_inj, inp_distr, T_or_E, m, n₀)

Sets the injected particle distributions for all species. Initially, particles are placed in
a Maxwell-Boltzmann (thermal) distribution based on supplied temperature and particle mass.
This is corrected at the end if a δ-function distribution was requested (not worried about
the wasted computation because this subroutine runs only rarely).

### Arguments

- `inj_weight`: whether each particle or each bin has equal weights
  (`true` for equal weight particles, `false` for equal weight bins)
- `n_pts_inj`: target # of particles for distribution
- `inp_distr`: thermal, δ-function, or some other distribution
- `T_or_E`: if using thermal distribution, this is temperature[K]; if δ function, it's injection energy[keV]
- `m`: mass for this particle species
- `n₀`: far upstream number density for this species

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
        # Can't say anything about total number of particles in this case,
        # because haven't split them into M-B distribution yet
        n_per_bin = -1
        n_pts_tot = -1
    else
        n_per_bin = n_pts_inj ÷ num_therm_bins  # Integer math loses excess particles, but that's intended behavior
        n_pts_tot = n_per_bin * num_therm_bins

        if n_per_bin < 5
            throw(ArgumentError("too few particles per bin ($n_per_bin; need at least 5). Increase n_pts_inj."))
        end
    end
    #------------------------------------------------------------------------
    # End administrative section

    p_range = create_inj_dist_momentum_range(m, T_or_E, num_therm_bins)
    Δp = step(p_range)

    E₀ = m * c^2
    kT = kB * T_or_E  # Working under assumption of thermal dist now
    # define energy over kT
    if (kT / E₀) < E_rel_pt # Does thermal energy go over relativistic cutoff?
        E_range = @. p_range ^ 2 / (2m * kT)
    else
        E_range = @. hypot(p_range * c, E₀) / kT
    end

    # Generate the Maxwell-Boltzmann distribution. The actual calculation of f1 and f2 does
    # not need modification for relativistic momenta. Such a modification would only affect
    # the normalization of the curve, not the dependence on momentum for a particular value
    # of kT. Since we only care about the relative fractional area of each bin, overall
    # normalization doesn't matter.
    #------------------------------------------------------------------------
    # Find total area under M-B curve
    area_tot = calc_MB_area(p_range, E_range)

    ptot_out = zeros(MomentumCGS, na_particles)
    weight_out = zeros(na_particles)

    # Fill the bins with particles. Needs to be done in two separate loops, despite the
    # similarities, since particle weights are handled differently based on value of inj_weight
    if inj_weight # First, particles have equal weight
        n_pts_tot = set_inj_dist_particle_equal_weight!(
            ptot_out, weight_out, p_range, E_range, area_tot, n_pts_inj, Δp, n₀
        )
    else # Bins have equal weight
        set_inj_dist_bin_equal_weight!(
            ptot_out, weight_out, p_range, E_range, area_tot, Δp, n₀, n_per_bin
        )
    end  # test of inj_weight

    n_pts_use = n_pts_tot
    #------------------------------------------------------------------------
    # End of Maxwell-Boltzmann section


    # If the particles are to use a δ-function distribution, set it here
    #--------------------------------------------------------------------------
    if inp_distr == 2
        n_pts_use = n_pts_inj
        set_δ_distr_inj_weights!(
            ptot_out, weight_out, 1:n_pts_inj, n_pts_tot, m, n₀, T_or_E, E_rel_pt,
        )
    end  # check of inp_distr
    #--------------------------------------------------------------------------
    # δ-function handled


    return (ptot_out, weight_out, n_pts_use)
end

"""
    calc_MB_area(p_range, E_range)

Area under the curve for a Maxwell-Boltzmann distribution,
up to arbitrary normalization.

### Arguments
- `p_range`
- `E_range`

### Returns
- `area_tot`
"""
function calc_MB_area(p_range::AbstractRange, E_range::AbstractVector)
    area_tot = 0.0g*cm/s
    Δp = step(p_range)
    for (i, p1) in enumerate(p_range[begin:(end - 1)])
        area_tot += calc_MB_area_single_bin(p1, p_range[i + 1], E_range[i], E_range[i + 1])
    end
    return area_tot
end

"""
    calc_MB_area_single_bin(p1, p2, E1, E2)

Area under the curve for a single bin of Maxwell-Boltzmann distribution,
up to arbitrary normalization. Uses trapezoid rule for integration.

### Arguments
- `p1`, `p2`: Momentum bounds of bin
- `E1`, `E2`: Energy bounds of bin (assumes rel/nonrel calculations are already
  done), in units of kT

### Returns

The area occupied on the bin under a Maxwell-Boltzmann distribution of
temperature T (the temperature is implicit in `E1` and `E2`).
"""
function calc_MB_area_single_bin(p1, p2, E1, E2)
    # Start working in log space because of potentially huge exponents
    log_f1 = 2 * log(p1 / (g*cm/s)) - E1
    log_f2 = 2 * log(p2 / (g*cm/s)) - E2
    f1 = exp(log_f1)
    f2 = exp(log_f2)
    # Integrate using the trapezoid rule
    return (p2 - p1) * (f1 + f2) / 2
end

"""
Create range of momenta over which the Maxwell-Boltzmann particle
distribution is set.

### Arguments
- `T`: distribution temperature[K]
- `m`: mass for this particle species

### Returns
- `p_range`: Array of momenta
"""
function create_inj_dist_momentum_range(m, T, nbins)
    E₀ = m * c^2

    kT = kB * T
    # Minimum, maximum extent of Maxwell-Boltzmann distribution
    kT_min = 2.0e-3 * kT
    kT_max = 10 * kT

    relativistic = (kT / E₀) ≥ E_rel_pt

    # Find min and max momenta of M-B curve
    if !relativistic
        # In non-relativistic case, kinetic energy ≈ thermal energy
        p_min = √(2m * kT_min)
        p_max = √(2m * kT_max)
    else
        # Once particles are relativistic, rest energy becomes important:
        #     E² = p²c² + m²c⁴  ≈  (kT + mc²)²
        p_min = √((kT_min + E₀)^2 - E₀^2) / c
        p_max = √((kT_max + E₀)^2 - E₀^2) / c
    end

    Δp = (p_max - p_min) / nbins
    p_range = range(start = p_min, step = Δp, length = nbins + 1)

    return p_range
end

function set_inj_dist_particle_equal_weight!(
        ptot_out::AbstractVector, weight_out::AbstractVector,
        p_range, E_range, area_tot, n_pts_inj::Integer, Δp, n₀
    )

    area_per_pt = area_tot / n_pts_inj # Area each particle gets if inj_weight = T
    n_pts_tot = 1  # Total number of particles in M-B distribution

    for (i, p1) in enumerate(@view p_range[begin:(end - 1)])
        p2 = p_range[i + 1]

        E1 = E_range[i]
        E2 = E_range[i + 1]

        bin_area = calc_MB_area_single_bin(p1, p2, E1, E2)

        area_frac = bin_area / area_per_pt

        # Rounded particle count for this bin
        n_pts_this_bin = round(Int, area_frac)

        slice = (n_pts_tot + 1):(n_pts_tot + n_pts_this_bin)

        # Geometric center of bin; particles in this bin will receive this momentum
        # particles in this bin will have momentum = the geometric center of the bin
        ptot_out[slice] .= √(p1 * p2)
        n_pts_tot += n_pts_this_bin
    end

    # If each particle has equal weight, then the total weight of the
    # distribution should be proportional to the density of the species
    weight_out[1:n_pts_tot] .= ustrip(cm^-3, n₀) / n_pts_tot

    return n_pts_tot
end

"""
    set_inj_dist_bin_equal_weight!(...)

### Arguments

- `ptot_out`: Modified in-place.
- `weight_out`: Modified in-place.
- `p_range`
- `E_range`
- `area_tot`
- `Δp`
- `n₀`
- `n_per_bin`

### Returns

- `ptot_out`
- `weight_out`
"""
function set_inj_dist_bin_equal_weight!(
        ptot_out, weight_out, p_range, E_range, area_tot, Δp, n₀, n_per_bin,
    )
    for (i, p1) in enumerate(@view p_range[begin:(end - 1)])
        p2 = p_range[i + 1]

        E1 = E_range[i]
        E2 = E_range[i + 1]

        bin_area = calc_MB_area_single_bin(p1, p2, E1, E2)

        area_frac = bin_area / area_tot

        # indices for start and end of bin
        jstart = (i - 1) * n_per_bin + 1
        jend = jstart + (n_per_bin - 1)

        # particles in this bin will have momentum = the geometric center of the bin
        ptot_out[jstart:jend] .= √(p1 * p2)
        weight_out[jstart:jend] .= area_frac / n_per_bin * n₀
    end  # loop over i
    return ptot_out, weight_out
end

function set_δ_distr_inj_weights!(
        ptot_out, weight_out, slice, n_pts_tot, m, n₀, E, E_rel_pt,
    )
    E₀ = m * c^2
    E_inj = ustrip(erg, E * keV)    # injection energy

    if E_inj / E₀ < E_rel_pt
        p = √(2m * E_inj)
    else
        p = √(E_inj^2 - E₀^2) / c
    end

    ptot_out[slice] .= p
    weight_out[slice] .= ustrip(cm^-3, n₀) / n_pts_tot

    return ptot_out, weight_out
end
end # module
