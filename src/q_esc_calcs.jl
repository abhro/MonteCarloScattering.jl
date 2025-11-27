using Unitful: k as kB, c
using LinearAlgebra: dot
include("parameters.jl")
using .parameters: β_rel_fl

"""
    q_esc_calcs(...)

Use the Rankine-Hugoniot relations to calculate the escaping momentum & energy flux.
"""
function q_esc_calcs(
        Γ, r_comp, r_RH, u₀, β₀, γ₀, species,
        γ₂, β₂, u₂)

    # Quick test first. If r_comp = r_RH, we expect no escaping flux.
    r_comp == r_RH && return 0.0, 0.0

    #--------------------------------------------------------------------------
    #  Two possibilities for R-H relations: relativistic/nonrelativistic.
    #  Determine which of the four to use. Cutoff for relativistic/nonrelativistic is set in
    #  module 'controls'
    #--------------------------------------------------------------------------
    relativistic = (β₀ ≥ β_rel_fl)

    Γ_fac = Γ / (Γ - 1)
    m_ion = mass.(species)

    T₀_ion = temperature.(species)
    n₀_ion = density.(species)
    if relativistic
        return q_esc_calcs_relativistic(n₀_ion, T₀_ion, m_ion, (u₀, β₀, γ₀), (u₂, β₂, γ₂), Γ_fac)
    else
        return q_esc_calcs_nonrelativistic(n₀_ion, T₀_ion, m_ion, (u₀, β₀, γ₀), (u₂, β₂, γ₂), Γ_fac)
    end
end


"""
    q_esc_calcs_nonrelativistic(...)

Solution comes from Ellison (1985) [1985JGR....90...29E] (Eqs 8-10).
Note assumption of zero escaping momentum flux, which is good to
within a couple percent for strong nonrelativistic shocks.
#TODO: check how much of a difference this assumption makes
"""
function q_esc_calcs_nonrelativistic(n₀_ion, T₀_ion, m_ion, (u₀, β₀, γ₀), (u₂, β₂, γ₂), Γ_fac)
    # Calculate thermal pressure of far upstream gas
    P₀ = dot(n₀_ion, T₀_ion) * kB   # pressure (thermal)
    ρ₀ = dot(n₀_ion, m_ion)         # mass density

    # Calculate upstream incoming energy flux   #assumecold
    F_pₓ_upstream_fl     = ρ₀ * u₀^2 + P₀
    F_energy_upstream_fl = ρ₀ * u₀^3 / 2 + 5//2 * P₀ * u₀

    # Calculate far downstream density (Eq 8) and pressure (Eq 9)
    ρ₂ = ρ₀ * γ₀*β₀ / (γ₂*β₂)           # mass density
    P₂ = F_pₓ_upstream_fl - ρ₂*u₂^2     # pressure

    # Calculate escaping energy flux using Eq (10)
    Q_en = F_energy_upstream_fl - ρ₀ * u₀ * u₂^2 / 2 - P₂ * u₂ * Γ_fac

    # Finally, put in units of F_en₀
    q_esc_cal_energy = Q_en / F_energy_upstream_fl
    q_esc_cal_pₓ = 0.0

    return q_esc_cal_energy, q_esc_cal_pₓ
end

@doc raw"""
    q_esc_calcs_relativistic(...)

Solution comes from Ellison+ (1990) [1991ApJ...378..214E]. Uses relativistic Rankine-Hugoniot
relations. See that paper for details of equations and associated quantities. Briefly,

```math
\begin{align}
g_0   n_0 b_0         &=  g_2   n_2 b_2                  \tag{RH1} \\
g_0^2 w_0 b_0^2 + P_0 &=  g_2^2 w_2 b_2^2 + P_2 + Q_{px} \tag{RH2} \\
g_0^2 w_0 b_0 c       &=  g_2^2 w_2 b_2 c       + Q_{en} \tag{RH3}
\end{align}
```
where
- ``w   = E_0 + E_k + P``   is the enthalpy as total energy density + pressure
- ``E_0 =       nmc^2``     is the rest mass energy density
- ``E_k = (γ-1) nmc^2``     is the kinetic energy density, with ``γ = \sqrt{1 + (p/mc)^2}``
- ``P   = \frac{1}{3} npv`` is the pressure

For closure, it is assumed that the two escaping fluxes are related by
```math
Q_{en} = \sqrt{(1+β₀)/2} ⋅ Q_{px} ⋅ c,
```
i.e. the geometric mean of the arithmetic mean of u₀ and c. This allows the solution to
smoothly join with the non-relativistic version. Use only fluid component of fluxes,
not fluid+EM, for now.
"""
function q_esc_calcs_relativistic(n₀_ion, T₀_ion, m_ion, (u₀, β₀, γ₀), (u₂, β₂, γ₂), Γ_fac)
    # Factor relating Q_en and Q_px
    q_fac = c * √((1 + β₀)/2)

    # Calculate thermal pressure of far upstream gas
    P₀ = dot(n₀_ion, T₀_ion) * kB   # pressure (thermal)
    ρ₀ = dot(n₀_ion, m_ion)         # mass density

    # Two terms to simplify the calculation of pressure₂.   #assumecold
    F_pₓ_upstream_fl     = γ₀^2 * β₀^2 * (ρ₀*c^2 + 5//2*P₀) + P₀
    F_energy_upstream_fl = γ₀^2 * u₀   * (ρ₀*c^2 + 5//2*P₀)
    term_aux = γ₂^2 * (q_fac * β₂^2 - u₂)

    # Calculate far downstream density and pressure
    ρ₂ = ρ₀ * γ₀*β₀ / (γ₂*β₂) # mass density
    P₂ = (q_fac * F_pₓ_upstream_fl - F_energy_upstream_fl - term_aux*ρ₂*c^2) / (q_fac + Γ_fac*term_aux)

    # Calculate Q_px & Q_en (the physical escaping momentum & energy fluxes)
    Q_px = F_pₓ_upstream_fl - (γ₂*β₂)^2 * (ρ₂ * c^2 + Γ_fac * P₂) - P₂
    Q_en = Q_px * q_fac

    # Convert Q_en & Q_px into q_esc_cal_**, which involves dividing by the far upstream values.
    # Subtract off mass-energy flux from F_energy_upstream_fl to bring results in line with
    # non-relativistic calculation.
    q_esc_cal_energy = Q_en / (F_energy_upstream_fl - γ₀ * u₀ * ρ₀*c^2)
    q_esc_cal_pₓ = Q_px / F_pₓ_upstream_fl

    return q_esc_cal_energy, q_esc_cal_pₓ
end
