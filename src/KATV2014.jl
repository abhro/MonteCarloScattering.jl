"""
Implementation of code from Kafexhiu et al. (2014)

### References

Kafexhiu, E., Aharonian, F., Taylor, A. M., & Vila, G. S. (2014).
*Parametrization of gamma-ray production cross sections for pp interactions in a
broad proton energy range from the kinematic threshold to PeV energies.*
Physical Review D, 90(12), 123014. https://doi.org/10.1103/PhysRevD.90.123014
"""
module KATV2014

using Unitful: mp
using ..constants: Tₜₕ, Mᵣₑₛ, Γᵣₑₛ, E₀_π⁰

"""
    get_σ_π(Tₚ, i_data, s_ECM)

Calculates the total inclusive π-0 production cross section using method
of Kafexhiu+ (2014) [PhRvD, V.90, 123014], section 4.
"""
function get_σ_π(Tₚ, i_data, s_ECM)
    # If Tₚ < 2 GeV, calculate cross sections using experimental data and equations (2) - (5)
    # Note that we have already screened for Tₜₕ < Tₚ prior to this subroutine
    #-------------------------------------------------------------------------
    if Tₚ < 2

        eq4_γ = Mᵣₑₛ * hypot(Mᵣₑₛ, Γᵣₑₛ)
        eq4_K = √8 * Mᵣₑₛ * Γᵣₑₛ * eq4_γ / (π * √(Mᵣₑₛ^2 + eq4_γ))
        f_BW = mp * eq4_K / (((√s_ECM - mp)^2 - Mᵣₑₛ^2)^2 + Mᵣₑₛ^2 * Γᵣₑₛ^2)
        eq3_η = √((s_ECM - E₀_π⁰^2 - 4 * mp^2)^2 - (4 * E₀_π⁰ * mp)^2) / (2E₀_π⁰ * √s_ECM)

        # Equation (2)
        σ_1π = 7.66e-3 * eq3_η^1.95 * (1 + eq3_η + eq3_η^5) * f_BW^1.86

        # Double pion production (Equation (5)) can only happen if Tₚ > 2Tₜₕ
        σ_2π = Tₚ < 2Tₜₕ ? 0.0 : 5.7 / (1 + exp(-9.3 * (Tₚ - 1.4)))

        # Combine σ_1π and σ_2π
        return σ_1π + σ_2π
    end


    # If 2 GeV < Tₚ < 5 GeV, then use Equations (1) and (6)
    #-------------------------------------------------------------------------
    if Tₚ < 5

        eq6_Q = (Tₚ - Tₜₕ) / mp
        n_π0 = -6.0e-3 + 0.237 * eq6_Q - 0.023 * eq6_Q^2  # Neutral pion multiplicity
        eq1_ratio = Tₚ / Tₜₕ
        log_ratio = log(eq1_ratio)
        σ_inel = (30.7 - 0.96 * log_ratio + 0.18 * log_ratio^2) * (1 - eq1_ratio^(-1.9))^3 # Total inelastic cross section
        # Combine to get σ_π
        return n_π0 * σ_inel
    end


    # If 5 GeV < Tₚ, then use Equations (1) and (7). Parameters chosen for Equation (7)
    # depend on values of Tₚ and i_data
    #-------------------------------------------------------------------------

    if i_data == 2 && Tₚ > 5.0e1
        # Use PYTHIA 8 parametrization
        a₁ = 0.652
        a₂ = 0.0016
        a₃ = 0.488
        a₄ = 0.1928
        a₅ = 0.483
    elseif i_data == 3 && Tₚ > 1.0e2
        # Use SIBYLL 2.1 parametrization
        a₁ = 5.436
        a₂ = 0.254
        a₃ = 0.072
        a₄ = 0.075
        a₅ = 0.166
    elseif i_data == 4 && Tₚ > 1.0e2
        # Use QGSJET-I parametrization
        a₁ = 0.908
        a₂ = 0.0009
        a₃ = 6.089
        a₄ = 0.176
        a₅ = 0.448
    else
        # Use GEANT 4 parametrization because of choice of i_data or because Tₚ too low for
        # validity range of other models
        a₁ = 0.728
        a₂ = 0.596
        a₃ = 0.491
        a₄ = 0.2503
        a₅ = 0.117
    end  # check on i_data and Tₚ

    ξ = (Tₚ - 3) / mp
    # Neutral pion multiplicity
    n_π0 = a₁ * ξ^a₄ * (1 + exp(-a₂ * ξ^a₅)) * (1 - exp(-a₃ * ξ^0.25))

    eq1_ratio = Tₚ / Tₜₕ
    log_ratio = log(eq1_ratio)
    σ_inel = (30.7 - 0.96 * log_ratio + 0.18 * log_ratio^2) * (1 - eq1_ratio^(-1.9))^3 # Total inelastic cross section
    # Combine to get σ_π
    return n_π0 * σ_inel
end

const get_sigma_pi = get_σ_π

@doc raw"""
Calculates F(Tₚ, Eᵧ) as defined in Equations (9) and (11) of Kafexhiu+ (2014) [PhRvD, V.90, 123014].

Equation 9:

```math
\begin{align*}
Y_γ        &= E_γ        + \frac{m_π^2}{4E_γ} , \\
Y_γ^{\max} &= E_γ^{\max} + \frac{m_π^2}{4E_γ^{\max}} , \\
X_γ        &= \frac{Y_γ - m_π}{Y_γ^{\max} - m_π}
\end{align*}
```

Equation 11:

```math
F(T_p, E_γ) = \frac{ ( 1 - X_γ^{α(Tₚ)} )^{β(Tₚ)} }{ ( 1 + X_γ / C )^{γ(T_p)} }
```

Also uses equation 14:

```math
\begin{align*}
b_0 &= 5.9, &
κ(T_p) &= 3.29 - \frac{1}{5} θ_p^{-3/2}
\end{align}
```

and equation 15:

```math
μ(T_p) = \frac{5}{4} q^{5/4} \exp\left(- \frac{5}{4} q\right)
```
"""
function get_Ffunc(Tₚ, Eᵧ, i_data, Eᵧ_max)

    # Compute Xᵧ (independent variable to be used in Equation (11)) using Equation (9)
    Yᵧ = Eᵧ + E₀_π⁰^2 / Eᵧ
    Y_max = Eᵧ_max + E₀_π⁰^2 / Eᵧ_max
    Xᵧ = (Yᵧ - E₀_π⁰) / (Y_max - E₀_π⁰)

    # Now use Xᵧ, Tₚ, and i_data to calculate F(Tₚ, Eᵧ)
    #-------------------------------------------------------------------------
    if Xᵧ < 0 || Xᵧ > 1
        return 0.0 # Current photon energy cannot be produced at this proton energy
    end
    if Tₚ < 1
        # Use parametrization based on experimental data (Eq 14)
        θ = Tₚ / mp
        κ = 3.29 - 0.2 * θ^(-1.5)
        β = κ
        return (1 - Xᵧ)^β
    end

    if Tₚ < 4
        # Use parametrization based on low energy GEANT 4 fit (equation 15)
        q = (Tₚ - 1) / mp
        μ = 1.25 * q^1.25 * exp(-1.25 * q)
        λ = 3.0
        α = 1.0
        β = μ + 2.45
        γ = μ + 1.45
    elseif Tₚ < 20
        # Use parametrization based on medium energy GEANT 4 fit (equation 15)
        q = (Tₚ - 1) / mp
        μ = 1.25 * q^1.25 * exp(-1.25 * q)
        λ = 3.0
        α = 1.0
        β = 1.5 * μ + 4.95
        γ = μ + 1.5
    elseif i_data == 1 && Tₚ > 100
        # Use parametrization based on highest energy GEANT 4 fit
        λ = 3.0
        α = 0.5
        β = 4.9
        γ = 1.0
    elseif i_data == 2 && Tₚ > 50
        # Use PYTHIA 8 parametrization
        λ = 3.5
        α = 0.5
        β = 4.0
        γ = 1.0
    elseif i_data == 3 && Tₚ > 100
        # Use SIBYLL 2.1 parametrization
        λ = 3.55
        α = 0.5
        β = 3.6
        γ = 1.0
    elseif i_data == 4 && Tₚ > 100
        # Use QGSJET-I parametrization
        λ = 3.55
        α = 0.5
        β = 4.5
        γ = 1.0
    else
        # Use GEANT 4 parametrization at moderate energies, because of choice of
        # i_data or because Tₚ too low/high for validity range of other models
        λ = 3.0
        α = 0.5
        β = 4.2
        γ = 1.0
    end

    C = λ * E₀_π⁰ / Y_max
    return (1 - Xᵧ^α)^β / (1 + Xᵧ / C)^γ
end

@doc raw"""
Calculates Amax as defined in Equation (12) of Kafexhiu+ (2014) [PhRvD, V.90, 123014].

```math
A_\max(Tₚ) = \begin{cases}
    b₀ ⋅ \frac{σ_π(Tₚ)}{E_π^\max}        & \text{ for } Tₚ^\mathrm{th} ≤ Tₚ < 1\,\mathrm{GeV} \\
    b₁ θₚ^{-b₂} \exp(b₃ \log^2 (θₚ)) ⋅ \frac{σ_π(Tₚ)}{mₚ} & \text{ for } Tₚ ≥ 1\,\mathrm{GeV}
\end{cases}
```
"""
function get_Amax(Tₚ, i_data, s_ECM, σ_π)

    # Determine Emax_π_LAB, the maximum pion energy (in the lab frame,
    # where the target nucleon is at rest) allowed by kinematics
    #----------------------------------------------------------------------
    # Maximum π-0 energy in center-of-mass frame
    E_π_CM = (s_ECM - 4mp^2 + E₀_π⁰^2) / (2 * √s_ECM)
    # Lorentz factor of proton in center-of-mass frame; associated β
    γ_CM = (Tₚ + 2mp) / √s_ECM
    β_CM = √(1 - 1 / γ_CM^2)
    # Maximum π-0 momentum in center-of-mass frame. Note the unusual units
    # of GeV, which will allow easier calculation of Emax_π_LAB
    P_π_CM = √(E_π_CM^2 - E₀_π⁰^2)

    Emax_π_LAB = γ_CM * (E_π_CM + P_π_CM * β_CM)
    #-------------------------------------------------------------------------
    # Emax_π_LAB calculated


    # Use Emax_π_LAB to calculate maximum photon energy in the laboratory frame, Eᵧ_max
    #-------------------------------------------------------------------------
    # Maximum π-0 Lorentz factor in laboratory frame; associated β
    γ_LAB = Emax_π_LAB / E₀_π⁰
    β_LAB = √(1 - 1 / γ_LAB^2)

    Eᵧ_max = E₀_π⁰ / 2 * γ_LAB * (1 + β_LAB)
    #-------------------------------------------------------------------------
    # Eᵧ_max found


    # Now calculate Amax, which depends on values of i_data and Tₚ;
    # we have already screened for Tₜₕ < Tₚ. # Use Equation (12) and Table VII.
    #-------------------------------------------------------------------------
    if Tₚ < 1
        # Based off of experimental data
        Amax = 5.9 * σ_π / Emax_π_LAB

    else
        if i_data == 1 && Tₚ < 5
            # Use GEANT 4 parametrization for low energies
            b₁ = 9.53
            b₂ = 0.52
            b₃ = 0.054
        elseif i_data == 2 && Tₚ > 50
            # Use PYTHIA 8 parametrization
            b₁ = 9.06
            b₂ = 0.3795
            b₃ = 0.01105
        elseif i_data == 3 && Tₚ > 100
            # Use SIBYLL 2.1 parametrization
            b₁ = 10.77
            b₂ = 0.412
            b₃ = 0.01264
        elseif i_data == 4 && Tₚ > 100
            # Use QGSJET-I parametrization
            b₁ = 13.16
            b₂ = 0.4419
            b₃ = 0.01439
        else
            # Use GEANT 4 parametrization at higher energies, because of choice of
            # i_data or because Tₚ too low for validity range of other models
            b₁ = 9.13
            b₂ = 0.35
            b₃ = 0.0097
        end
        θ = Tₚ / mp
        Amax = b₁ * θ^(-b₂) * σ_π / mp * exp(b₃ * log(θ)^2)
    end

    #-------------------------------------------------------------------------
    # Amax determined

    return Eᵧ_max, Amax
end
end # module
