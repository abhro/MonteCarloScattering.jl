"""
Implementation of code from Kafexhiu et al. (2014)

### References

Kafexhiu, E., Aharonian, F., Taylor, A. M., & Vila, G. S. (2014).
*Parametrization of gamma-ray production cross sections for pp interactions in a
broad proton energy range from the kinematic threshold to PeV energies.*
Physical Review D, 90(12), 123014. https://doi.org/10.1103/PhysRevD.90.123014
"""
module KATV2014

using ..constants: Tₜₕ

"""
Calulates the total inclusive π-0 production cross section using method
of Kafexhiu+ (2014) [PhRvD, V.90, 123014], section 4.
"""
function get_σ_π(Tₚ, i_data, s_ECM)
    # If Tₚ < 2 GeV, calculate cross sections using experimental data and equations (2) - (5)
    # Note that we have already screened for Tₜₕ < Tₚ prior to this subroutine
    #-------------------------------------------------------------------------
    if Tₚ < 2

        eq4_γ = M_res * hypot(M_res, Γ_res)
        eq4_K = √8 * M_res * Γ_res * eq4_γ / ( π * √( M_res^2 + eq4_γ ) )
        f_BW  = rmp * eq4_K  /  ( ((√s_ECM - rmp)^2 - M_res^2)^2 +  M_res^2 * Γ_res^2  )
        eq3_η = √( (s_ECM - rmpi^2 - 4*rmp^2)^2 -  (4*rmpi*rmp)^2 ) /  (2rmpi * √s_ECM)

        # Equation (2)
        σ_1π = 7.66e-3 * eq3_η^1.95 * (1 + eq3_η + eq3_η^5) * f_BW^1.86

        # Double pion production (Equation (5)) can only happen if Tₚ > 2T_th
        σ_2π = Tₚ < 2T_th ? 0.0 : 5.7  /  ( 1 + exp(-9.3 * (Tₚ - 1.4)) )

        # Combine σ_1π and σ_2π
        σ_π = σ_1π  +  σ_2π


    # If 2 Gev < Tₚ < 5 GeV, then use Equations (1) and (6)
    #-------------------------------------------------------------------------
    elseif Tₚ < 5

        eq6_Q = ( Tₚ - Tₜₕ ) / rmp
        n_π0 = -6e-3 + 0.237*eq6_Q - 0.023*eq6_Q^2  # Neutral pion multiplicity
        eq1_ratio  = Tₚ / Tₜₕ
        log_ratio  = log( eq1_ratio )
        σ_inel = ( 30.7  -  0.96*log_ratio  +  0.18*log_ratio^2 ) * ( 1 - eq1_ratio^(-1.9) )^3 # Total inelastic cross section
        # Combine to get σ_π
        σ_π   = n_π0 * σ_inel


    # If 5 GeV < Tₚ, then use Equations (1) and (7). Parameters chosen for
    # Equation (7) depend on values of Tₚ and i_data
    #-------------------------------------------------------------------------
    else

        if i_data == 2 && Tₚ > 5e1
            # Use PYTHIA 8 parametrization
            eq7_a1 = 0.652
            eq7_a2 = 0.0016
            eq7_a3 = 0.488
            eq7_a4 = 0.1928
            eq7_a5 = 0.483
        elseif i_data == 3 && Tₚ > 1e2
            # Use SIBYLL 2.1 parametrization
            eq7_a1 = 5.436
            eq7_a2 = 0.254
            eq7_a3 = 0.072
            eq7_a4 = 0.075
            eq7_a5 = 0.166
        elseif i_data == 4 && Tₚ > 1e2
            # Use QGSJET-I parametrization
            eq7_a1 = 0.908
            eq7_a2 = 0.0009
            eq7_a3 = 6.089
            eq7_a4 = 0.176
            eq7_a5 = 0.448
        else
            # Use GEANT 4 parametrization because of choice of i_data or because
            #   Tₚ too low for validity range of other models
            eq7_a1 = 0.728
            eq7_a2 = 0.596
            eq7_a3 = 0.491
            eq7_a4 = 0.2503
            eq7_a5 = 0.117
        end  # check on i_data and Tₚ

        eq7_ξ = (Tₚ - 3) / rmp
        # Neutral pion multiplicity
        n_π0 = eq7_a1 * eq7_ξ^eq7_a4 * ( 1 + exp( -eq7_a2 * eq7_ξ^eq7_a5 ) ) * ( 1 - exp( -eq7_a3 * eq7_ξ^0.25 ) )

        eq1_ratio  = Tₚ / Tₜₕ
        log_ratio  = log( eq1_ratio )
        σ_inel = ( 30.7  -  0.96*log_ratio  +  0.18*log_ratio^2 ) * ( 1 - eq1_ratio^(-1.9) )^3 # Total inelastic cross section
        # Combine to get σ_π
        σ_π   = n_π0 * σ_inel

    end  # check on Tₚ

    return σ_π
end

const get_sigma_pi = get_σ_π

@doc raw"""
Calulates F(Tₚ, E_γ) as defined in Equations (9) and (11) of Kafexhiu+ (2014) [PhRvD, V.90, 123014].

Equation 9:

```math
Y_γ        = E_γ        + \frac{m_π²}{4E_γ} , \;
Y_γ^{\max} = E_γ^{\max} + \frac{m_π²}{4E_γ^{\max}} , \;
X_γ        = \frac{Y_γ - m_π}{Y_γ^{\max} - m_π}
```

Equation 11:

```math
F(Tₚ, E_γ) = \frac{
    ( 1 - X_γ^{α(Tₚ)} )^{β(Tₚ)}
}{
    ( 1 + X_γ / C )^{γ(Tₚ)}
}
```

Also uses equation 14:

```math
b₀ = 5.9, \;
κ(Tₚ) = 3.29 - \frac{1}{5} θₚ^{-3/2}
```

and equation 15:

```math
μ(Tₚ) = \frac{5}{4} q^{5/4} \exp(- \frac{5}{4} q)
```
"""
function get_Ffunc(Tₚ, Eγ, i_data, Eγ_max)

    # Compute X_γ (independent variable to be used in Equation (11)) using Equation (9)
    Y_γ   = Eγ      +  rmpi^2 / Eγ
    Y_max = Eγ_max  +  rmpi^2 / Eγ_max

    X_γ   = ( Y_γ - rmpi )  /  ( Ymax - rmpi )


    # Now use X_γ, Tₚ, and i_data to calculate F(Tₚ, E_γ)
    #-------------------------------------------------------------------------
    if X_γ < 0 || X_γ > 1
        # Current photon energy cannot be produced at this proton energy
        F_func = 0.0

    elseif Tₚ < 1
        # Use parametrization based on experimental data
        eq14_θ = Tₚ / rmp
        eq14_κ = 3.29 - 0.2 * eq14_θ^(-1.5)

        β = eq14_κ

        F_func = (1 - X_γ)^β
    else
        if Tₚ < 4
            # Use parametrization based on low energy GEANT 4 fit
            eq15_q = (Tₚ - 1) / rmp
            eq15_μ = 1.25 * eq15_q^1.25 * exp( -1.25 * eq15_q )
            λ = 3.0
            α = 1.0
            β = eq15_μ + 2.45
            γ = eq15_μ + 1.45
        elseif Tₚ < 20
            # Use parametrization based on medium energy GEANT 4 fit
            eq15_q = (Tₚ - 1) / rmp
            eq15_μ = 1.25 * eq15_q^1.25 * exp( -1.25 * eq15_q )
            λ = 3.0
            α = 1.0
            β = 1.50*eq15_μ + 4.95
            γ = eq15_μ + 1.5
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

        C = λ * rmpi / Y_max
        F_func = ( 1 - X_γ^α )^β / ( 1 + X_γ / C )^γ
    end
    return Ffunc
end

@doc raw"""
Calculates Amax as defined in Equation (12) of Kafexhiu+ (2014) [PhRvD, V.90, 123014].

```math
A_\max(Tₚ) = \begin{cases}
b₀ ⋅ \frac{σ_π(Tₚ)}{E_π^\max}  & \text{ for } Tₚ^\mathrm{th} ≤ Tₚ < 1\,\mathrm{GeV} \\
b₁ θₚ^{-b₂} \exp(b₃ \log^2 (θₚ)) ⋅ \frac{σ_π(Tₚ)}{mₚ} & \text{ for } Tₚ ≥ 1\,\mathrm{GeV}
\end{cases}
```
"""
function get_Amax(Tₚ, i_data, s_ECM, σ_π)

    # Determine Emax_π_LAB, the maximum pion energy (in the lab frame,
    # where the target nucleon is at rest) allowed by kinematics
    #----------------------------------------------------------------------
    # Maximum π-0 energy in center-of-mass frame
    E_π_CM = (s_ECM - 4rmp^2 + rmpi^2) / (2 * √s_ECM)
    # Lorentz factor of proton in center-of-mass frame; associated β
    γ_CM  = ( Tₚ + 2rmp ) / √s_ECM
    β_CM = √( 1 - 1 / γ_CM^2 )
    # Maximum π-0 momentum in center-of-mass frame. Note the unusual units
    #  of GeV, which will allow easier calculation of Emax_π_LAB
    P_π_CM = √( E_π_CM^2 - rmpi^2 )

    Emax_π_LAB = γ_CM * ( E_π_CM  +  P_π_CM*β_CM )
    #-------------------------------------------------------------------------
    # Emax_π_LAB calculated


    # Use Emax_π_LAB to calculate maximum photon energy in the laboratory frame, Eγ_max
    #-------------------------------------------------------------------------
    # Maximum π-0 Lorentz factor in laboratory frame; associated β
    γ_LAB = Emax_π_LAB / rmpi
    β_LAB = √( 1 - 1/γ_LAB^2 )

    Eγ_max = rmpi/2 * γ_LAB * (1 + β_LAB)
    #-------------------------------------------------------------------------
    # Eγ_max found


    # Now calculate Amax, which depends on values of i_data and Tₚ;
    # we have already screened for Tₜₕ < Tₚ. # Use Equation (12) and Table VII.
    #-------------------------------------------------------------------------
    if Tₚ < 1
        # Based off of experimental data
        Amax = 5.9 * σ_π / Emax_π_LAB

    else
        if i_data == 1 && Tₚ < 5
            # Use GEANT 4 parametrization for low energies
            eq12_b1 = 9.53
            eq12_b2 = 0.52
            eq12_b3 = 0.054
        elseif i_data == 2 && Tₚ > 50
            # Use PYTHIA 8 parametrization
            eq12_b1 = 9.06
            eq12_b2 = 0.3795
            eq12_b3 = 0.01105
        elseif i_data == 3 && Tₚ > 100
            # Use SIBYLL 2.1 parametrization
            eq12_b1 = 10.77
            eq12_b2 = 0.412
            eq12_b3 = 0.01264
        elseif i_data == 4 && Tₚ > 100
            # Use QGSJET-I parametrization
            eq12_b1 = 13.16
            eq12_b2 = 0.4419
            eq12_b3 = 0.01439
        else
            # Use GEANT 4 parametrization at higher energies, because of choice of
            # i_data or because Tₚ too low for validity range of other models
            eq12_b1 = 9.13
            eq12_b2 = 0.35
            eq12_b3 = 0.0097
        end
        eq12_θ = Tₚ / rmp
        log_θ  = log(eq12_θ)
        Amax = eq12_b1 * eq12_θ^(-eq12_b2) * σ_π / rmp * exp( eq12_b3 * log_θ^2 )
    end

    #-------------------------------------------------------------------------
    # Amax determined

    return Eγ_max, Amax
end
end # module
