function parse_shock_speed(skspd, skspd_unit)
    skspd > 0 || error("Shock speed must be positive")

    if skspd_unit ∈ ("gamma", "γ")
        skspd > 1 || error("shock-speed: Lorentz factor must be > 1")
        γ = skspd
        β = √(1 - 1/γ^2)
        u = β * c |> cm/s
    else
        if skspd_unit == "km/s"
            0 < skspd < ustrip(km/s, Unitful.c0) || error("shock-speed: u must be between 0 and c")
            u = (skspd * 1e5) * cm/s
            β = u / c
        elseif skspd_unit == "c"
            0 < skspd < 1 || error("shock-speed: β must be between 0 and 1")
            β = skspd
            u = β * c |> cm/s
        else
            error("shock-speed: unknown units provided with shock-speed-unit")
        end
        γ = lorentz(β)
    end

    return (u, β, γ)
end

function parse_maximum_energy(energy_max)
    if energy_max[1] > 0      # All species have same max energy
        Emax        = energy_max[1]
        Emax_per_aa = 0.0
        pmax        = 0.0

    elseif energy_max[2] > 0  # Max energy depends on aa
        Emax        = 0.0
        Emax_per_aa = energy_max[2]
        pmax        = 0.0

    elseif energy_max[3] > 0  # All species have same max momentum
        Emax        = 0.0
        Emax_per_aa = 0.0
        pmax        = energy_max[3]
    else
        error("ENMAX: at least one choice must be non-zero.")
    end
    return (Emax*keV, Emax_per_aa*keV, pmax*mp*c)
end

function parse_electron_critical_energy(Eₑ_crit)
    if isnothing(Eₑ_crit) || Eₑ_crit ≤ 0
        return (-1.0*me*c, -1.0)
    end

    Eₑ_crit *= keV # attach units
    # Convert input energy to momentum and Lorentz factor
    Eₑ_crit_rm = Eₑ_crit / E₀ₑ

    # Different forms for nonrelativistic and relativstic momenta
    if Eₑ_crit_rm < 1e-2
        pₑ_crit = √(2 * me * Eₑ_crit)
        γₑ_crit = 1.0
    else
        pₑ_crit = me*c * √((Eₑ_crit_rm + 1)^2 - 1)
        γₑ_crit = Eₑ_crit_rm + 1
    end
    return (pₑ_crit |> g*cm/s, γₑ_crit)
end

function check_shock_angle(θ)
    if θ > 0
        error("program cannot currently handle oblique shocks. Adjust theta-B0.")
    elseif θ < 0
        error("unphysical value for theta-B0. Must be at least 0.")
    end
end

function check_x_grid_limits(x_grid_start_rg, x_grid_stop_rg)
    x_grid_start_rg ≥ 0 && error("x_grid_limits: x_grid_start must be negative.")
    x_grid_stop_rg  ≤ 0 && error("x_grid_limits: x_grid_stop must be positive.")
end

function check_pcuts_in(pcuts_in, Emax, Emax_per_aa, pmax)
    n_pcuts = length(pcuts_in)
    n_pcuts+1 > na_c && error("momentum-cutoffs: parameter na_c smaller than desired number of pcuts.")

    if Emax > 0keV
        # Convert from momentum[mₚc/aa] to energy[keV]
        Emax_eff = 56 * pcuts_in[n_pcuts-1] * ustrip(keV, E₀ₚ*erg)

        if Emax > Emax_eff
            error("PCUTS: max energy exceeds highest pcut. Add more pcuts or lower Emax. ",
                  "Emax (assuming Fe) = $Emax; Emax_eff = $Emax_eff")
        end
    elseif Emax_per_aa > 0keV   # Limit was on energy per nucleon
        # Convert from momentum[mₚc/aa] to energy[keV/aa]
        Emax_eff_per_aa = pcuts_in[n_pcuts-1] * ustrip(keV, E₀ₚ*erg)

        if Emax_per_aa > Emax_eff_per_aa
            error("PCUTS: max energy per aa exceeds highest pcut. Add more pcuts or lower Emax_per_aa. ",
                  "Emax_per_aa = $Emax_per_aa; Emax_eff/aa = $Emax_eff_per_aa")
        end

    elseif pmax > 0mp*c # Limit was on total momentum. Assume Fe for strictest limit on mom/nuc.
        pmax_eff = 56mp*c * pcuts_in[n_pcuts-1]
        if pmax > pmax_eff
            error("PCUTS: max momentum exceeds highest pcut. Add more pcuts or lower pmax. ",
                  "pmax[m_pc] = $pmax; pmax_eff (for Fe) = $pmax_eff")
        end
    else   # Something unexpected has happened
        error("Unexpected result when comparing pcut max to energy/momentum max")
    end
end
