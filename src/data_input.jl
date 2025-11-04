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

function check_pcuts(pcuts, Emax, Emax_per_aa, pmax)
    length(pcuts) > na_c && error("momentum-cutoffs: parameter na_c smaller than desired number of pcuts.")

    if Emax > 0keV
        # Convert from momentum[mₚc/aa] to energy[keV]
        Emax_eff = 56 * pcuts[end] * c

        if Emax > Emax_eff
            error("PCUTS: max energy exceeds highest pcut. Add more pcuts or lower Emax. ",
                  "Emax (assuming Fe) = $Emax; Emax_eff = $Emax_eff")
        end
    elseif Emax_per_aa > 0keV   # Limit was on energy per nucleon
        # Convert from momentum[mₚc/aa] to energy[keV/aa]
        Emax_eff_per_aa = pcuts[end]*c

        if Emax_per_aa > Emax_eff_per_aa
            error("PCUTS: max energy per aa exceeds highest pcut. Add more pcuts or lower Emax_per_aa. ",
                  "Emax_per_aa = $Emax_per_aa; Emax_eff/aa = $Emax_eff_per_aa")
        end

    elseif pmax > 0g*c # Limit was on total momentum. Assume Fe for strictest limit on mom/nuc.
        pmax_eff = 56 * pcuts[end]
        if pmax > pmax_eff
            error("PCUTS: max momentum exceeds highest pcut. Add more pcuts or lower pmax. ",
                  "pmax = $pmax; pmax_eff (for Fe) = $pmax_eff")
        end
    else   # Something unexpected has happened
        error("Unexpected result when comparing pcut max to energy/momentum max")
    end
end

function get_feb(febup, febdw, x_grid_start_rg, rg₀)
    if isnothing(febup)
        feb_upstream = x_grid_start_rg * rg₀ # default value
    else
        if febup[1] < 0
            feb_upstream = febup[1] * rg₀
        elseif febup[2] < 0
            feb_upstream = uconvert(cm, febup[2] * pc)
        else
            error("FEB-upstream: at least one choice must be negative.")
        end
        (feb_upstream/rg₀ < x_grid_start_rg) && error("FEB-upstream: upstream FEB must be within x_grid_start")
    end

    use_prp = false
    if isnothing(febdw)
        feb_downstream = -1 # default value
    else
        if febdw[1] > 0
            feb_downstream = febdw[1] * rg₀
        elseif febdw[2] > 0
            feb_downstream = uconvert(cm, febdw[2] * pc)
        else
            feb_downstream = 0.0cm
            use_prp = true
        end
    end
    return (feb_upstream, feb_downstream, use_prp)
end

function parse_jet_frac(jetfr, do_photons=false)
    if isnothing(jetfr) # default behavior, handled differently based on PHOTNS
        do_photons && error("If calculating photons, 'JETFR' must be specified manually.")
        jet_sph_frac     = 0.0
        jet_open_ang_deg = 0.0
    elseif 0 < jetfr[1] ≤ 1
        jet_sph_frac     = jetfr[1]
        jet_open_ang_deg = acosd(1 - 2jet_sph_frac)
    elseif 0 < jetfr[2] ≤ 180
        jet_open_ang_deg = jetfr[2]
        jet_sph_frac     = (1 - cosd(jet_open_ang_deg)) / 2
    else
        error("JETFR: Unphysical values entered.")
    end
    return (jet_sph_frac, jet_open_ang_deg)
end

function parse_species(cfg)
    masses = cfg["AA_ION"] # species mass in units of proton mass
    electron_index = findfirst(isnan, masses)
    masses[electron_index] = NoUnits(me/mp) # electron mass over proton mass

    charges = cfg["ZZ_ION"]
    charges[electron_index] = -1

    temperatures = cfg["TZ_ION"] # temperature of each species
    densities = cfg["DENZ_ION"] # number density of each species

    if !(length(masses) == length(charges) == length(temperatures) == length(densities))
        error("Inconsistent number of ion parameters given (AA_ION, ZZ_ION, TZ_ION, DENZ_ION)")
    end

    return Species.(masses*mp, charges*qcgs, temperatures*K, densities/cm^3)
end
