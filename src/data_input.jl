using TOML

# FIXME the purpose of this file is to be include()-d in mc_cr.jl
# This populates the global variables.
# XXX THIS IS A BAD DESIGN.
# TODO Design another way to pass data around -> maybe all const?

cfg_toml = TOML.parsefile("mc_in.toml")

const u₀, β₀, γ₀ = let

    skspd = cfg_toml["SKSPD"]
    skspd > 0 || error("Shock speed must be positive")

    skspd_unit = cfg_toml["SKSPD_UNIT"]
    if skspd_unit ∈ ("gamma", "γ")
        skspd > 1 || error("SKSPD: Lorentz factor must be > 1")
        γ = skspd
        β = √(1 - 1/γ^2)
        u = β * c |> cm/s
    else
        if skspd_unit == "km/s"
            0 < skspd < ustrip(km/s, Unitful.c0) || error("SKSPD: u must be between 0 and c")
            u = (skspd * 1e5) * cm/s
            β = u / c
        elseif skspd_unit == "c"
            0 < skspd < 1 || error("SKSPD: β must be between 0 and 1")
            β = skspd
            u = β * c |> cm/s
        else
            error("SKSPD: unknown units provided with SKSPD_UNIT")
        end
        γ = lorentz(β)
    end

    (u, β, γ)
end

const species = let
    masses = cfg_toml["AA_ION"] # species mass in units of proton mass
    electron_index = findfirst(isnan, masses)
    masses[electron_index] = NoUnits(me/mp) # electron mass over proton mass

    charges = cfg_toml["ZZ_ION"]
    charges[electron_index] = -1

    temperatures = cfg_toml["TZ_ION"] # temperature of each species
    densities = cfg_toml["DENZ_ION"] # number density of each species

    if !(length(masses) == length(charges) == length(temperatures) == length(densities))
        error("Inconsistent number of ion parameters given (AA_ION, ZZ_ION, TZ_ION, DENZ_ION)")
    end

    Species.(masses*mp, charges*qcgs, temperatures*K, densities/cm^3)
end
const n_ions = length(species)

const inp_distr = cfg_toml["INDST"]
const energy_inj = cfg_toml["ENINJ"] * keV
const inj_weight = get(cfg_toml, "INJWT", true)

const Emax_keV, Emax_keV_per_aa, pmax_cgs = let
    energy_max = cfg_toml["ENMAX"]
    if energy_max[1] > 0      # All species have same max energy
        Emax_keV        = energy_max[1]
        Emax_keV_per_aa = 0.0
        pmax_cgs        = 0.0

    elseif energy_max[2] > 0  # Max energy depends on aa
        Emax_keV        = 0.0
        Emax_keV_per_aa = energy_max[2]
        pmax_cgs        = 0.0

    elseif energy_max[3] > 0  # All species have same max momentum
        Emax_keV        = 0.0
        Emax_keV_per_aa = 0.0
        pmax_cgs        = energy_max[3]
    else
        error("ENMAX: at least one choice must be non-zero.")
    end
    (Emax_keV*keV, Emax_keV_per_aa*keV, pmax_cgs*mp*c)
end

const η_mfp = get(cfg_toml, "GYFAC", 1)


const bmag₀ = cfg_toml["BMAGZ"]*G
# rg₀ below is the gyroradius of a proton whose speed is u₀ that is gyrating in a field
# of strength bmag₀. Note that this formula is relativistically correct
const rg₀ = (γ₀ * E₀_proton * β₀) / (qcgs * bmag₀) |> cm


begin
    const θ_B₀ = cfg_toml["THTBZ"] # must be zero
    if θ_B₀ > 0
        error("program cannot currently handle oblique shocks. Adjust THTBZ.")
    elseif θ_B₀ < 0
        error("unphysical value for THTBZ. Must be at least 0.")
    end
end

begin
    const x_grid_start_rg = cfg_toml["XGDUP"]
    const x_grid_stop_rg  = cfg_toml["XGDDW"]
    x_grid_start_rg ≥ 0 && error("XGDUP: x_grid_start must be negative.")
    x_grid_stop_rg  ≤ 0 && error("XGDDW: x_grid_stop must be positive.")
end

const feb_UpS = let
    febup = get(cfg_toml, "FEBUP", nothing)
    if isnothing(febup)
        feb_UpS = x_grid_start_rg * rg₀ # default value
        return
    end
    if febup[1] < 0
        feb_UpS = febup[1] * rg₀
    elseif febup[2] < 0
        feb_UpS = uconvert(cm, febup[2] * pc)
    else
        error("FEBUP: at least one choice must be negative.")
    end
    ((feb_UpS/rg₀) < x_grid_start_rg) && error("FEBUP: UpS FEB must be within x_grid_start")

    feb_UpS
end

const feb_DwS, use_prp = let
    febdw = get(cfg_toml, "FEBDW", nothing)
    use_prp = false
    if isnothing(febdw)
        feb_DwS = -1 # default value
        return
    end
    if febdw[1] > 0
        feb_DwS = febdw[1] * rg₀
    elseif febdw[2] > 0
        feb_DwS = uconvert(cm, febdw[2] * pc)
    else
        feb_DwS = 0.0cm
        use_prp = true
    end
    (feb_DwS, use_prp)
end

begin
    const x_spec = get(cfg_toml, "XSPEC", Float64[])
    const n_xspec = length(x_spec)
end

const n_itrs = cfg_toml["NITRS"]
const xn_per_coarse = cfg_toml["XN_PER_COARSE"]
const xn_per_fine = cfg_toml["XN_PER_FINE"]

begin
    const n_pts_inj = cfg_toml["N_PTS_INJ"]
    const n_pts_pcut = cfg_toml["N_PTS_PCUT"]
    max(n_pts_inj,n_pts_pcut) > na_particles && error("Array size na_particles too small.")
end

begin
    const n_pts_pcut_hi = cfg_toml["N_PTS_PCUT_HI"]
    const energy_pcut_hi = cfg_toml["EN_PCUT_HI"]
    n_pts_pcut_hi > na_particles && error("Array size na_particles too small.")
end

begin
    const pcuts_in = cfg_toml["PCUTS"]
    const n_pcuts = length(pcuts_in)
    n_pcuts+1 > na_c && error("PCUTS: parameter na_c smaller than desired number of pcuts.")

    if Emax_keV > 0keV
        # Convert from momentum[mₚc/aa] to energy[keV]
        Emax_eff = 56 * pcuts_in[n_pcuts-1] * ustrip(keV, E₀_proton*erg)

        if Emax_keV > Emax_eff
            error("PCUTS: max energy exceeds highest pcut. Add more pcuts or lower Emax_keV. ",
                  "Emax_keV (assuming Fe) = $Emax_keV; Emax_eff = $Emax_eff")
        end
    elseif Emax_keV_per_aa > 0keV   # Limit was on energy per nucleon
        # Convert from momentum[mₚc/aa] to energy[keV/aa]
        Emax_eff = pcuts_in[n_pcuts-1] * ustrip(keV, E₀_proton*erg)

        if Emax_keV_per_aa > Emax_eff
            error("PCUTS: max energy per aa exceeds highest pcut. Add more pcuts or lower Emax_keV_per_aa. ",
                  "Emax_keV_per_aa = $Emax_keV_per_aa; Emax_eff/aa = $Emax_eff")
        end

    elseif pmax_cgs > 0mp*c # Limit was on total momentum. Assume Fe for strictest limit on mom/nuc.
        pmax_eff = 56mp*c * pcuts_in[n_pcuts-1]
        if pmax_cgs > pmax_eff
            error("PCUTS: max momentum exceeds highest pcut. Add more pcuts or lower pmax. ",
                  "pmax[m_pc] = $pmax_cgs; pmax_eff (for Fe) = $pmax_eff")
        end
    else   # Something unexpected has happened
        error("Unexpected result when comparing pcut max to energy/momentum max")
    end
end

const dont_shock = (get(cfg_toml, "NOSHK", 0) == 66)

const dont_scatter = (get(cfg_toml, "NOSCT", 0) == 66)

const dont_DSA = (get(cfg_toml, "NODSA", 0) == 66)

const do_smoothing = (cfg_toml["SMSHK"] != 66)

const prof_weight_fac = get(cfg_toml, "SMIWT", 1.0)

const do_prof_fac_damp = (get(cfg_toml, "SMVWT", 0) == 66)

begin
    const smooth_mom_energy_fac = get(cfg_toml, "SMMOE", 0.0)
    if smooth_mom_energy_fac < 0 || smooth_mom_energy_fac > 1
        throw(DomainError(smooth_mom_energy_fac, "smooth_mom_energy_fac/SMMOE must be in [0, 1]"))
    end
end

begin
    const smooth_pressure_flux_psd_fac = get(cfg_toml, "SMPFP", 0)
    if smooth_pressure_flux_psd_fac < 0 || smooth_pressure_flux_psd_fac > 1
        throw(DomainError(smooth_pressure_flux_psd_fac,
                          "smooth_pressure_flux_psd_fac/SMPFP must be in [0, 1]"))
    end
    # TODO: actually get pressure calculation working properly
    if smooth_pressure_flux_psd_fac > 0
        error("SMPFP: code does not properly calculate pressure from PSD. ",
              "Set to 0 or get this code working")
    end
end

const r_comp, r_RH, Γ₂_RH = let
    r_comp = cfg_toml["RCOMP"]
    r_RH, Γ₂_RH = calc_rRH(u₀, β₀, γ₀, species)
    if r_comp == -1
        r_comp = r_RH
    end
    (r_comp, r_RH, Γ₂_RH) # shadowed variables
end

begin
    β₂, γ₂, bmag₂, θ_B₂, θᵤ₂ = calc_DwS(bmag₀, r_comp, β₀)
    const u₂ = β₂*c
    @debug("Results from calc_DwS()", u₂, β₂, γ₂, bmag₂, θ_B₂, θᵤ₂)
end

begin
    const do_old_prof = (get(cfg_toml, "OLDIN", 0) == 66)
    if do_old_prof
        n_old_skip, n_old_profs, n_old_per_prof = cfg_toml["OLDDT"]
    else
        n_old_skip, n_old_profs, n_old_per_prof = [0, 0, 0]
    end
end

const age_max = let
    age_max = get(cfg_toml, "AGEMX", -1.0)
    if age_max < 0
        age_max = -1.0
    end
    age_max * s
end
# default behavior of do_retro is dependent on age_max
const do_retro = (get(cfg_toml, "RETRO", age_max > 0s ? 66 : 0) == 66)

const do_fast_push = (get(cfg_toml, "FPUSH", 0) == 66)
const x_fast_stop_rg = do_fast_push ? cfg_toml["FPSTP"] : 0.0

const x_art_start_rg, x_art_scale = let
    art = get(cfg_toml, "ARTSM", nothing)
    if isnothing(art)
        x_art_start_rg = 0.0
        x_art_scale = 0.0
    else
        x_art_start_rg, x_art_scale = art
    end
    (x_art_start_rg, x_art_scale)
end

const pₑ_crit, γₑ_crit = let

    energyₑ_crit_keV = get(cfg_toml, "EMNFP", nothing) * keV
    # If needed, convert input energy to momentum and Lorentz factor
    if !isnothing(energyₑ_crit_keV) && energyₑ_crit_keV > 0keV
        energyₑ_crit_rm = energyₑ_crit_keV / E₀_electron

        # Different forms for nonrelativistic and relativstic momenta
        if energyₑ_crit_rm < 1e-2
            pₑ_crit = me*c * √(2energyₑ_crit_rm)
            γₑ_crit = 1.0
        else
            pₑ_crit = me*c * √((energyₑ_crit_rm + 1)^2 - 1)
            γₑ_crit = energyₑ_crit_rm + 1.0
        end
    else
        energyₑ_crit_keV = -1.0
        pₑ_crit = -1.0
        γₑ_crit = -1.0
    end
    (pₑ_crit, γₑ_crit)
end

const do_rad_losses = (get(cfg_toml, "NORAD", 0) != 66)

const do_photons = (get(cfg_toml, "PHOTN", 0) == 66)

# JETRD only mandatory if doing photons
const jet_rad_pc = do_photons ? cfg_toml["JETRD"] : get(cfg_toml, "JETRD", 0.0)

const jet_sph_frac, jet_open_ang_deg = let

    jetfr = get(cfg_toml, "JETFR", nothing)
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
    (jet_sph_frac, jet_open_ang_deg)
end

begin
    jet_dist_kpc = get(cfg_toml, "JETDS", 1.0)
    redshift = get(cfg_toml, "RDSHF", 0.0)
    if jet_dist_kpc > 0 && redshift > 0
        error("JETDS: At most one of 'JETDS' and 'RDSHF' may be non-zero.")
    end
end
begin
    # The following option is not in the Fortran program
    const cosmo_var = cfg_toml["COSMO_VAR"]
    cosmo_var ≠ 1 && cosmo_var ≠ 2 && error("Invalid value for cosmo_var")
end

begin
    const energy_transfer_frac = float(get(cfg_toml, "ENXFR", 0.0))
    if energy_transfer_frac < 0 || energy_transfer_frac > 1
        error("ENXFR: energy_transfer_frac must be in [0,1]")
    end
end

const num_UpS_shells, num_DwS_shells = cfg_toml["NSHLS"]

begin
    const bturb_comp_frac = get(cfg_toml, "BTRBF", 0.0)
    const bfield_amp = get(cfg_toml, "BAMPF", 1.0)
    bfield_amp < 1 && error("BAMPF: must be ≥ 1.d0")
    if bfield_amp > 1 && iszero(bturb_comp_frac)
        error("BTRBF: bfield_amp > 1 has no effect if BTRBF = 0")
    end
end

const psd_bins_per_dec_mom, psd_bins_per_dec_θ = let
    psd_bins = get(cfg_toml, "PSDBD", [10, 10])
    psd_bins_per_dec_mom::Int = psd_bins[1]
    psd_bins_per_dec_θ::Int   = psd_bins[2]
    if psd_bins_per_dec_mom ≤ 0 || psd_bins_per_dec_θ ≤ 0
        error("PSDBD: both values must be positive.")
    end
    (psd_bins_per_dec_mom::Int, psd_bins_per_dec_θ::Int)
end

const psd_lin_cos_bins, psd_log_θ_decs = let
    psd_bins = get(cfg_toml, "PSDTB", [119, 4])
    psd_lin_cos_bins = psd_bins[1]
    psd_log_θ_decs   = psd_bins[2]
    if psd_lin_cos_bins ≤ 0 || psd_log_θ_decs ≤ 0
        error("PSDTB: both values must be positive.")
    end
    (psd_lin_cos_bins::Int, psd_log_θ_decs::Int)
end

const use_custom_frg = (get(cfg_toml, "NWFRG", 0) == 66)

const emin_therm_fac = get(cfg_toml, "EMNFC", 0.01)

const do_multi_dNdps = (get(cfg_toml, "DNDPS", 0) == 66)

const do_tcuts, tcuts, n_tcuts = let
    do_tcuts = haskey(cfg_toml, "TCUTS")
    if do_tcuts
        tcuts = cfg_toml["TCUTS"] * s
        n_tcuts = length(tcuts)

        age_max < 0s && error("tcut tracking must be used with anaccel time limit. Adjust keyword 'AGEMX'.")
        # Check to make sure we haven't used more tcuts than allowed by na_c
        (n_tcuts+1) > na_c && error("TCUTS: parameter na_c smaller than desired number of tcuts.")
        # Check to make sure final tcut is much larger than age_max so that
        #   we never have to worry about exceeding it
        tcuts[end] ≤ 10age_max && error("TCUTS: final tcut must be much (10x) larger than age_max.")
    else
        tcuts = TimeCGS[]
        n_tcuts = 0
    end
    (do_tcuts, tcuts, n_tcuts)
end

begin
    const inj_fracs = get(cfg_toml,"INJFR", fill(1.0, n_ions))
    length(inj_fracs) == n_ions || error("Number of injection probabilities must match NIONS")
end

const use_custom_εB = (get(cfg_toml, "NWEPB", 0) == 66)
