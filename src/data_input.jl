using TOML

# FIXME the purpose of this file is to be include()-d in mc_cr.jl
# This populates the global variables.
# XXX THIS IS A BAD DESIGN.
# TODO Design another way to pass data around

cfg_toml = TOML.parsefile("mc_in.toml")

let skspd = cfg_toml["SKSPD"]
    global u₀, γ₀, β₀
    if 0 < skspd[1] < c_cgs*1e-5    # km/sec
        u₀ = skspd[1] * 1e5
        β₀ = u₀ / c_cgs
        γ₀ = 1 / √(1 - γ₀^2)
    elseif skspd[2] > 1             # Lorentz factor
        γ₀ = skspd[2]
        β₀ = √(1 - 1/γ₀^2)
        u₀ = β₀ * c_cgs
    elseif 0 < skspd[3] < 1         # in units of c
        β₀ = skspd[3]
        u₀ = β₀ * c_cgs
        γ₀ = 1 / √(1 - β₀^2)
    else
        error("SKSPD: at least one choice must be non-zero and physically reasonable.")
    end
end

let
    global aa_ion = cfg_toml["AA_ION"] # species mass in units of proton mass
    global n_ions = length(aa_ion)

    electron_index = findfirst(isnan, aa_ion)
    global m_ion = aa_ion * mₚ_cgs # species mass in grams

    aa_ion[electron_index] = mₑ_cgs / mₚ_cgs
    m_ion[electron_index]  = mₑ_cgs

    global zz_ion = cfg_toml["ZZ_ION"]
    zz_ion[aa_ion .< 1] .= 1

    global T₀_ion = cfg_toml["TZ_ION"] # temperature of each species
    global ρ_N₀_ion = cfg_toml["DENZ_ION"] # number density of each species

    if any(!=(n_ions), (length(zz_ion), length(T₀_ion), length(ρ_N₀_ion)))
        error("Inconsistent number of ion parameters given (AA_ION, ZZ_ION, TZ_ION, DENZ_ION)")
    end

    @debug("Species parameters", n_ions, aa_ion, m_ion, zz_ion, T₀_ion, ρ_N₀_ion)
end

inp_distr = cfg_toml["INDST"]

energy_inj = cfg_toml["ENINJ"]

inj_weight = get(cfg_toml, "INJWT", true)

let energy_max = cfg_toml["ENMAX"]
    global Emax_keV, Emax_keV_per_aa, pmax_cgs
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
        pmax_cgs        = energy_max[3] * mₚ_cgs * c_cgs
    else
        error("ENMAX: at least one choice must be non-zero.")
    end
end

η_mfp = get(cfg_toml, "GYFAC", 1)

begin # just for grouping and better understanding of code, no actual scoping rules here
    bmag₀ = cfg_toml["BMAGZ"]
    # rg₀ below is the gyroradius of a proton whose speed is u₀ that is gyrating in a field
    # of strength bmag₀. Note that this formula is relativistically correct
    rg₀ = (γ₀ * mₚ_cgs * c_cgs^2 * β₀) / (qₚ_cgs * bmag₀)
end

begin
    θ_B₀ = cfg_toml["THTBZ"]
    if θ_B₀ > 0
        oblique = true
        error("program cannot currently handle oblique shocks. Adjust THTBZ.")
    elseif θ_B₀ < 0
        error("unphysical value for THTBZ. Must be at least 0.")
    else
        oblique = false
    end
end

begin
    x_grid_start_rg = cfg_toml["XGDUP"]
    x_grid_stop_rg  = cfg_toml["XGDDW"]
    x_grid_start_rg ≥ 0 && error("XGDUP: x_grid_start must be negative.")
    x_grid_stop_rg  ≤ 0 && error("XGDDW: x_grid_stop must be positive.")
end

let febup = get(cfg_toml, "FEBUP", nothing)
    global feb_UpS
    if isnothing(febup)
        feb_UpS = x_grid_start_rg * rg₀ # default value
        return
    end
    if febup[1] < 0
        feb_UpS = febup[1] * rg₀
    elseif febup[2] < 0
        feb_UpS = ustrip(u"cm", febup[2] * u"pc")
    else
        error("FEBUP: at least one choice must be negative.")
    end
    ((feb_UpS/rg₀) < x_grid_start_rg) && error("FEBUP: UpS FEB must be within x_grid_start")
end

let febdw = get(cfg_toml, "FEBDW", nothing)
    global feb_DwS
    global use_prp = false
    if isnothing(febdw)
        feb_DwS = -1 # default value
        return
    end
    if febdw[1] > 0
        feb_DwS = febdw[1] * rg₀
    elseif febdw[2] > 0
        feb_DwS = ustrip(u"cm", febdw[2] * u"pc")
    else
        feb_DwS = 0.0
        use_prp = true
    end
end

begin
    x_spec = get(cfg_toml, "XSPEC", Float64[])
    n_xspec = length(x_spec)
end

n_itrs = cfg_toml["NITRS"]
xn_per_coarse = cfg_toml["XN_PER_COARSE"]
xn_per_fine = cfg_toml["XN_PER_FINE"]

begin
    n_pts_inj = cfg_toml["N_PTS_INJ"]
    n_pts_pcut = cfg_toml["N_PTS_PCUT"]
    max(n_pts_inj,n_pts_pcut) > na_particles && error("Array size na_particles too small.")
end

begin
    n_pts_pcut_hi = cfg_toml["N_PTS_PCUT_HI"]
    energy_pcut_hi = cfg_toml["EN_PCUT_HI"]
    n_pts_pcut_hi > na_particles && error("Array size na_particles too small.")
end

begin
    pcuts_in = cfg_toml["PCUTS"]
    n_pcuts = length(pcuts_in)
    n_pcuts+1 > na_c && error("PCUTS: parameter na_c smaller than desired number of pcuts.")

    if Emax_keV > 0
        # Convert from momentum[mₚc/aa] to energy[keV]
        Emax_eff = 56 * pcuts_in[n_pcuts-1] * ustrip(u"keV", E₀_proton*u"erg")

        if Emax_keV > Emax_eff
            error("PCUTS: max energy exceeds highest pcut. Add more pcuts or lower Emax_keV. ",
                  "Emax_keV (assuming Fe) = $Emax_keV; Emax_eff = $Emax_eff")
        end
    elseif Emax_keV_per_aa > 0   # Limit was on energy per nucleon
        # Convert from momentum[mₚc/aa] to energy[keV/aa]
        Emax_eff = pcuts_in[n_pcuts-1] * ustrip(u"keV", E₀_proton*u"erg")

        if Emax_keV_per_aa > Emax_eff
            error("PCUTS: max energy per aa exceeds highest pcut. Add more pcuts or lower Emax_keV_per_aa. ",
                  "Emax_keV_per_aa = $Emax_keV_per_aa; Emax_eff/aa = $Emax_eff")
        end

    elseif pmax_cgs > 0 # Limit was on total momentum. Assume Fe for strictest limit on mom/nuc.
        pmax_eff = 56*mₚ_cgs*c_cgs * pcuts_in[n_pcuts-1]
        if pmax_cgs > pmax_eff
            error("PCUTS: max momentum exceeds highest pcut. Add more pcuts or lower pmax. ",
                  "pmax[m_pc] = $pmax_cgs; pmax_eff (for Fe) = $pmax_eff")
        end
    else   # Something unexpected has happened
        error("Unexpected result when comparing pcut max to energy/momentum max")
    end
end

dont_shock = (get(cfg_toml, "NOSHK", 0) == 66)

dont_scatter = (get(cfg_toml, "NOSCT", 0) == 66)

dont_DSA = (get(cfg_toml, "NODSA", 0) == 66)

do_smoothing = (cfg_toml["SMSHK"] != 66)

prof_weight_fac = get(cfg_toml, "SMIWT", 1.0)

do_prof_fac_damp = (get(cfg_toml, "SMVWT", 0) == 66)

begin
    smooth_mom_energy_fac = get(cfg_toml, "SMMOE", 0.0)
    if smooth_mom_energy_fac < 0 || smooth_mom_energy_fac > 1
        throw(DomainError(smooth_mom_energy_fac, "smooth_mom_energy_fac/SMMOE must be in [0, 1]"))
    end
end

begin
    smooth_pressure_flux_psd_fac = get(cfg_toml, "SMPFP", 0)
    if smooth_pressure_flux_psd_fac < 0 || smooth_pressure_flux_psd_fac > 1
        throw(DomainError(smooth_pressure_flux_psd_fac, "smooth_pressure_flux_psd_fac/SMPFP must be in [0, 1]"))
    end
    # TODO: actually get pressure calculation working properly
    if smooth_pressure_flux_psd_fac > 0
        error("SMPFP: code does not properly calculate pressure from PSD. ",
              "Set to 0 or get this code working")
    end
end

begin #let
    #global r_comp, r_RH, Γ₂_RH, β₂, γ₂, bmag₂, θ_B₂, θᵤ₂, u₂
    r_comp = cfg_toml["RCOMP"]
    r_RH, Γ₂_RH = calc_rRH(β₀, γ₀, n_ions, aa_ion, ρ_N₀_ion, T₀_ion, oblique)
    if r_comp == -1
        r_comp = r_RH
    end
    β₂, γ₂, bmag₂, θ_B₂, θᵤ₂ = calc_DwS(oblique, bmag₀, r_comp, β₀)
    u₂ = β₂ * c_cgs
    @debug("Results from calc_DwS()", u₂, β₂, γ₂, bmag₂, θ_B₂, θᵤ₂)
end

begin
    do_old_prof = (get(cfg_toml, "OLDIN", 0) == 66)
    if do_old_prof
        n_old_skip, n_old_profs, n_old_per_prof = cfg_toml["OLDDT"]
    else
        n_old_skip, n_old_profs, n_old_per_prof = [0, 0, 0]
    end
end

begin
    age_max = get(cfg_toml, "AGEMX", -1.0)
    if age_max < 0
        age_max = -1.0
    end
    # default behavior of do_retro is dependent on age_max
    do_retro = (get(cfg_toml, "RETRO", age_max > 0 ? 66 : 0) == 66)
end

begin
    do_fast_push = (get(cfg_toml, "FPUSH", 0) == 66)
    x_fast_stop_rg = do_fast_push ? cfg_toml["FPSTP"] : 0.0
end

let art = get(cfg_toml, "ARTSM", nothing)
    global x_art_start_rg, x_art_scale
    if isnothing(art)
        x_art_start_rg = 0.0
        x_art_scale = 0.0
    else
        x_art_start_rg, x_art_scale = art
    end
end

let energy_electron_crit_keV = get(cfg_toml, "EMNFP", nothing)
    global p_electron_crit, γ_electron_crit
    # If needed, convert input energy[keV] to momentum and Lorentz factor
    if !isnothing(energy_electron_crit_keV) && energy_electron_crit_keV > 0
        energy_electron_crit_rm = ustrip(u"erg", energy_electron_crit_keV*u"keV") / E₀_electron

        # Different forms for nonrel and rel momenta
        if energy_electron_crit_rm < 1e-2
            p_electron_crit = (mₑ_cgs*c_cgs) * √(2energy_electron_crit_rm)
            γ_electron_crit = 1.0
        else
            p_electron_crit = (mₑ_cgs*c_cgs) * √((energy_electron_crit_rm + 1)^2 - 1)
            γ_electron_crit = energy_electron_crit_rm + 1.0
        end
    else
        energy_electron_crit_keV = -1.0
        p_electron_crit = -1.0
        γ_electron_crit = -1.0
    end
end

do_rad_losses = (get(cfg_toml, "NORAD", 0) != 66)

do_photons = (get(cfg_toml, "PHOTN", 0) == 66)

# JETRD only mandatory if doing photons
jet_rad_pc = do_photons ? cfg_toml["JETRD"] : get(cfg_toml, "JETRD", 0.0)

let jetfr = get(cfg_toml, "JETFR", nothing)
    global jet_sph_frac, jet_open_ang_deg
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
end

begin
    jet_dist_kpc = get(cfg_toml, "JETDS", 1.0)
    redshift = get(cfg_toml, "RDSHF", 0.0)
    if jet_dist_kpc > 0 && redshift > 0
        error("JETDS: At most one of 'JETDS' and 'RDSHF' may be non-zero.")
    end
    # The following option is not in the Fortran program
    cosmo_var = cfg_toml["COSMO_VAR"]
    if cosmo_var ≠ 1 && cosmo_var ≠ 2
        error("Invalid value for cosmo_var")
    end
end

begin
    energy_transfer_frac = float(get(cfg_toml, "ENXFR", 0.0))
    if energy_transfer_frac < 0 || energy_transfer_frac > 1
        error("ENXFR: energy_transfer_frac must be in [0,1]")
    end
end

num_UpS_shells, num_DwS_shells = cfg_toml["NSHLS"]

begin
    bturb_comp_frac = get(cfg_toml, "BTRBF", 0.0)
    bfield_amp = get(cfg_toml, "BAMPF", 1.0)
    bfield_amp < 1 && error("BAMPF: must be ≥ 1.d0")
    if bfield_amp > 1 && iszero(bturb_comp_frac)
        error("BTRBF: bfield_amp > 1 has no effect if BTRBF = 0")
    end
end

let psd_bins = get(cfg_toml, "PSDBD", [10, 10])
    global psd_bins_per_dec_mom::Int = psd_bins[1]
    global psd_bins_per_dec_θ::Int   = psd_bins[2]
    if psd_bins_per_dec_mom ≤ 0 || psd_bins_per_dec_θ ≤ 0
        error("PSDBD: both values must be positive.")
    end
end

let psd_bins = get(cfg_toml, "PSDTB", [119, 4])
    global psd_lin_cos_bins, psd_log_θ_decs
    psd_lin_cos_bins::Int = psd_bins[1]
    psd_log_θ_decs::Int   = psd_bins[2]
    if psd_lin_cos_bins ≤ 0 || psd_log_θ_decs ≤ 0
        error("PSDTB: both values must be positive.")
    end
end

use_custom_frg = (get(cfg_toml, "NWFRG", 0) == 66)

emin_therm_fac = get(cfg_toml, "EMNFC", 0.01)

do_multi_dNdps = (get(cfg_toml, "DNDPS", 0) == 66)

begin
    do_tcuts = haskey(cfg_toml, "TCUTS")
    if do_tcuts
        tcuts = cfg_toml["TCUTS"]
        n_tcuts = length(tcuts)

        age_max < 0 && error("tcut tracking must be used with anaccel time limit. Adjust keyword 'AGEMX'.")
        # Check to make sure we haven't used more tcuts than allowed by na_c
        (n_tcuts+1) > na_c && error("TCUTS: parameter na_c smaller than desired number of tcuts.")
        # Check to make sure final tcut is much larger than age_max so that
        #   we never have to worry about exceeding it
        tcuts[end] ≤ 10age_max && error("TCUTS: final tcut must be much (10x) larger than age_max.")
    else
        tcuts = Float64[]
        n_tcuts = 0
    end
end

begin
    inj_fracs = get(cfg_toml,"INJFR", fill(1.0, n_ions))
    length(inj_fracs) == n_ions || error("Number of injection probabilities must match NIONS")
end

use_custom_εB = (get(cfg_toml, "NWEPB", 0) == 66)
