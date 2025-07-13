module MonteCarloScattering
import Random
using Dates
using JLD2
using Unitful, UnitfulAstro, UnitfulGaussian, UnitfulEquivalences
using Unitful: g, K, cm, s, dyn, erg, keV, GeV
using Unitful: mp, me, c, q, k as kB, h, ħ    # physical constants
using UnitfulGaussian: Fr, G, qcgs
using Cosmology
using Statistics: mean
using StaticArrays
using OffsetArrays
using JuliaInterpreter
using TOML

module CGSTypes
using Unitful: g, K, cm, s, dyn, erg, keV, GeV
using UnitfulGaussian: Fr, G, qcgs
export LengthCGS, TimeCGS, MomentumCGS, BFieldCGS, EnergyCGS, MomentumFluxCGS,
       MomentumDensityFluxCGS, EnergyFluxCGS, EnergyDensityFluxCGS
const LengthCGS              = typeof(1.0 * cm)
const TimeCGS                = typeof(1.0 * s)
const MomentumCGS            = typeof(1.0 * g*cm/s)
const BFieldCGS              = typeof(1.0 * G)
const EnergyCGS              = typeof(1.0 * erg)
const MomentumFluxCGS        = typeof(1.0 * g*cm^2/s^2) # (g*cm/s) * cm/s
const MomentumDensityFluxCGS = typeof(1.0 * erg/cm^3)
const EnergyFluxCGS          = typeof(1.0 * erg*cm/s)
const EnergyDensityFluxCGS   = typeof(1.0 * erg/(cm^2*s))
end
using .CGSTypes

include("parameters.jl"); using .parameters
include("constants.jl"); using .constants
include("utils.jl")
include("initializers.jl"); using .initializers
include("io.jl")
include("transformers.jl"); using .transformers
include("cosmo_calc.jl"); using .cosmo_calc
include("all_flux.jl")
include("prob_return.jl")
include("scattering.jl")
include("cuts.jl")
include("get_psd_bins.jl")
include("particle_finish.jl")
begin
    include("iter_init.jl")
    include("iter_finalize.jl")
    include("particle_loop.jl")
    include("ion_init.jl")
    include("ion_finalize.jl")
end

zero!(A::AbstractArray{T}) where T = fill!(A, zero(T))

include("particle_counter.jl")
include("smoothers.jl")
include("data_input.jl")

function (@main)()
    # Start the wall clock for this run
    t_start = now()

    # Get input, control variables, etc.
    cfg_toml = TOML.parsefile("mc_in.toml")

    (u₀, β₀, γ₀) = parse_shock_speed(cfg_toml["shock-speed"], cfg_toml["shock-speed-unit"])

    species = let
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
    n_ions = length(species)

    inp_distr = cfg_toml["input-distribution"]
    energy_inj = cfg_toml["ENINJ"] * keV
    inj_weight = get(cfg_toml, "INJWT", true)

    Emax, Emax_per_aa, pmax = parse_maximum_energy(cfg_toml["ENMAX"])

    η_mfp = get(cfg_toml, "gyrofactor", 1)

    bmag₀ = cfg_toml["BMAGZ"]*G
    # rg₀ below is the gyroradius of a proton whose speed is u₀ that is gyrating in a field
    # of strength bmag₀. Note that this formula is relativistically correct
    rg₀ = (γ₀ * E₀ₚ * β₀) / (qcgs * bmag₀) |> cm

    θ_B₀ = cfg_toml["THTBZ"] # must be zero
    check_shock_angle(θ_B₀)

    x_grid_start_rg, x_grid_stop_rg = cfg_toml["x_grid_limits"]
    check_x_grid_limits(x_grid_start_rg, x_grid_stop_rg)

    feb_upstream = let
        febup = get(cfg_toml, "FEBUP", nothing)
        if isnothing(febup)
            feb_upstream = x_grid_start_rg * rg₀ # default value
            return feb_upstream
        end
        if febup[1] < 0
            feb_upstream = febup[1] * rg₀
        elseif febup[2] < 0
            feb_upstream = uconvert(cm, febup[2] * pc)
        else
            error("FEBUP: at least one choice must be negative.")
        end
        (feb_upstream/rg₀ < x_grid_start_rg) && error("FEBUP: upstream FEB must be within x_grid_start")

        feb_upstream
    end

    feb_downstream, use_prp = let
        febdw = get(cfg_toml, "FEBDW", nothing)
        use_prp = false
        if isnothing(febdw)
            feb_downstream = -1 # default value
            return
        end
        if febdw[1] > 0
            feb_downstream = febdw[1] * rg₀
        elseif febdw[2] > 0
            feb_downstream = uconvert(cm, febdw[2] * pc)
        else
            feb_downstream = 0.0cm
            use_prp = true
        end
        (feb_downstream, use_prp)
    end

    begin
        x_spec = get(cfg_toml, "XSPEC", Float64[])
        n_xspec = length(x_spec)
    end

    n_itrs = cfg_toml["num-iterations"]
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

        if Emax > 0keV
            # Convert from momentum[mₚc/aa] to energy[keV]
            Emax_eff = 56 * pcuts_in[n_pcuts-1] * ustrip(keV, E₀ₚ*erg)

            if Emax > Emax_eff
                error("PCUTS: max energy exceeds highest pcut. Add more pcuts or lower Emax. ",
                      "Emax (assuming Fe) = $Emax; Emax_eff = $Emax_eff")
            end
        elseif Emax_per_aa > 0keV   # Limit was on energy per nucleon
            # Convert from momentum[mₚc/aa] to energy[keV/aa]
            Emax_eff = pcuts_in[n_pcuts-1] * ustrip(keV, E₀ₚ*erg)

            if Emax_per_aa > Emax_eff
                error("PCUTS: max energy per aa exceeds highest pcut. Add more pcuts or lower Emax_per_aa. ",
                      "Emax_per_aa = $Emax_per_aa; Emax_eff/aa = $Emax_eff")
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

    dont_shock = get(cfg_toml, "no-shock", false)
    dont_scatter = get(cfg_toml, "no-scatter", false)
    dont_DSA = get(cfg_toml, "no-DSA", false)
    do_smoothing = cfg_toml["smooth-shocks"]

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
            throw(DomainError(smooth_pressure_flux_psd_fac,
                              "smooth_pressure_flux_psd_fac/SMPFP must be in [0, 1]"))
        end
        # TODO: actually get pressure calculation working properly
        if smooth_pressure_flux_psd_fac > 0
            error("SMPFP: code does not properly calculate pressure from PSD. ",
                  "Set to 0 or get this code working")
        end
    end

    r_comp, r_RH, Γ₂_RH = let
        r_comp = cfg_toml["RCOMP"]
        r_RH, Γ₂_RH = calc_rRH(u₀, β₀, γ₀, species)
        if r_comp == -1
            r_comp = r_RH
        end
        (r_comp, r_RH, Γ₂_RH) # shadowed variables
    end

    begin
        β₂, γ₂, bmag₂, θ_B₂, θᵤ₂ = calc_downstream(bmag₀, r_comp, β₀)
        u₂ = β₂*c
        @debug("Results from calc_downstream()", u₂, β₂, γ₂, bmag₂, θ_B₂, θᵤ₂)
    end

    begin
        do_old_prof = (get(cfg_toml, "OLDIN", 0) == 66)
        if do_old_prof
            n_old_skip, n_old_profs, n_old_per_prof = cfg_toml["OLDDT"]
        else
            n_old_skip, n_old_profs, n_old_per_prof = [0, 0, 0]
        end
    end

    age_max = let
        age_max = get(cfg_toml, "AGEMX", -1.0)
        if age_max < 0
            age_max = -1.0
        end
        age_max * s
    end
    # default behavior of do_retro is dependent on age_max
    do_retro = (get(cfg_toml, "RETRO", age_max > 0s ? 66 : 0) == 66)

    do_fast_push = get(cfg_toml, "fast-upstream-transport", false)
    x_fast_stop_rg = do_fast_push ? cfg_toml["FPSTP"] : 0.0

    x_art_start_rg, x_art_scale = let
        art = get(cfg_toml, "ARTSM", nothing)
        if isnothing(art)
            x_art_start_rg = 0.0
            x_art_scale = 0.0
        else
            x_art_start_rg, x_art_scale = art
        end
        (x_art_start_rg, x_art_scale)
    end

    pₑ_crit, γₑ_crit = parse_electron_critical_energy(get(cfg_toml, "EMNFP", nothing))

    do_rad_losses = get(cfg_toml, "radiation-losses", true)
    do_photons = get(cfg_toml, "calculate-photon-production", false)

    # jet-shock-radius only mandatory if doing photons
    jet_rad_pc = do_photons ? cfg_toml["jet-shock-radius"] : get(cfg_toml, "jet-shock-radius", 0.0)

    jet_sph_frac, jet_open_ang_deg = let

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
        jet_dist_kpc = get(cfg_toml, "jet-distance", 1.0)
        redshift = get(cfg_toml, "RDSHF", 0.0)
        if jet_dist_kpc > 0 && redshift > 0
            error("jet-distance: At most one of 'jet-distance' and 'RDSHF' may be non-zero.")
        end
    end
    begin
        # The following option is not in the Fortran program
        cosmo_var = cfg_toml["COSMO_VAR"]
        cosmo_var ≠ 1 && cosmo_var ≠ 2 && error("Invalid value for cosmo_var")
    end

    begin
        energy_transfer_frac = float(get(cfg_toml, "ENXFR", 0.0))
        if energy_transfer_frac < 0 || energy_transfer_frac > 1
            error("ENXFR: energy_transfer_frac must be in [0,1]")
        end
    end

    num_upstream_shells, num_downstream_shells = cfg_toml["num-shells"]

    begin
        bturb_comp_frac = get(cfg_toml, "BTRBF", 0.0)
        bfield_amp = get(cfg_toml, "BAMPF", 1.0)
        bfield_amp < 1 && error("BAMPF: must be ≥ 1.d0")
        if bfield_amp > 1 && iszero(bturb_comp_frac)
            error("BTRBF: bfield_amp > 1 has no effect if BTRBF = 0")
        end
    end

    psd_bins_per_dec_mom, psd_bins_per_dec_θ = let
        psd_bins = get(cfg_toml, "PSDBD", [10, 10])
        psd_bins_per_dec_mom::Int = psd_bins[1]
        psd_bins_per_dec_θ::Int   = psd_bins[2]
        if psd_bins_per_dec_mom ≤ 0 || psd_bins_per_dec_θ ≤ 0
            error("PSDBD: both values must be positive.")
        end
        (psd_bins_per_dec_mom::Int, psd_bins_per_dec_θ::Int)
    end

    psd_lin_cos_bins, psd_log_θ_decs = let
        psd_bins = get(cfg_toml, "PSDTB", [119, 4])
        psd_lin_cos_bins = psd_bins[1]
        psd_log_θ_decs   = psd_bins[2]
        if psd_lin_cos_bins ≤ 0 || psd_log_θ_decs ≤ 0
            error("PSDTB: both values must be positive.")
        end
        (psd_lin_cos_bins::Int, psd_log_θ_decs::Int)
    end

    use_custom_frg = (get(cfg_toml, "NWFRG", 0) == 66)

    emin_therm_fac = get(cfg_toml, "EMNFC", 0.01)

    do_multi_dNdps = (get(cfg_toml, "DNDPS", 0) == 66)

    do_tcuts, tcuts, n_tcuts = let
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
        inj_fracs = get(cfg_toml, "INJFR", fill(1.0, n_ions))
        length(inj_fracs) == n_ions || error("Number of injection probabilities must match NIONS")
    end

    use_custom_εB = (get(cfg_toml, "NWEPB", 0) == 66)

    # Set up the computational grid
    x_grid_rg, x_grid_start, x_grid_stop = setup_grid(x_grid_start_rg, x_grid_stop_rg, use_prp, feb_downstream, rg₀)
    n_grid = length(x_grid_rg) - 2
    grid_axis = axes(x_grid_rg, 1)
    x_grid_cm = x_grid_rg * rg₀ # Convert everything from rg₀ units to cgs units
    Γ_grid = zeros(n_grid, 2)

    # debug variables
    energy_density = zeros(n_grid, n_ions)
    therm_energy_density = zeros(n_grid, n_ions)

    zone_vol = zeros(n_grid)
    energy_pool = zeros(n_grid)

    # Set quantities related to the phase space distribution, including the bins
    psd_cos_fine = 1 - 2 / (psd_lin_cos_bins+1)
    psd_θ_fine = acos(psd_cos_fine)
    psd_θ_min  = psd_θ_fine / 10^psd_log_θ_decs
    @debug("Phase space angle parameters", psd_cos_fine, psd_θ_fine, psd_θ_min)

    if inp_distr == 1
        # Set minimum PSD energy using thermal distribution for upstream plasma
        Emin = uconvert(keV, minimum(temperature.(species)), Thermal())

        # Allow for a few extra zones below the thermal peak
        Emin *= emin_therm_fac

    elseif inp_distr == 2
        # Set minimum PSD energy using δ-function dist for upstream plasma;
        # allow for a few extra zones below the location of the distribution
        Emin = energy_inj/5
    end

    # Determine minimum momentum associated with the given energy, which will occur for
    # the lightest particle species. Use a cutoff of 0.1% of the rest-mass energy for
    # relativistic/nonrelativistic calculation.
    psd_mom_min = let
        rest_mass_min = minimum(mass.(species))
        rest_energy_min = uconvert(erg, rest_mass_min, MassEnergy())
        if Emin < 1e-3*rest_energy_min
            √(2 * rest_mass_min * Emin)
        else
            γ = 1 + Emin/rest_energy_min
            rest_mass_min*c * √(γ^2 - 1)
        end
    end

    # Now find the maximum momentum for the PSD (this will be adjusted due to SF->PF Lorentz
    # transformation). How to actually calculate it depends on the user-specified maximum energy
    # "ENMAX"
    psd_mom_max = let
        rest_mass_max = maximum(mass.(species))
        rest_energy_max = uconvert(erg, rest_mass_max, MassEnergy())
        if Emax > 0keV
            γ = 1 + Emax/rest_energy_max
            rest_mass_max*c * √(γ^2 - 1)
        elseif Emax_per_aa > 0keV
            γ = 1 + Emax_per_aa/E₀ₚ
            rest_mass_max*c * √(γ^2 - 1)
        elseif pmax > 0g*cm/s
            pmax
        else
            # Something has gone very wrong.
            error("Max CR energy not set in data_input, so can not set PSD bins.")
        end
    end |> g*cm/s

    # Adjust max momentum based on a SF->PF Lorentz transform
    psd_mom_max *= 2γ₀
    num_psd_mom_bins, psd_mom_bounds = set_psd_mom_bins(psd_mom_min, psd_mom_max, psd_bins_per_dec_mom)
    psd_mom_axis = axes(psd_mom_bounds, 1)
    Δcos, psd_θ_bounds = set_psd_angle_bins(psd_bins_per_dec_θ, psd_lin_cos_bins, psd_cos_fine, psd_θ_min)
    num_psd_θ_bins = length(psd_θ_bounds)-2
    psd_θ_axis = axes(psd_θ_bounds, 1)
    @debug("Setting PSD parameters",
           psd_mom_max, num_psd_mom_bins, psd_mom_axis, psd_mom_bounds,
           Δcos, num_psd_θ_bins, psd_θ_axis, psd_θ_bounds)

    # Set the boundaries of the shells to use for photon calculation
    if do_photons
        x_shell_midpoints, x_shell_endpoints = set_photon_shells(num_upstream_shells, num_downstream_shells, use_prp,
                                                                 feb_upstream, feb_downstream, rg₀, x_grid_stop_rg)
    end

    begin # "module" iteration_vars
        # units of momentum flux density: [p] * [v] * [n] = [energy density]
        pxx_flux = Vector{MomentumDensityFluxCGS}(undef, n_grid)
        pxz_flux = Vector{MomentumDensityFluxCGS}(undef, n_grid)
        energy_flux = Vector{Float64}(undef, n_grid)
        esc_flux = zeros(n_ions)
        pₓ_esc_feb = zeros(n_ions, n_itrs)
        energy_esc_feb = zeros(n_ions, n_itrs)
        esc_energy_eff = zeros(0:psd_max, n_ions)
        esc_num_eff = zeros(0:psd_max, n_ions)

        # Arrays for holding thermal distribution information; they're set at the start of
        # the run, but included here because of chance they could change due to fast push
        n_pts_MB   = zeros(Int, n_ions)
        ptot_inj   = zeros(MomentumCGS, na_particles, n_ions)
        weight_inj = zeros(na_particles, n_ions)
        # Arrays for holding information about particle counts and spectra at various tcuts
        weight_coupled  = Matrix{Float64}(undef, na_c, n_ions)
        spectra_coupled = zeros(0:psd_max, na_c, n_ions)
    end # "module" iteration_vars

    outfile = open("mc_out.dat", "w")

    # Check x_spec data to make sure it falls within the grid, and add it to the output file
    if n_xspec > 0
        for i in 1:n_xspec
            if x_spec[i] < x_grid_start || x_spec[i] > x_grid_stop
                throw(DomainError("x_spec position $i falls outside grid start/stop bounds."))
            end

            println(outfile, "x position for spectrum calculation:  ",
                    i, x_spec[i]/rg₀, " rg₀ ", uconvert(pc, x_spec[i] * cm))
        end

        println(outfile)
    end


    # Get grid zone numbers for the boundaries between photon shells, and for
    # the location of the upstream FEB. Add the photon shells to the output file
    n_shell_endpoints = zeros(Int, num_upstream_shells+num_downstream_shells+1)
    if do_photons
        let i_tmp = 1
            for i in 1:n_grid
                if x_grid_cm[i] ≤ x_shell_endpoints[i_tmp] && x_grid_cm[i+1] > x_shell_endpoints[i_tmp]
                    n_shell_endpoints[i_tmp] = i
                    i_tmp += 1
                end
            end
        end

        for i in 1:(num_upstream_shells+num_downstream_shells)
            println(outfile, "x boundaries[rg₀] for photon shell ", i, ":",
                    x_shell_endpoints[i]/rg₀, " -> ", x_shell_endpoints[i+1]/rg₀, "  (",
                    n_shell_endpoints[i], "->", n_shell_endpoints[i+1], ")")
        end

        println(outfile)
    end

    i_grid_feb = findfirst(>(feb_upstream), x_grid_cm) - 1

    # Because redshift will be needed to compute radiative losses (it affects both the energy
    # and density of CMB photons), calculate it here. If redshift was provided during
    # data_input, cosmo_calc will return distance instead
    # Cosmo_calc expects distance in megaparsecs, so convert from value read in during data_input
    jet_dist_Mpc = jet_dist_kpc * 1e-3
    # TODO use cosmology.jl for this instead?
    #jet_dist_Mpc, redshift = cosmo_calc_wrapper(jet_dist_Mpc, redshift)
    if cosmo_var == 1
        # FIXME redshift is defined as a const in data_input.jl, so it should not be redefined here.
        # figure out scoping issues
        redshift = get_redshift(jet_dist_Mpc)
    else
        jet_dist_Mpc = ustrip(Mpc, comoving_radial_distance(cosmo_calc_params.cosmo, redshift))
    end


    # Set a handful of constants related to radiative losses. electron_rm will not be the same
    # as E₀ₑ in module "constants" unless (a) there are electrons in the run, and
    # (b) they are true electrons, with aa = mₑ/mₚ. (electron_rm may in fact = E₀ₚ,
    # but in that case it won't be used because radiative losses will never be calculated.)
    # To convert rad_loss_fac to dp/dt, multiply by p²B², both in cgs. Prefactor of rad_loss_fac
    # comes from average over pitch and is Eq (16) of Sturner+ (1997) [1997ApJ...490..619S].
    # Note the extra factor of c in the denominator, because code tracks dp/dt, not dE/dt as
    # given in Sturner+ (1997).
    # unit shenanigans because of overflow
    B_CMBz = B_CMB0 * (1 + redshift)^2
    @debug "Calculated radiation loss constants" B_CMBz


    # Zero out total escaping fluxes and calculate the far upstream fluxes
    #pₓ_esc_flux_upstream_tot = 0.0
    #energy_esc_flux_upstream_tot = 0.0
    pₓ_esc_flux_upstream     = zeros(n_itrs)
    energy_esc_flux_upstream = zeros(n_itrs)
    flux_px_upstream, flux_pz_upstream, flux_energy_upstream = upstream_fluxes(
        number_density.(species), temperature.(species), mass.(species),
        bmag₀, θ_B₀, u₀, β₀, γ₀)

    # Determine upstream Mach numbers (sonic & Alfvén)
    mach_sonic, mach_alfven = upstream_machs(β₀, species, bmag₀)


    # Set up the initial shock profile, or read it in from a file
    if ! do_old_prof
        (
         uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
         β_ef_grid, γ_ef_grid, btot_grid, θ_grid, εB_grid, bmag₂,
        ) = setup_profile(
                          u₀, β₀, γ₀, bmag₀, θ_B₀, r_comp, bturb_comp_frac, bfield_amp, use_custom_εB,
                          n_ions, species, flux_px_upstream, flux_energy_upstream, grid_axis,
                          x_grid_cm, x_grid_rg,
                         )
    else
        error("Reading old profiles not yet supported")
        (
         x_grid_rg, x_grid_cm, uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid, β_ef_grid, γ_ef_grid,
         btot_grid, εB_grid, θ_grid, n_grid, u₀, γ₀, rg₀, r_comp, r_RH, β₀, bmag₀,
         u₂, β₂, γ₂, θᵤ₂, bmag₂, θ_B₀, θ_B₂, flux_px_upstream, flux_pz_upstream, flux_energy_upstream
        ) = read_old_prof(n_old_skip, n_old_profs, n_old_per_prof)

        # Must set far upstream and downstream limits manually, since they won't be read in from the file
        x_grid_rg[begin] = -1e30
        x_grid_rg[end]   =  1e30
        x_grid_cm[begin] = -1e30 * rg₀
        x_grid_cm[end]   =  1e30 * rg₀
    end


    # Find the location of the shock
    i_shock = findlast(≤(0), x_grid_rg)
    isnothing(i_shock) && error("Shock location not found")


    # How frequently will the code print during initial particle propagation?
    n_print_pt = n_pts_inj < 500 ? 25 : 250


    # The random number seed depends on the specific particle.
    # Compute a necessary quantity to determine that seed
    n_pts_max = max(n_pts_pcut, n_pts_pcut_hi)


    # Because electrons might have different Monte Carlo weights than protons, set that ratio here
    # WARNING: assumes electrons are last ion species of input file
    electron_weight_fac = 1 / density(species[end])


    # Print a bunch of data about the run to screen/file
    print_input(n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,
                NaN, NaN, n_xspec, n_pcuts, n_grid, r_RH, r_comp,
                u₀, β₀, γ₀, u₂, β₂, γ₂, species, bmag₀, bmag₂, θ_B₀, θ_B₂, θᵤ₂,
                mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine,
                feb_upstream, feb_downstream, rg₀, age_max, energy_pcut_hi, do_fast_push, bturb_comp_frac,
                outfile)

    weights_file = jldopen("mc_coupled_weights.hdf5", "a+")
    spectra_file = jldopen("mc_coupled_spectra.hdf5", "a+")

    pressure_psd_par   = Vector{Float64}(undef, n_grid)
    pressure_psd_perp  = Vector{Float64}(undef, n_grid)
    energy_density_psd = Vector{Float64}(undef, n_grid)

    energy_transfer_pool = Vector{Float64}(undef, n_grid)
    energy_recv_pool     = Vector{Float64}(undef, n_grid)

    energy_density       = Matrix{Float64}(undef, na_c, n_ions)
    therm_energy_density = Matrix{Float64}(undef, na_c, n_ions)

    psd = OffsetArray{Float64}(undef, (psd_mom_axis, psd_θ_axis, n_grid))

    begin # "module" species_vars
        # Arrays will hold crossing data for thermal particles; pₓ and pt are shock frame values
        num_crossings = zeros(Int, n_grid)
        therm_grid = zeros(Int, na_cr)
        therm_pₓ_sk = zeros(MomentumCGS, na_cr)
        therm_ptot_sk = zeros(MomentumCGS, na_cr)
        therm_weight = zeros(typeof(1.0s/cm), na_cr)

        # Spectra at x_spec locations
        spectra_sf = zeros(0:psd_max, n_grid)
        spectra_pf = zeros(0:psd_max, n_grid)

        # Escaping spectra upstream and downstream from shock; 2-D arrays store angular information
        esc_spectra_feb_upstream = zeros(0:psd_max)
        esc_spectra_feb_downstream = zeros(0:psd_max)
        # should these be inverse velocity?
        esc_psd_feb_upstream = zeros(0:psd_max, 0:psd_max)
        esc_psd_feb_downstream = zeros(0:psd_max, 0:psd_max)
    end # "module" species_vars

    begin # "module" photons
        # Number of different sources for emission mechanism (i.e. ion species for
        # pion decay, photon fields for IC)
        #integer :: n_pion_specs, n_IC_specs

        # Energy bins for photon production by a specific mechanism to ensure
        # uniformity across different sources, in units of MeV
        energy_pion_MeV = zeros(na_photons)
        energy_IC_MeV = zeros(na_photons)

        # Arrays to hold summed spectra due to the various sources
        pion_photon_sum = zeros(na_photons, n_grid)
        ic_photon_sum = zeros(na_photons, n_grid)
    end # "module" photons

    weight_new = zeros(na_particles)
    ptot_pf_new = zeros(MomentumCGS, na_particles)
    pb_pf_new = zeros(MomentumCGS, na_particles)
    x_PT_cm_new = zeros(LengthCGS, na_particles)

    pcuts_use = zeros(MomentumCGS, na_c)

    begin # "module" pcut_vars
        l_save = zeros(Bool, na_particles) # Whether or not to save particle for next pcut
        grid_sav        = zeros(Int, na_particles)
        tcut_sav        = zeros(Int, na_particles)
        downstream_sav  = zeros(Bool, na_particles)
        inj_sav         = zeros(Bool, na_particles)
        weight_sav      = zeros(na_particles)
        ptot_pf_sav     = zeros(MomentumCGS, na_particles)
        pb_pf_sav       = zeros(MomentumCGS, na_particles)
        x_PT_cm_sav     = zeros(LengthCGS, na_particles)
        xn_per_sav      = zeros(na_particles)
        #zz_sav         = zeros(na_particles)
        prp_x_cm_sav    = zeros(LengthCGS, na_particles)
        acctime_sec_sav = zeros(TimeCGS, na_particles)
        φ_rad_sav       = zeros(na_particles)

        grid_new        = zeros(Int, na_particles)
        tcut_new        = zeros(Int, na_particles)
        downstream_new  = zeros(Bool, na_particles)
        inj_new         = zeros(Bool, na_particles)
        xn_per_new      = zeros(na_particles)
        #zz_new         = zeros(na_particles)
        prp_x_cm_new    = zeros(LengthCGS, na_particles)
        acctime_sec_new = zeros(TimeCGS, na_particles)
        φ_rad_new       = zeros(na_particles)
    end

    ε_target = zeros(n_grid)

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #  Main computational loops
    #
    #  Loops are set in the following sequence of nesting
    #   1. loop_itr:   i_iter   Iteration number
    #   2. loop_ion:   i_ion    Particle species number
    #   3. loop_pcut:  i_cut    Particle splitting (pcut) number
    #   4. loop_pt:    i_prt    Individual particle number
    #
    #
    # Start of loop over iterations
    #--------------------------------------------------------------------------
    for i_iter in 1:n_itrs # loop_itr

        @info("Starting loop", i_iter)

        # Zero out numerous quantities that will be modified over the course of this iteration.
        # Minimally positive number is used to prevent errors when taking logarithms later
        # XXX uses the type-pirated version of Base.fill!
        fill!((pxx_flux, pxz_flux), 1e-99erg/cm^3)
        fill!((energy_flux, pressure_psd_par, pressure_psd_perp, energy_density_psd), 1e-99)
        fill!((esc_spectra_feb_upstream, esc_spectra_feb_downstream), 1e-99)
        fill!(weight_coupled, 1e-99)


        # Additionally, set/reset scalar quantities that will change
        ∑P_downstream  = 1e-99erg/cm^3         # total downstream pressure
        ∑KEdensity_downstream = 1e-99erg/cm^3  # total downstream kinetic energy density

        energy_esc_upstream    = 1e-99
        pₓ_esc_upstream        = 1e-99


        # To facilitate energy transfer from ions to electrons, calculate here the target energy
        # density fraction for electrons at each grid zone, and zero out the pool of
        # plasma-frame energy that will be taken from ions and donated to electrons
        # Per Ardaneh+ (2015) [10.1088/0004-637X/811/1/57], ε_electron is proportional to √ε_b.
        # ε_b is itself roughly proportional to density² -- B² is proportional to z² (z being the
        # density compression factor) -- so ε_electron should vary roughly linearly with density.
        z_max = γ₀ * β₀ / (γ₂ * β₂)
        populate_ε_target!(ε_target, z_max, γ_sf_grid, uₓ_sk_grid, u₀, γ₀, energy_transfer_frac)

        fill!((energy_transfer_pool, energy_recv_pool), 0.0)
        fill!((energy_density, therm_energy_density), 0.0)

        #------------------------------------------------------------------------
        #  Start of loop over particle species
        #
        #  First species always protons. If electrons present they MUST be last species.
        #  Each species has mass "aa" in units of proton mass and charge "zz" in units of
        #  elementary charge.
        #------------------------------------------------------------------------
        for i_ion in 1:n_ions # loop_ion

            zz = charge(species[i_ion])
            m = mass(species[i_ion])
            aa = m/mp |> NoUnits
            mc = m*c

            # At the start of each ion, print a glyph to the screen
            @info("Starting species iteration", i_iter, i_ion)

            pmax_cutoff = get_pmax_cutoff(Emax, Emax_per_aa, pmax)

            # Zero out the phase space distributions and set variables related to
            # tracking thermal particles
            n_cr_count = clear_psd!(num_crossings, therm_grid, therm_pₓ_sk, therm_ptot_sk,
                                    therm_weight, psd, esc_psd_feb_upstream, esc_psd_feb_downstream)

            # In addition to initializing the phase space distribution, open the scratch (i.e.
            # temporary) file to which we will write information about thermal particle grid crossings
            nc_unit = open("mc_crossings.dat", "a")

            # To maintain identical results between OpenMP and serial versions,
            # set RNG seed based on current iteration/ion/pcut/particle number
            iseed_mod = (i_iter - 1)*n_ions + (i_ion - 1)
            Random.seed!(iseed_mod)


            # Initialize the particle populations that will be propagated through
            # the shock structure
            (n_pts_use, i_grid_in, weight_in, ptot_pf_in, pb_pf_in, x_PT_cm_in,
             pxx_flux, pxz_flux, energy_flux) = init_pop(
                do_fast_push, inp_distr, i_ion, m,
                temperature.(species), energy_inj, inj_weight, n_pts_inj,
                density.(species), x_grid_start, rg₀, η_mfp, x_fast_stop_rg,
                β₀, γ₀, u₀, n_ions, mass.(species),
                n_grid, x_grid_rg, uₓ_sk_grid, γ_sf_grid,
                ptot_inj, weight_inj, n_pts_MB,
            )
            @info("Finished init_pop on",
                  i_iter, i_ion,
                  n_pts_use, i_grid_in, weight_in, ptot_pf_in,
                  pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux)

            # Assign the various particle properties to the population
            assign_particle_properties_to_population!(
                n_pts_use, xn_per_fine, x_grid_stop,
                weight_new, weight_in, ptot_pf_new, ptot_pf_in,
                pb_pf_new, pb_pf_in, x_PT_cm_new, x_PT_cm_in, grid_new, i_grid_in,
                downstream_new, inj_new, xn_per_new, prp_x_cm_new, acctime_sec_new, tcut_new, φ_rad_new)

            # Weight of remaining particles, printed after each pcut; note that this will not be
            # correct for all particles if they were originally created so each thermal bin
            # would have equal weight
            weight_running = weight_in[1]

            # When using OpenMP, the array energy_transfer_pool can't be conditionally assigned
            # shared or reduction status, so it can't be used for both the ion and electron loops.
            # To get around this, use one array to hold the donated energy, and another to hold
            # the received energy.
            energy_recv_pool .= energy_transfer_pool


            # The array of pcuts read in by data_input has units momentum/mc.
            # Convert to momentum for this species
            pcuts_use[1:n_pcuts] .= pcuts_in[1:n_pcuts] * aa*mp*c

            # Determine the pcut at which to switch from low-E particle counts to high-E
            # particle counts. Recall that energy_pcut_hi has units of keV per aa, so when
            # dividing by particle mass the factor of aa is already present in the denominator.
            # Also set the maximum momentum cutoff based on the values given in keyword "ENMAX"
            p_pcut_hi = pcut_hi(energy_pcut_hi, energy_rel_pt, mass(species[i_ion]))

            #----------------------------------------------------------------------
            #  Start of loop over pcuts
            #
            #  Particle splitting and resetting of n_pts_use is handled by the call
            #  to "new_pcut" at the end of each iteration of loop_pcut.
            #----------------------------------------------------------------------
            for i_cut in 1:n_pcuts # loop_pcut

                # Initialize all of the *_sav arrays to help prevent bleeding over between pcuts or ion species
                l_save .= false  # Whole array must be initialized in case number of particles changes from pcut to pcut

                zero!(weight_sav)
                zero!(ptot_pf_sav)
                zero!(pb_pf_sav)
                zero!(x_PT_cm_sav)
                zero!(grid_sav)
                zero!(downstream_sav)
                zero!(inj_sav)
                zero!(xn_per_sav)
                zero!(prp_x_cm_sav)
                zero!(acctime_sec_sav)
                zero!(φ_rad_sav)
                zero!(tcut_sav)

                # A separate variable tracks the number of finished particles,
                # so that race conditions can be avoided in OMP mode
                i_fin = 0

                # For high-energy electrons in a strong magnetic field, need to know
                # previous cutoff momentum for calculating new PRP downstream
                pcut_prev = i_cut > 1 ? pcuts_use[i_cut-1] : 0.0

                #--------------------------------------------------------------------
                #  Start of loop over particles
                #
                #  Quick explanation of particle properties. All units cgs where a unit exists.
                #
                #  weight: weighting value (used in momentum splitting)
                #  ptot_pf: total plasma frame momentum
                #  pb_pf: momentum along B field in plasma frame. NOT along x-axis unless
                #         upstream orientation is parallel
                #  r_PT_cm: current particle position
                #  i_grid: current grid zone number
                #  l_downstream: whether particle has been downstream
                #  inj: whether particle has been back upstream (i.e. is a CR)
                #  xn_per: time steps per gyro period, Δt = T_g/xn_per
                #  prp_x_cm: downstream position of PRP; adjusted to allow all particles,
                #            regardless of momentum, to isotropize before reaching
                #  acctime_sec: time since crossing shock for first time; not started until l_downstream = true
                #  φ_rad: phase angle of gyration relative to z axis
                #  tcut_curr: next time at which particle tracking takes place
                #--------------------------------------------------------------------
                #$omp parallel for default(none), schedule(dynamic,1), num_threads(6)
                for i_prt in 1:n_pts_use # loop_pt

                    (i_fin, ∑P_downstream, ∑KEdensity_downstream) = particle_loop(
                        i_iter, i_ion, i_cut, i_prt, i_grid_feb, i_shock,
                        n_ions, n_pcuts, n_pts_max, n_xspec, n_grid,
                        γ₀, β₀, u₀, u₂, bmag₂,
                        pₑ_crit, γₑ_crit, η_mfp,
                        energy_transfer_frac, energy_recv_pool,
                        psd, num_crossings, x_spec, feb_upstream, feb_downstream,
                        energy_esc_upstream, pₓ_esc_upstream,
                        pcut_prev, i_fin, ∑P_downstream, ∑KEdensity_downstream,
                        aa, zz, m, mc,
                        nc_unit, n_cr_count, pmax_cutoff,
                        B_CMBz,
                        weight_new, ptot_pf_new, pb_pf_new, grid_new, downstream_new, inj_new,
                        xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new, tcut_new, x_PT_cm_new,
                        use_custom_εB, x_grid_stop,
                        uₓ_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
                        γ_ef_grid, β_ef_grid, btot_grid, θ_grid,
                        pxx_flux, pxz_flux, energy_flux,
                        do_rad_losses, do_retro, do_tcuts, dont_DSA, dont_scatter, use_custom_frg,
                        inj_fracs, xn_per_fine, xn_per_coarse,
                        x_grid_cm, therm_grid, therm_ptot_sk, therm_pₓ_sk, therm_weight,
                        spectra_pf, spectra_sf, tcuts, age_max,
                       )

                end # loop_pt
                #$omp end parallel do
                #--------------------------------------------------------------------
                # Conclusion of particle loop

                break_pcut = pcut_finalize(i_iter, i_ion, i_cut, p_pcut_hi, n_pts_pcut, n_pts_pcut_hi)
                break_pcut && break


            end # loop_pcut
            #----------------------------------------------------------------------
            # Conclusion of pcuts loop

            ion_finalize()

        end # loop_ion
        #------------------------------------------------------------------------
        # Conclusion of species loop

        #DEBUGLINE
        for i in 1:n_grid
            @debug("",
                   i_iter, i,
                   #x_grid_rg[i],
                   therm_energy_density[i,1:n_ions]/zone_vol[i],
                   therm_energy_density[i,n_ions]/max(sum(therm_energy_density[i,1:n_ions]),1e-99),
                   sum(therm_energy_density[i,1:n_ions]/zone_vol[i]),
                   energy_density[i,1:n_ions]/zone_vol[i],
                   energy_density[i,n_ions]/max(sum(energy_density[i,1:n_ions]),1e-99),
                   sum(energy_density[i,1:n_ions]/zone_vol[i]))
        end
        #print_plot_vals(888)

        iter_finalize()

    end # loop_itr
    #--------------------------------------------------------------------------
    # Conclusion of iteration loop

    close(weights_file)
    close(spectra_file)

    # End the run
    t_end = now()

    run_time = round(t_end - t_start, Second)

    println()
    @info(" Finished. Run time = $(run_time) sec, $(round(run_time, Minute)) min")

    println(outfile)
    println(outfile, " Finished. Run time = ", run_time, ", ", round(run_time, Minute))
    println(outfile)
end # main
end # module
# vim: set textwidth=92:shiftwidth=4:
