module MonteCarloScattering
import Random
using Dates
using JLD2
using CSV
using Unitful, UnitfulAstro, UnitfulGaussian, UnitfulEquivalences
using Unitful: g, K, km, cm, s, dyn, erg, keV, GeV
using Unitful: mp, me, c, q, k as kB, h, ħ    # physical constants
using UnitfulAstro: Mpc
using UnitfulGaussian: Fr, G, qcgs
using Cosmology
using Statistics: mean
using StaticArrays
using OffsetArrays
using JuliaInterpreter
using TOML
using Logging

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
const EnergyDensityCGS       = typeof(1.0 * erg/cm^3)
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
include("identify_corners.jl")
include("transformers.jl"); using .transformers
include("cosmo_calc.jl"); using .cosmo_calc
include("smoothers.jl"); using .smoothers
include("all_flux.jl")
include("prob_return.jl")
include("scattering.jl")
include("cuts.jl")
include("get_psd_bins.jl")
include("particle_finish.jl")
include("q_esc_calcs.jl")
include("particle_counter.jl")

include("iter_init.jl")
include("iter_finalize.jl")
include("particle_loop.jl")
include("ion_init.jl")
include("ion_finalize.jl")

# to shut the linter up
using .parameters: na_c, na_photons, psd_max
using .constants: B_CMB0, E₀ₑ
using .initializers: set_photon_shells
using .particle_counter: get_normalized_dNdp

zero!(A::AbstractArray{T}) where T = fill!(A, zero(T))

include("data_input.jl")
include("main_loops.jl")

function @main(args)
    # Start the wall clock for this run
    t_start = now()

    global_logger(ConsoleLogger(show_limited=false))

    # Get input, control variables, etc.
    cfg_toml = TOML.parsefile("mc_in.toml")

    (u₀, β₀, γ₀) = parse_shock_speed(cfg_toml["shock-speed"], cfg_toml["shock-speed-unit"])

    species = parse_species(cfg_toml)
    n_ions = length(species)

    inp_distr = cfg_toml["input-distribution"]
    energy_inj = cfg_toml["injection-energy"] * keV
    inj_weight = get(cfg_toml, "injection-weights", true)

    Emax, Emax_per_aa, pmax = parse_maximum_energy(cfg_toml["maximum-energy"])

    η_mfp = get(cfg_toml, "gyrofactor", 1)

    bmag₀ = cfg_toml["B-mag-upstream"]*G
    # rg₀ below is the gyroradius of a proton whose speed is u₀ that is gyrating in a field
    # of strength bmag₀. Note that this formula is relativistically correct
    rg₀ = (γ₀ * E₀ₚ * β₀) / (qcgs * bmag₀) |> cm

    θ_B₀ = cfg_toml["theta-B0"] # must be zero
    check_shock_angle(θ_B₀)

    x_grid_start_rg, x_grid_stop_rg = cfg_toml["x_grid_limits"]
    check_x_grid_limits(x_grid_start_rg, x_grid_stop_rg)

    feb_upstream, feb_downstream, use_prp = get_feb(
            get(cfg_toml, "FEB-upstream", nothing),
            get(cfg_toml, "FEB-downstream", nothing),
            x_grid_start_rg, rg₀)

    x_spec = get(cfg_toml, "XSPEC", Float64[])
    n_xspec = length(x_spec)

    n_itrs = cfg_toml["num-iterations"]
    xn_per_coarse = cfg_toml["coarse-scattering-Ng"]
    xn_per_fine = cfg_toml["fine-scattering-Ng"]

    n_pts_inj = cfg_toml["N_PTS_INJ"]
    n_pts_pcut = cfg_toml["N_PTS_PCUT"]
    max(n_pts_inj,n_pts_pcut) > na_particles && error("Array size na_particles too small.")

    n_pts_pcut_hi = cfg_toml["N_PTS_PCUT_HI"]
    energy_pcut_hi = cfg_toml["EN_PCUT_HI"]
    n_pts_pcut_hi > na_particles && error("Array size na_particles too small.")

    pcuts = cfg_toml["momentum-cutoffs"] * mp * c .|> (g*cm/s)
    check_pcuts(pcuts, Emax, Emax_per_aa, pmax)

    dont_shock = get(cfg_toml, "no-shock", false)
    dont_scatter = get(cfg_toml, "no-scatter", false)
    dont_DSA = get(cfg_toml, "no-DSA", false)
    do_smoothing = cfg_toml["smooth-shocks"]

    prof_weight_fac = get(cfg_toml, "old-profile-weight", 1.0)

    do_prof_fac_damp = get(cfg_toml, "increase-old-profile-weighting", false)

    smooth_mom_energy_fac = get(cfg_toml, "SMMOE", 0.0)
    if smooth_mom_energy_fac < 0 || smooth_mom_energy_fac > 1
        throw(DomainError(smooth_mom_energy_fac, "smooth_mom_energy_fac/SMMOE must be in [0, 1]"))
    end

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

    r_comp, r_RH, Γ₂_RH = let
        r_comp = cfg_toml["target-compression-ratio"]
        r_RH, Γ₂_RH = calc_rRH((u₀, β₀, γ₀), species)
        if r_comp == -1
            r_comp = r_RH
        end
        (r_comp, r_RH, Γ₂_RH) # shadowed variables
    end

    β₂, γ₂, bmag₂, θ_B₂, θᵤ₂ = calc_downstream(bmag₀, r_comp, β₀)
    u₂ = β₂*c |> cm/s
    @debug("Results from calc_downstream()", u₂, β₂, γ₂, bmag₂, θ_B₂, θᵤ₂)


    do_old_prof = get(cfg_toml, "read-old-profile", false)
    if do_old_prof
        d = cfg_toml["old-profile-config"]
        n_old_skip = d["lines-to-skip"]
        n_old_profs = d["profiles-to-average"]
        n_old_per_prof = d["lines-per-profile"]
    else
        n_old_skip, n_old_profs, n_old_per_prof = [0, 0, 0]
    end

    age_max = let
        age_max = get(cfg_toml, "maximum-age", -1.0)
        if age_max < 0
            age_max = -1.0
        end
        age_max * s
    end
    # default behavior of do_retro is dependent on age_max
    do_retro = get(cfg_toml, "use-retro", age_max > 0s ? true : false)

    do_fast_push = get(cfg_toml, "fast-upstream-transport", false)
    x_fast_stop_rg = do_fast_push ? cfg_toml["proton-fast-transport-stop"] : 0.0

    x_art_start_rg, x_art_scale = get(cfg_toml, "artificial-smoothing", (0.0, 0.0))

    pₑ_crit, γₑ_crit = parse_electron_critical_energy(get(cfg_toml, "electron-energy-mfp-threshold", nothing))

    do_rad_losses = get(cfg_toml, "radiation-losses", true)
    do_photons = get(cfg_toml, "calculate-photon-production", false)

    # jet-shock-radius only mandatory if doing photons
    jet_rad_pc = do_photons ? cfg_toml["jet-shock-radius"] : get(cfg_toml, "jet-shock-radius", 0.0)

    jet_sph_frac, jet_open_ang_deg = parse_jet_frac(get(cfg_toml, "JETFR", nothing), do_photons)

    jet_dist_kpc = get(cfg_toml, "jet-distance", 1.0)
    redshift = get(cfg_toml, "RDSHF", 0.0)
    if jet_dist_kpc > 0 && redshift > 0
        error("jet-distance: At most one of 'jet-distance' and 'RDSHF' may be non-zero.")
    end

    # The following option is not in the Fortran program
    cosmo_var = cfg_toml["COSMO_VAR"]
    cosmo_var ≠ 1 && cosmo_var ≠ 2 && error("Invalid value for cosmo_var")

    energy_transfer_frac = float(get(cfg_toml, "energy-transfer-frac", 0.0))
    0 ≤ energy_transfer_frac ≤ 1 || error("energy_transfer_frac must be in [0,1]")

    num_upstream_shells, num_downstream_shells = cfg_toml["num-shells"]

    bturb_comp_frac = get(cfg_toml, "b-field-turbulence", 0.0)
    bfield_amp = get(cfg_toml, "b-field-amplify", 1.0)
    bfield_amp < 1 && error("b-field-amplify: must be ≥ 1")
    if bfield_amp > 1 && iszero(bturb_comp_frac)
        error("b-field-turbulence: bfield_amp > 1 has no effect if b-field-turbulence = 0")
    end

    psd_bins_per_dec_mom, psd_bins_per_dec_θ = let
        psd_bins = get(cfg_toml, "num-psd-bins-per-decade", [10, 10])
        psd_bins_per_dec_mom::Int = psd_bins[1]
        psd_bins_per_dec_θ::Int   = psd_bins[2]
        if psd_bins_per_dec_mom ≤ 0 || psd_bins_per_dec_θ ≤ 0
            error("num-psd-bins-per-decade: both values must be positive.")
        end
        (psd_bins_per_dec_mom::Int, psd_bins_per_dec_θ::Int)
    end

    psd_lin_cos_bins = get(cfg_toml, "psd-linear-cosine-bins", 119)
    psd_lin_cos_bins > 0 || error("psd-linear-cosine-bins must be positive")

    psd_log_θ_decs = get(cfg_toml, "psd-log-theta-decs", 4)
    psd_log_θ_decs > 0 || error("psd-log-theta-decs must be positive")

    use_custom_frg = get(cfg_toml, "use-custom-frg", false)

    emin_therm_fac = get(cfg_toml, "EMNFC", 0.01)

    do_multi_dNdps = get(cfg_toml, "separate-dNdp-write", false)

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

    inj_fracs = get(cfg_toml, "INJFR", fill(1.0, n_ions))
    length(inj_fracs) == n_ions || error("Number of injection probabilities must match NIONS")

    use_custom_εB = get(cfg_toml, "use-custom-epsB", false)

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
    psd_θ_min  = psd_θ_fine / exp10(psd_log_θ_decs)
    @debug("Phase space angle parameters", psd_cos_fine, psd_θ_fine, psd_θ_min)

    if inp_distr == 1
        # Set minimum PSD energy using thermal distribution for upstream plasma
        Emin = uconvert(keV, minimum(temperature.(species)), Thermal())
        Emin *= emin_therm_fac # Allow for a few extra zones below the thermal peak
    elseif inp_distr == 2
        # Set minimum PSD energy using δ-function dist for upstream plasma;
        # allow for a few extra zones below the location of the distribution
        Emin = energy_inj/5
    else
        error("Unknown input distribution ", inp_distr)
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
    # "maximum-energy"
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
    @debug("Setting PSD momentum parameters",
           psd_mom_max, num_psd_mom_bins, psd_mom_axis, psd_mom_bounds)
    Δcos, psd_θ_bounds = set_psd_angle_bins(psd_bins_per_dec_θ, psd_lin_cos_bins, psd_cos_fine, psd_θ_min)
    num_psd_θ_bins = length(psd_θ_bounds)-2
    psd_θ_axis = axes(psd_θ_bounds, 1)
    @debug("Setting PSD angle parameters",
           Δcos, num_psd_θ_bins, psd_θ_axis, psd_θ_bounds)

    # Set the boundaries of the shells to use for photon calculation
    if do_photons
        x_shell_midpoints, x_shell_endpoints = set_photon_shells(num_upstream_shells, num_downstream_shells, use_prp,
                                                                 feb_upstream, feb_downstream, rg₀, x_grid_stop_rg)
    else
        x_shell_midpoints, x_shell_endpoints = nothing, nothing
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
        jet_dist_Mpc = ustrip(Mpc, comoving_radial_dist(cosmo_calc.cosmo, redshift))
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
                num_psd_mom_bins, num_psd_θ_bins, n_xspec, length(pcuts), n_grid, r_RH, r_comp,
                u₀, β₀, γ₀, u₂, β₂, γ₂, species, bmag₀, bmag₂, θ_B₀, θ_B₂, θᵤ₂,
                mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine,
                feb_upstream, feb_downstream, rg₀, age_max, energy_pcut_hi, do_fast_push, bturb_comp_frac)

    weights_file = "mc_coupled_weights.csv"
    spectra_file = jldopen("mc_coupled_spectra.hdf5", "a+")

    pressure_psd_par   = Vector{Float64}(undef, n_grid)
    pressure_psd_perp  = Vector{Float64}(undef, n_grid)
    energy_density_psd = Vector{Float64}(undef, n_grid)

    energy_transfer_pool = Vector{EnergyCGS}(undef, n_grid)
    energy_recv_pool     = Vector{EnergyCGS}(undef, n_grid)

    energy_density       = Matrix{Float64}(undef, na_c, n_ions)
    therm_energy_density = Matrix{Float64}(undef, na_c, n_ions)

    psd = OffsetArray{Float64}(undef, (psd_mom_axis, psd_θ_axis, n_grid))

    # "module" species_vars
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
    # end "module" species_vars

    # "module" photons
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
    # end "module" photons

    weight_new = zeros(na_particles)
    ptot_pf_new = zeros(MomentumCGS, na_particles)
    pb_pf_new = zeros(MomentumCGS, na_particles)
    x_PT_cm_new = zeros(LengthCGS, na_particles)

    # "module" pcut_vars
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
    # end "module" pcut_vars

    ε_target = zeros(n_grid)
    q_esc_cal_pₓ = zeros(n_itrs)
    q_esc_cal_energy = zeros(n_itrs)

    # Adiabatic index of particles that were lost downstream
    Γ_downstream = zeros(n_itrs)

    m_ion = mass.(species)
    aa_ion = m_ion / mp .|> NoUnits
    zz_ion = charge.(species)
    T₀_ion = temperature.(species)
    n₀_ion = density.(species)
    main_loops(
        n_itrs, n_ions, n_grid, n_pts_inj, n_tcuts, species,
        (u₀, β₀, γ₀), (u₂, β₂, γ₂),
        Emax, Emax_per_aa, energy_pcut_hi, pmax,
        pxx_flux, pxz_flux, energy_flux,
        pressure_psd_par, pressure_psd_perp, energy_density_psd,
        esc_spectra_feb_upstream, esc_spectra_feb_downstream,
        weight_coupled,
        ε_target, Γ₂_RH, εB_grid, Γ_grid, γ_sf_grid, uₓ_sk_grid, uz_sk_grid, utot_grid, energy_transfer_frac,
        energy_transfer_pool, energy_recv_pool, energy_density, therm_energy_density,
        num_crossings, therm_grid, therm_pₓ_sk, therm_ptot_sk,
        therm_weight, psd, esc_psd_feb_upstream, esc_psd_feb_downstream,
        do_fast_push, inp_distr, energy_inj, inj_weight, ptot_inj, weight_inj,
        x_grid_start, x_grid_stop, rg₀, η_mfp, x_fast_stop_rg, x_grid_rg, x_grid_cm, n_pts_MB,
        xn_per_fine, xn_per_coarse, xn_per_sav, xn_per_new,
        weight_sav, weight_new, ptot_pf_sav, ptot_pf_new, pb_pf_sav, pb_pf_new,
        x_PT_cm_sav, x_PT_cm_new, grid_sav, grid_new, inj_sav, inj_new, downstream_sav, downstream_new,
        prp_x_cm_sav, prp_x_cm_new, acctime_sec_sav, acctime_sec_new, tcut_sav, tcut_new, φ_rad_sav, φ_rad_new,
        pcuts, l_save, i_grid_feb, i_shock,
        n_pts_max, n_xspec, n_print_pt, num_psd_mom_bins, num_psd_θ_bins,
        psd_lin_cos_bins, psd_cos_fine, Δcos, psd_mom_bounds, psd_θ_bounds, psd_θ_min,
        psd_mom_min, psd_bins_per_dec_mom, psd_bins_per_dec_θ,
        bmag₂, zone_vol, pₑ_crit, γₑ_crit,
        x_spec, feb_upstream, feb_downstream, B_CMBz, use_custom_εB,
        γ_ef_grid, β_ef_grid, btot_grid, θ_grid, pₓ_esc_feb, energy_esc_feb,
        do_rad_losses, do_retro, do_tcuts, dont_DSA, dont_scatter, use_custom_frg,
        inj_fracs, spectra_pf, spectra_sf, tcuts, age_max, spectra_coupled,
        esc_energy_eff, esc_num_eff, esc_flux, electron_weight_fac,
        n_pts_pcut, n_pts_pcut_hi, t_start, weights_file, spectra_file, outfile,
        pₓ_esc_flux_upstream, flux_px_upstream, flux_energy_upstream, energy_esc_flux_upstream,
        Γ_downstream, q_esc_cal_pₓ, q_esc_cal_energy,
        do_multi_dNdps, do_photons,
        jet_rad_pc, jet_sph_frac, m_ion, aa_ion, zz_ion, T₀_ion, n₀_ion,
        r_comp, r_RH, n_shell_endpoints,
    )

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
