import Random
using Dates
using JLD2
using Unitful, UnitfulAstro
using Cosmology
using Statistics: mean
using StaticArrays
using OffsetArrays

include("parameters.jl"); using .parameters
include("constants.jl"); using .constants
include("initializers.jl"); using .initializers
include("io.jl"); using .io
include("transformers.jl"); using .transformers
include("cosmo_calc.jl"); using .cosmo_calc
include("all_flux.jl")
include("prob_return.jl")
include("scattering.jl")
include("new_pcut.jl")
include("tcut_track.jl")
include("get_psd_bins.jl")
include("particle_finish.jl")
include("helix_loop.jl")

# Start the wall clock for this run
t_start = now()

γ_adiab_grid = zeros(na_grid, 2)

#outfileunit = Ref{Int}()
#outfd = Ref{Int}()
#@ccall "./particle_scattering.so".__io_MOD_openergy_log_file(outfileunit::Ref{Int}, outfd::Ref{Int})::Cvoid
#outfile = fdio(outfd[])
outfile = open("mc_out.dat", "w")

# Get input, control variables, etc.
include("data_input.jl") # FIXME

# Set quantities related to the phase space distribution, including the bins
psd_cos_fine = 1  -  2 / (psd_lin_cos_bins+1)
psd_θ_fine = acos(psd_cos_fine)
psd_θ_min  = psd_θ_fine / 10^psd_log_θ_decs

if inp_distr == 1
    # Set minimum PSD energy using thermal distribution for UpS plasma
    Emin_keV = kB_cgs * minimum(tZ_ion[1:n_ions]) * erg2keV

    # Allow for a few extra zones below the thermal peak
    Emin_keV = emin_therm_fac * Emin_keV

elseif inp_distr == 2
    # Set minimum PSD energy using δ-function dist for UpS plasma;
    # allow for a few extra zones below the location of the distribution
    Emin_keV = 0.2 * energy_inj
end

# Determine minimum momentum associated with the given energy, which will
# occur for the lightest particle species. Use a cutoff of 0.1% of the
# rest-mass energy for rel/non-rel calculation.
aa_min = minimum(aa_ion)
if Emin_keV*keV2erg < 1e-3*aa_min*E₀_proton
    psd_mom_min = √( 2 * aa_min*mp_cgs * Emin_keV * keV2erg )
else
    γ           = 1  +  (Emin_keV*keV2erg)/(aa_min*E₀_proton)
    psd_mom_min = aa_min * mp_cgs * c_cgs * √( γ^2 - 1 )
end

# Now find the maximum momentum for the PSD (this will be adjusted due to
# SF->PF Lorentz transformation). How to actually calculate it depends
# on the user-specified maximum energy "ENMAX"
aa_max = maximum(aa_ion)
if Emax_keV > 0
    _γ          = 1  +  (Emax_keV*keV2erg)/(aa_max*E₀_proton)
    psd_mom_max = aa_max * mp_cgs * c_cgs * √( _γ^2 - 1 )
elseif Emax_keV_per_aa > 0
    _γ          = 1  +  (Emax_keV_per_aa*keV2erg/E₀_proton)
    psd_mom_max = aa_max * mp_cgs * c_cgs * √( _γ^2 - 1 )
elseif pmax_cgs > 0
    psd_mom_max = pmax_cgs
else
    # Something has gone very wrong.
    error("Max CR energy not set in data_input, so can not set PSD bins.")
end

# Adjust max momentum based on a SF->PF Lorentz transform
psd_mom_max *= 2γ_Z
num_psd_mom_bins, psd_mom_bounds = set_psd_mom_bins(
    psd_mom_min, psd_mom_max, psd_bins_per_dec_mom)
num_psd_θ_bins, Δcos, psd_θ_bounds = set_psd_angle_bins(
    psd_bins_per_dec_θ, psd_lin_cos_bins, psd_cos_fine, psd_θ_min) # Set angle bins

# Set the boundaries of the shells to use for photon calculation
if do_photons
    x_shell, x_shell_end_points = set_photon_shells(
        num_UpS_shells, num_DwS_shells, use_prp,
        feb_UpS, feb_DwS, rg0, x_grid_stop_rg)
end

begin # "module" iteration_vars
    pxx_flux = Vector{Float64}(undef, na_grid)
    pxz_flux = Vector{Float64}(undef, na_grid)
    energy_flux  = Vector{Float64}(undef, na_grid)
    esc_flux = zeros(na_ions)
    px_esc_feb = zeros(na_ions, na_itrs)
    energy_esc_feb = zeros(na_ions, na_itrs)
    esc_energy_eff = zeros(0:psd_max, na_ions)
    esc_num_eff = zeros(0:psd_max, na_ions)
    # Arrays for holding thermal distribution information; they're set at the start of
    # the run, but included here because of chance they could change due to fast push
    n_pts_MB = zeros(Int, na_ions)
    ptot_inj = zeros(na_particles, na_ions)
    wt_inj   = zeros(na_particles, na_ions)
    # Arrays for holding information about particle counts and spectra at various tcuts
    wt_coupled   = Matrix{Float64}(undef, na_c, na_ions)
    spectra_coupled = zeros(0:psd_max, na_c, na_ions)
end # "module" iteration_vars


# Set up the computational grid
(n_grid, x_grid_start, x_grid_stop, x_grid_rg, x_grid_cm) = setup_grid(
    outfile, x_grid_start_rg, use_prp, feb_DwS, x_grid_stop_rg, rg0)


# Check x_spec data to make sure it falls within the grid, and add it to
# the output file
if n_xspec > 0
    for i in 1:n_xspec
        if x_spec[i] < x_grid_start || x_spec[i] > x_grid_stop
            throw(DomainError("x_spec position $i falls outside grid start/stop bounds."))
        end

        println(outfile, "x position[rg0,pc] for spectrum calculation:  ", i, x_spec[i]/rg0, x_spec[i]/pc2cm)
    end

    println(outfile)
end


# Get grid zone numbers for the boundaries between photon shells, and for
# the location of the UpS FEB. Add the photon shells to the output file
n_shell_end_points = zeros(Int, na_grid)
if do_photons
    local i_tmp = 1
    for i in 1:n_grid
        if x_grid_cm[i] ≤ x_shell_end_points[i_tmp] && x_grid_cm[i+1] > x_shell_end_points[i_tmp]
            n_shell_end_points[i_tmp] = i
            i_tmp += 1
        end
    end

    for i in 1:(num_UpS_shells+num_DwS_shells)
        println(outfile, "x boundaries[rg0] for photon shell ", i, ":",
                x_shell_end_points[i]/rg0, " -> ", x_shell_end_points[i+1]/rg0, "  (",
                n_shell_end_points[i], "->", n_shell_end_points[i+1], ")")
    end

    println(outfile)
end

i_grid_feb = findfirst(>(feb_UpS), x_grid_cm) - 1

# Because redshift will be needed to compute radiative losses (it affects both the
# energy and density of CMB photons), calculate it here. If redshift was provided
# during data_input, cosmo_calc will return distance instead
# Cosmo_calc expects distance in megaparsecs, so convert from value read in during data_input
jet_dist_Mpc = jet_dist_kpc * 1e-3
# TODO use cosmology.jl for this instead?
#jet_dist_Mpc, redshift = cosmo_calc_wrapper(jet_dist_Mpc, redshift)
if cosmo_var == 1
    redshift = get_redshift(jet_dist_Mpc)
else
    jet_dist_Mpc = ustrip(u"Mpc", comoving_radial_distance(cosmo_calc_params.cosmo, redshift))
end


# Set a handful of constants related to radiative losses. electron_rm will not be the same
# as E₀_electron in module "constants" unless (a) there are electrons in the run, and
# (b) they are true electrons, with aa = xme/mp_cgs. (electron_rm may in fact = E₀_proton,
# but in that case it won't be used because radiative losses will never be calculated.)
# To convert rad_loss_fac to dp/dt, multiply by p^2*B^2, both in cgs. Prefactor of
# rad_loss_fac comes from average over pitch and is Eq (16) of Sturner+ (1997)
# [1997ApJ...490..619S]. Note the extra factor of c in the denominator, because code tracks
# dp/dt, not dE/dt as given in Sturner+ (1997).
electron_mass_ratio = minimum(aa_ion)
electron_rm = electron_mass_ratio * E₀_proton
σ_T = 8π/3 * (qp_cgs^2 / electron_rm)^2 # Thomson scattering cross-section
rad_loss_fac = 4/3 * c_cgs * σ_T /  (c_cgs^3 * (mp_cgs*electron_mass_ratio)^2 * 8π)
B_CMBz = B_CMB0 * (1 + redshift)^2


# Zero out total escaping fluxes and calculate the far UpS fluxes
#px_esc_flux_UpS_tot = 0.0
#energy_esc_flux_UpS_tot = 0.0
px_esc_flux_UpS     = zeros(na_itrs)
energy_esc_flux_UpS = zeros(na_itrs)
(flux_px_UpS, flux_pz_UpS, flux_energy_UpS) = upstream_fluxes(
    oblique, n_ions, denZ_ion, tZ_ion, aa_ion, zz_ion, sc_electron, tZ_electron,
    bmag_Z, θ_BZ, γ_Z, β_Z, u_Z
)

# Determine upstream Mach numbers (sonic & Alfven)
mach_sonic, mach_alfven = upstream_machs(
    β_Z, n_ions, denZ_ion, tZ_ion, aa_ion, zz_ion, sc_electron, tZ_electron, bmag_Z
)


# Set up the initial shock profile, or read it in from a file
if ! do_old_prof
    (ux_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
     β_ef_grid, γ_ef_grid, btot_grid, θ_grid, εB_grid, bmag_2) = setup_profile(
        u_Z, β_Z, γ_Z, bmag_Z, θ_BZ, r_comp, bturb_comp_frac, bfield_amp, use_custom_εB,
        n_ions, aa_ion, denZ_ion, sc_electron, zz_ion, flux_px_UpS, flux_energy_UpS, n_grid,
        x_grid_cm, x_grid_rg,
)
else
    error("Reading old profiles not yet supported")
    (
     x_grid_rg, x_grid_cm, ux_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid, β_ef_grid, γ_ef_grid,
     btot_grid, εB_grid, θ_grid, n_grid, u_Z, γ_Z, rg0, r_comp, rRH, β_Z, bmag_Z,
     u_2, β_2, γ_2, θ_u2, bmag_2, θ_BZ, θ_B2, flux_px_UpS, flux_pz_UpS, flux_energy_UpS
    ) = read_old_prof(n_old_skip, n_old_profs, n_old_per_prof)

    # Must set far UpS and DwS limits manually, since they won't be read in from the file
    x_grid_rg[0]        = -1e30
    x_grid_rg[n_grid+1] =  1e30
    x_grid_cm[0]        = -1e30 * rg0
    x_grid_cm[n_grid+1] =  1e30 * rg0
end


# Find the location of the shock
i_shock = findlast(≤(0), x_grid_rg)
isnothing(i_shock) && error("Shock location not found")


# How frequently will the code print during initial particle propagation?
n_print_pt = n_pts_inj < 500 ? 25 : 250


# The random number seed depends on the specific particle.
# Compute a necessary quantity to determine that seed
n_pts_max = max(n_pts_pcut, n_pts_pcut_hi)


# Because electrons might have different Monte Carlo weights than protons,
# protons, set that ratio here
# WARNING: assumes electrons are last ion species of input file
electron_wt_fac = sc_electron ? 1 / denZ_ion[n_ions] : 0.0


# Print a bunch of data about the run to screen/file
print_input(n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions, num_psd_mom_bins,
            num_psd_θ_bins, n_xspec, n_pcuts, n_grid, rRH, r_comp,
            u_Z, β_Z, γ_Z, u_2, β_2, γ_2, denZ_ion, bmag_Z, bmag_2, θ_BZ, θ_B2, θ_u2,
            aa_ion, tZ_ion, tZ_electron, mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine,
            feb_UpS, feb_DwS, rg0, age_max, energy_pcut_hi, do_fast_push, bturb_comp_frac,
            sc_electron, outfile)

wts_file = jldopen("mc_coupled_wts.hdf5", "a+")
spectra_file = jldopen("mc_coupled_spectra.hdf5", "a+")

pressure_psd_par   = Vector{Float64}(undef, na_grid)
pressure_psd_perp  = Vector{Float64}(undef, na_grid)
energy_density_psd = Vector{Float64}(undef, na_grid)

energy_transfer_pool = Vector{Float64}(undef, na_grid)
energy_recv_pool     = Vector{Float64}(undef, na_grid)

energy_density       = Matrix{Float64}(undef, na_c, na_ions)
therm_energy_density = Matrix{Float64}(undef, na_c, na_ions)

psd = OffsetArray{Float64}(undef, (0:psd_max, 0:psd_max, na_grid))

begin # "module" species_vars
    # Arrays will hold crossing data for thermal particles; px and pt are shock frame values
    #integer :: n_cr_count
    num_crossings = zeros(Int, na_grid)
    therm_grid = zeros(Int, na_cr)
    therm_px_sk = zeros(na_cr)
    therm_pt_sk = zeros(na_cr)
    therm_xwt = zeros(na_cr)

    # Spectra at x_spec locations
    spectra_sf = zeros(0:psd_max, na_grid)
    spectra_pf = zeros(0:psd_max, na_grid)

    # Escaping spectra UpS and DwS from shock; 2-D arrays store angular information
    esc_spectra_feb_UpS = zeros(0:psd_max)
    esc_spectra_feb_DwS = zeros(0:psd_max)
    esc_psd_feb_UpS = zeros(0:psd_max, 0:psd_max)
    esc_psd_feb_DwS = zeros(0:psd_max, 0:psd_max)
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
    pion_photon_sum = zeros(na_photons, na_grid)
    ic_photon_sum = zeros(na_photons, na_grid)
end # "module" photons

xwt_new = zeros(na_particles)
pt_pf_new = zeros(na_particles)
pb_pf_new = zeros(na_particles)
x_PT_cm_new = zeros(na_particles)

pcuts_use = zeros(na_c)

begin # "module" pcut_vars
    l_save = zeros(Bool, na_particles) # Whether or not to save particle for next pcut
    grid_sav        = zeros(Int, na_particles)
    tcut_sav        = zeros(Int, na_particles)
    DwS_sav         = zeros(Bool, na_particles)
    inj_sav         = zeros(Bool, na_particles)
    xwt_sav         = zeros(na_particles)
    pt_pf_sav        = zeros(na_particles)
    pb_pf_sav        = zeros(na_particles)
    x_PT_cm_sav     = zeros(na_particles)
    xn_per_sav      = zeros(na_particles)
    zz_sav          = zeros(na_particles)
    prp_x_cm_sav    = zeros(na_particles)
    acctime_sec_sav = zeros(na_particles)
    φ_rad_sav       = zeros(na_particles)

    grid_new        = zeros(Int, na_particles)
    tcut_new        = zeros(Int, na_particles)
    DwS_new         = zeros(Bool, na_particles)
    inj_new         = zeros(Bool, na_particles)
    xn_per_new      = zeros(na_particles)
    zz_new          = zeros(na_particles)
    prp_x_cm_new    = zeros(na_particles)
    acctime_sec_new = zeros(na_particles)
    φ_rad_new       = zeros(na_particles)
end

global r_PT_old = SVector(0.0, 0.0, 0.0)

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#  Main computational loops
#
#  Loops are set in the following sequence of nesting
#   1. loop_itr:   i_iter    Iteration number
#   2. loop_ion:   i_ion    Particle species number
#   3. loop_pcut:  i_cut    Particle splitting (pcut) number
#   4. loop_pt:    i_prt    Individual particle number
#
#
# Start of loop over iterations
#--------------------------------------------------------------------------
for i_iter in 1:n_itrs # loop_itr

    # Zero out numerous quantities that will be modified over the course of
    # this iteration. Minimally positive number is used to prevent errors
    # when taking logarithms later
    fill!(pxx_flux, 1e-99)
    fill!(pxz_flux, 1e-99)
    fill!(energy_flux, 1e-99)

    fill!(esc_spectra_feb_UpS, 1e-99)
    fill!(esc_spectra_feb_DwS, 1e-99)

    fill!(pressure_psd_par, 1e-99)
    fill!(pressure_psd_perp, 1e-99)
    fill!(energy_density_psd, 1e-99)

    fill!(wt_coupled, 1e-99)


    # Additionally, set/reset scalar quantities that will change
    sum_pressure_DwS  = 1e-99
    sum_KEdensity_DwS = 1e-99

    energy_esc_UpS    = 1e-99
    px_esc_UpS        = 1e-99


    # To facilitate energy transfer from ions to electrons, calculate here the target energy
    # density fraction for electrons at each grid zone, and zero out the pool of
    # plasma-frame energy that will be taken from ions and donated to electrons
    # Per Ardaneh+ (2015) [10.1088/0004-637X/811/1/57], ε_electron is proportional to √ε_b.
    # ε_b is itself roughly proportional to density² -- B² is proportional to z² (z being the
    # density compression factor) -- so ε_electron should vary roughly linearly with density.
    z_max = γ_Z * β_Z / (γ_2 * β_2)
    ε_target = zeros(na_grid)
    for i in 1:n_grid
        if ux_sk_grid[i] != u_Z
            z_curr = γ_Z * u_Z / (γ_sf_grid[i] * ux_sk_grid[i])
            ε_target[i] = energy_transfer_frac * (z_curr - 1) / (z_max - 1)
        end
    end

    fill!(energy_transfer_pool, 0.0)
    fill!(energy_recv_pool, 0.0)
    fill!(energy_density, 0.0)
    fill!(therm_energy_density, 0.0)

    #------------------------------------------------------------------------
    #  Start of loop over particle species
    #
    #  First species always protons. If electrons present they MUST be last species.
    #  Each species has mass number "aa" and charge number "zz" in units of
    #  proton mass and charge, respectively.
    #------------------------------------------------------------------------
    for i_ion in 1:n_ions # loop_ion

        aa = aa_ion[i_ion]
        zz = zz_ion[i_ion]

        o_o_mc = 1 / (aa*mp_cgs * c_cgs)

        # At the start of each ion, print a glyph to the screen
        @info("Starting species iteration", i_iter, i_ion)

        # Determine the pcut at which to switch from low-E particle counts to high-E
        # particle counts. Recall that energy_pcut_hi has units of keV per aa, so when
        # dividing by particle mass the factor of aa is already present in the denominator.
        # Also set the maximum momentum cutoff based on the values given in keyword "ENMAX"
        if (energy_pcut_hi * keV2erg / E₀_proton) < energy_rel_pt
            p_pcut_hi = √( 2 * energy_pcut_hi * keV2erg / E₀_proton )
        else
            p_pcut_hi = aa * mp_cgs * c_cgs * √( (energy_pcut_hi*keV2erg/E₀_proton + 1)^2  -  1 )
        end

        if Emax_keV > 0
            γ = 1 + (Emax_keV*keV2erg)/(aa*E₀_proton)
            pmax_cutoff = aa*mp_cgs * c_cgs * √( γ^2 - 1 )
        elseif Emax_keV_per_aa > 0
            γ = 1 + (Emax_keV_per_aa*keV2erg/E₀_proton)
            pmax_cutoff = aa*mp_cgs * c_cgs * √( γ^2 - 1 )
        elseif pmax_cgs > 0
            pmax_cutoff = pmax_cgs
        else
            # Something has gone very wrong.
            error("Max CR energy not set in data_input, so can't set pmax_cutoff.")
        end

        # Zero out the phase space distributions and set variables related to
        # tracking thermal particles
        n_cr_count = 0
        fill!(num_crossings,   0)
        fill!(therm_grid,      0)
        fill!(psd,             1e-99)
        fill!(esc_psd_feb_UpS, 1e-99)
        fill!(esc_psd_feb_DwS, 1e-99)
        fill!(therm_px_sk,     0.0)
        fill!(therm_pt_sk,     0.0)
        fill!(therm_xwt,       0.0)

        # In addition to initializing the phase space distribution, open the scratch (i.e.
        # temporary) file to which we will write information about thermal particle grid crossings
        nc_unit = open("mc_crossings.dat", "a")

        # To maintain identical results between OpenMP and serial versions,
        # set RNG seed based on current iteration/ion/pcut/particle number
        iseed_mod = (i_iter - 1)*n_ions  +  (i_ion - 1)
        Random.seed!(iseed_mod)


        # Initialize the particle populations that will be propagated through
        # the shock structure
        global n_pts_use, i_grid_in, xwt_in, pt_pf_in, pb_pf_in, x_PT_cm_in, globals... = init_pop(
            do_fast_push, inp_distr, i_ion, aa,
            tZ_ion, energy_inj, inj_wt, n_pts_inj, denZ_ion, x_grid_start, rg0, η_mfp,
            x_fast_stop_rg, β_Z, γ_Z, u_Z, n_ions, aa_ion, zz_ion, sc_electron, tZ_electron, oblique,
            n_grid, x_grid_rg, ux_sk_grid, γ_sf_grid,
            ptot_inj, wt_inj, n_pts_MB,
        )
        global pxx_flux = globals[1]
        global pxz_flux = globals[2]
        global energy_flux = globals[3]
        @info("Finished init_pop on",
               i_iter, i_ion,
               n_pts_use, i_grid_in, xwt_in, pt_pf_in,
               pb_pf_in, x_PT_cm_in, pxx_flux, pxz_flux, energy_flux)

        # Assign the various particle properties to the population
        xwt_new[1:n_pts_use]     .= xwt_in[1:n_pts_use]
        pt_pf_new[1:n_pts_use]   .= pt_pf_in[1:n_pts_use]
        pb_pf_new[1:n_pts_use]   .= pb_pf_in[1:n_pts_use]
        x_PT_cm_new[1:n_pts_use] .= x_PT_cm_in[1:n_pts_use]
        grid_new[1:n_pts_use]    .= i_grid_in[1:n_pts_use]

        DwS_new[1:n_pts_use]         .= false
        inj_new[1:n_pts_use]         .= false
        xn_per_new[1:n_pts_use]      .= xn_per_fine
        prp_x_cm_new[1:n_pts_use]    .= x_grid_stop
        acctime_sec_new[1:n_pts_use] .= 0.0
        tcut_new[1:n_pts_use]        .= 1

        φ_rad_new[1:n_pts_use] .= 2π*Random.rand(n_pts_use)

        # Weight of remaining particles, printed after each pcut; note that this will not be
        # correct for all particles if they were originally created so each thermal bin
        # would have equal weight
        global wt_running = xwt_in[1]


        # When using OpenMP, the array energy_transfer_pool can't be conditionally assigned
        # shared or reduction status, so it can't be used for both the ion and electron loops.
        # To get around this, use one array to hold the donated energy, and another to hold
        # the received energy.
        energy_recv_pool .= energy_transfer_pool


        # The array of pcuts read in by data_input has units momentum/mc.
        # Convert to momentum for this species
        pcuts_use[1:n_pcuts] .= pcuts_in[1:n_pcuts] * aa*mp_cgs*c_cgs

        #----------------------------------------------------------------------
        #  Start of loop over pcuts
        #
        #  Particle splitting and resetting of n_pts_use is handled by the call
        #  to "new_pcut" at the end of each iteration of loop_pcut.
        #----------------------------------------------------------------------
        for i_cut in 1:n_pcuts # loop_pcut

            # Initialize all of the *_sav arrays to help prevent bleedover between pcuts or ion species
            l_save .= false  # Whole array must be initialized in case number of particles changes from pcut to pcut
            xwt_sav         .= 0.0
            pt_pf_sav       .= 0.0
            pb_pf_sav       .= 0.0
            x_PT_cm_sav     .= 0.0
            grid_sav        .= 0
            DwS_sav         .= false
            inj_sav         .= false
            xn_per_sav      .= 0.0
            prp_x_cm_sav    .= 0.0
            acctime_sec_sav .= 0.0
            φ_rad_sav       .= 0.0
            tcut_sav        .= 0

            # A separate variable tracks the number of finished particles, so
            # that race conditions can be avoided in OMP mode
            i_fin = 0

            # For high-energy electrons in a strong magnetic field, need to know
            # previous cutoff momentum for calculating new PRP downstream
            pcut_prev = i_cut > 1 ? pcuts_use[i_cut-1] : 0.0


            #--------------------------------------------------------------------
            #  Start of loop over particles
            #
            #  Quick explanation of particle properties. All units cgs where a unit exists.
            #
            #  xwt: weighting value (used in momentum splitting)
            #  pt_pf: total plasma frame momentum
            #  pb_pf: momentum along B field in plasma frame. NOT along x-axis
            #         unless upstream orientation is parallel
            #  r_PT_cm: current particle position
            #  i_grid: current grid zone number
            #  l_DwS: whether particle has been downstream
            #  inj: whether particle has been back upstream (i.e. is a CR)
            #  xn_per: time steps per gyro period, Δt = T_g/xn_per
            #  prp_x_cm: DwS position of PRP; adjusted to allow all particles,
            #            regardless of momentum, to isotropize before reaching
            #  acctime_sec: time since crossing shock for first time; not started until l_DwS = true
            #  φ_rad: phase angle of gyration rel. to z axis
            #  tcut_curr: next time at which particle tracking takes place
            #--------------------------------------------------------------------
            #$omp parallel for default(none), schedule(dynamic,1), num_threads(6)
            for i_prt in 1:n_pts_use # loop_pt

                # To maintain identical results between OpenMP and serial versions,
                # set RNG seed based on current iteration/ion/pcut/particle number
                iseed_mod =  (  (i_iter - 1)*n_ions*n_pcuts*n_pts_max
                              + (i_ion - 1)       *n_pcuts*n_pts_max
                              + (i_cut - 1)               *n_pts_max
                              +  i_prt)
                Random.seed!(iseed_mod)

                # Reset the counter for number of times through the main loop
                helix_count = 0


                # Get the properties of the particle we're about to treat
                xwt         = xwt_new[i_prt]
                pt_pf       = pt_pf_new[i_prt]
                pb_pf       = pb_pf_new[i_prt]
                i_grid      = grid_new[i_prt]
                i_grid_old  = i_grid  # Needed for energy transfer
                l_DwS       = DwS_new[i_prt]
                inj         = inj_new[i_prt]
                xn_per      = xn_per_new[i_prt]
                prp_x_cm    = prp_x_cm_new[i_prt]
                acctime_sec = acctime_sec_new[i_prt]
                φ_rad       = φ_rad_new[i_prt]
                tcut_curr   = tcut_new[i_prt]

                r_PT_cm = SVector{3,Float64}(
                    x_PT_cm_new[i_prt], # x
                    # Not currently tracked, but could be added in at later date
                    0.0, # y
                    0.0) # z

                γ_pt_pf     = √( 1  +  (pt_pf * o_o_mc)^2 )

                # Constant that will be used repeatedly during loop
                # Square root corresponds to Blandford-McKee solution, where
                #   e ∝ 1/χ ∝ 1/r
                gyro_denom = 1 / (zz*qp_cgs * btot_grid[i_grid])
                if use_custom_εB && (r_PT_cm.x > x_grid_stop)
                    gyro_denom *=  √(r_PT_cm.x / x_grid_stop)
                end

                # Gyroradius assuming all motion is perpendicular to B field;
                # pitch-angle-correct gyroradius is gyro_rad_cm
                global gyro_rad_tot_cm = pt_pf * c_cgs * gyro_denom

                # Gyroperiod in seconds
                gyro_period_sec = 2π * γ_pt_pf * aa*mp_cgs * c_cgs * gyro_denom


                # Get the properties of the grid zone the particle's in
                ux_sk = ux_sk_grid[i_grid]
                uz_sk = uz_sk_grid[i_grid]
                utot  = utot_grid[i_grid]
                γ_usf = γ_sf_grid[i_grid]
                γ_uef = γ_ef_grid[i_grid]
                β_uef = β_ef_grid[i_grid]
                bmag  = btot_grid[i_grid]
                bθ    = θ_grid[i_grid]

                b_sinθ = sin(bθ)
                b_cosθ = cos(bθ)


                #------------------------------------------------------------------
                # Helix loop: The following loop is a large-scale restructuring of
                # the original code's loops 5702 and 2001.
                # The original code was slightly clearer in intent, but the goto
                # statements made it impossible to parallelize for CPU computing;
                # this code has sacrificed some of the clarity in exchange for
                # potential speed gains down the road.
                #
                #    ORIGINAL                       THIS CODE
                #  5702 continue                 first_iter = true
                #
                #   Code Block 1                 while keep_looping
                #                                    if first_iter || condition 1
                #  2001 continue                         first_iter = false
                #                                        Code Block 1
                #   Code Block 2                     else
                #                                        Code Block 3
                #  if (condition 1) goto 5702        end
                #
                #   Code Block 3                     Code Block 2
                #                                end
                #  goto 2001
                #
                # Code Block 1: start of helix loop; minor setup
                # Code Block 2: particle movement; check on DwS PRP return
                # Code Block 3: grid zone changes; non-DwS escape (exception: if no
                #   scattering); radiative losses; momentum splitting; scattering
                #
                # loop_helix in this code is the old loop 2001. The first part of
                # loop 5702 has been subsumed into the if statement, which
                # determines whether Code Block 3 would be run or whether the
                # code would cycle back to run Code Block 1 again.
                # Following every combination of logical possibilities should
                # result in the same path through the two loops.
                #------------------------------------------------------------------
                keep_looping = true
                first_iter   = true
                # Control variable for "condition 1":
                #  -1: default
                #   0: particle escapes
                #   1: particle returns from DwS via convection
                #   2: particle didn't enter return calcs
                i_return = -1
                i_reason = 0
                lose_pt  = false

                while keep_looping # loop_helix
                    helix_count, keep_looping, first_iter, i_grid_old = helix_loop!(
                        helix_count, keep_looping, first_iter, n_cr_count, i_return, i_prt, i_grid, i_grid_old,
                        pt_pf, pb_pf, γ_pt_pf,
                        γ_usf,
                        gyro_denom, gyro_period_sec,
                        r_PT_cm, φ_rad, xn_per,
                        ux_sk, uz_sk, utot,
                        b_cosθ, b_sinθ,
                        aa, zz, o_o_mc,
                        prp_x_cm, acctime_sec, pcut_prev, tcut_curr, pmax_cutoff,
                        l_DwS, inj,
                        xwt, energy_esc_UpS, px_esc_UpS,
                        nc_unit,
                    )
                    @debug("helix_loop! returned with", helix_count, keep_looping, first_iter, i_grid_old)
                end # loop_helix
                #------------------------------------------------------------------
                # End of loop moving/tracking particles on/off grid


                if ! l_save[i_prt]
                    ENV["JULIA_DEBUG"] = "all"
                    particle_finish(aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad, ux_sk, uz_sk,
                                    utot, γ_usf, b_cosθ, b_sinθ, i_reason, xwt,
                                    esc_psd_feb_DwS, esc_psd_feb_UpS, esc_flux,
                                    px_esc_feb, energy_esc_feb, esc_energy_eff, esc_num_eff,
                                    i_iter, i_ion, oblique, o_o_mc)
                    ENV["JULIA_DEBUG"] = ""
                end


                # Particle counting
                if i_cut == 1 && ( i_prt == 1 || i_prt % n_print_pt == 0 )
                    @info("Particle = $i_prt")
                end

                i_fin += 1
                #if i_fin % 16 == 0
                #    print_progress_bar(i_fin, n_pts_use)
                #end

            end # loop_pt
            #$omp end parallel do
            #--------------------------------------------------------------------
            # Conclusion of particle loop

            n_saved = count(l_save)

            t_end = now()
            run_time = t_end - t_start

            @info("", i_iter, i_ion, i_cut,
                  pcuts_in[i_cut], pcuts_use[i_cut]/(mp_cgs*c_cgs),
                  n_saved, n_pts_use, wt_running, run_time)
            println(outfile, " itr=", i_iter, " ion=", i_ion, " icut=", i_cut,
                    pcuts_in[i_cut], pcuts_use[i_cut]/(mp_cgs*c_cgs),
                    "  n_sav=", n_saved, "/", n_pts_use, wt_running, run_time)


            # If no particles saved, don't bother with remaining pcuts
            n_saved == 0 && break # loop_pcut

            # Prepare population for next pcut
            n_pts_target = pcuts_use[i_cut] < p_pcut_hi ? n_pts_pcut : n_pts_pcut_hi

            global (
                grid_new, tcut_new, DwS_new, inj_new, xwt_new, pt_pf_new, pb_pf_new,
                x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new,
                n_pts_use, wt_running) = new_pcut(
                    n_pts_target, n_saved, l_save, grid_sav, DwS_sav, inj_sav,
                    xwt_sav, pt_pf_sav, pb_pf_sav, x_PT_cm_sav, xn_per_sav,
                    prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav,
                    n_pts_use, wt_running)


        end # loop_pcut
        #----------------------------------------------------------------------
        # Conclusion of pcuts loop

        println()
        println(outfile)

        # Obtain the dN/dp arrays for the current species: one 2-D array for each grid zone.
        get_normalized_dNdp(dNdp_therm, dNdp_therm_pvals, dNdp_cr, nc_unit, zone_pop)

        # Get pressure (both components) and kinetic energy density everywhere in the shock profile
        thermo_calcs(num_crossings, n_cr_count, therm_grid, therm_px_sk, therm_pt_sk,
                     therm_xwt, nc_unit, psd, zone_pop,
                     pressure_psd_par, pressure_psd_perp, energy_density_psd)

        # Transform just the electron PSD into the explosion frame, since we
        # will need it shortly to calculate inverse Compton emission
        get_dNdp_2D(nc_unit, zone_pop, d2N_dpdcos_ef)


        # Output the spectra associated with this ion species
        #TODO: spectrum_plot

        # Print out escaping particle population for this species
        print_dNdp_esc(esc_psd_feb_UpS, esc_psd_feb_DwS)

        # Handle photon calculations
        if do_photons
            photon_calcs(dNdp_therm_pvals, dNdp_therm, dNdp_cr, aa, n_shell_end_points, d2N_dpdcos_ef)
        end

    end # loop_ion
    #DEBUGLINE
    for i in 1:n_grid
        @debug(i_iter, i,
               #x_grid_rg[i],
               therm_energy_density[i,1:n_ions]/zone_vol[i],
               therm_energy_density[i,n_ions]/max(sum(therm_energy_density[i,1:n_ions]),1e-99),
               sum(therm_energy_density[i,1:n_ions]/zone_vol[i]),
               energy_density[i,1:n_ions]/zone_vol[i],
               energy_density[i,n_ions]/max(sum(energy_density[i,1:n_ions]),1e-99),
               sum(energy_density[i,1:n_ions]/zone_vol[i]))
    end
    #print_plot_vals(888)
    #------------------------------------------------------------------------
    # Conclusion of species loop

    # Compute the escaping flux for this iteration
    px_esc_flux_UpS[i_iter] = px_esc_UpS / flux_px_UpS
    energy_esc_flux_UpS[i_iter] = energy_esc_UpS / flux_energy_UpS

    # Compute the adiabatic index everywhere on grid now that all particles
    # are accounted for. First store value from previous iteration, since
    # both the pre- and post-iteration adiabatic indices will be used in
    # the profile smoothing subroutine smooth_grid
    if i_iter == 1
        index_mask = (x_grid_cm .≤ 0)
        γ_adiab_grid[ index_mask,1] .= 5//3         #    #assumecold
        γ_adiab_grid[~index_mask,1] .= γ_adiab_2_RH
    else
        γ_adiab_grid[:,1] .= γ_adiab_grid[:,2]
    end

    γ_adiab_grid[:,2] .= @. 1 + (pressure_psd_par + pressure_psd_perp) / energy_density_psd

    index_mask = (energy_density_psd .== 1e-99)
    γ_adiab_grid[index_mask,2] .= 1e-99

    # Also compute the adiabatic index of particles that were lost DwS
    γ_adiab_DwS[i_iter] = 1 + sum_pressure_DwS/sum_KEdensity_DwS

    # Calculate expected escaping fluxes, now that adiabatic index is known
    # far DwS. Also average them so that the smoothing subroutine treats
    # calculated and actual escaping fluxes identically
    q_esc_cal_px[i_iter], q_esc_cal_energy[i_iter] = q_esc_calcs(γ_adiab_DwS[i_iter],)
    n_avg = min(i_iter, 4)
    q_esc_cal_px_avg = mean(q_esc_cal_px[i_iter-n_avg+1:i_iter])
    q_esc_cal_energy_avg = mean(q_esc_cal_energy[i_iter-n_avg+1:i_iter])

    # In runs with OpenMP, race conditions may cause roundoff error at the
    # least significant decimal place of the **_flux variables. This causes
    # slight variations in the smoothing algorithm across different runs.
    # Correct for that here by rounding off the last two decimal places
    # (the second place is a fudge factor).
    pxx_flux[1:n_grid] .= round(pxx_flux[1:n_grid], digits=13)
    # Commented out because it's not necessary for smoothing parallel shocks
    #pxz_flux[1:n_grid] = round(pxz_flux[1:n_grid], digits=13)
    energy_flux[1:n_grid] = round(energy_flux[1:n_grid], digits=13)

    # Output grid data for this iteration, and smooth the grid for the next iteration
    if oblique
        @error("ERROR: grid smoothing not coded for oblique shocks.")
        @error("Stopping program now.")
        exit()
    else
        smooth_grid_par(i_iter, i_shock, n_grid, x_grid_rg, x_grid_cm,
                        γ_adiab_grid, uz_sk_grid, θ_grid,
                        pressure_psd_par, pressure_psd_perp, flux_px_UpS,
                        flux_energy_UpS, γ_adiab_2_RH, q_esc_cal_px_avg,
                        q_esc_cal_energy_avg, pxx_flux, energy_flux, ux_sk_grid,
                        γ_sf_grid, btot_grid, utot_grid, γ_ef_grid,
                        β_ef_grid, εB_grid)
    end


    # Compute average escaping flux over last four iterations and write to
    # file; average over all iterations if i_iter < 4.
    n_avg          = min(i_iter, 4)
    px_esc_avg     = mean(px_esc_flux_UpS[1:n_avg])
    energy_esc_avg = mean(energy_esc_flux_UpS[1:n_avg])

    println(outfile, " Parallel shock q_esc from Double et al (2004) equations:")
    println(outfile, "     Esc. energy flux/UpS    = ", q_esc_cal_energy_avg)
    println(outfile, "     Esc. momentum flux/UpS  = ", q_esc_cal_px_avg)
    energy_esc_flux_UpS[i_iter] = max(energy_esc_flux_UpS[i_iter], 1e-99)
    energy_esc_avg             = max(energy_esc_avg,             1e-99)
    px_esc_flux_UpS[i_iter] = max(px_esc_flux_UpS[i_iter], 1e-99)
    px_esc_avg             = max(px_esc_avg,             1e-99)
    println(outfile,
            " Esc. en flux FEB/UpS  for i_iter = ", i_iter, ":   en esc = ",
            energy_esc_flux_UpS[i_iter], "   Avg. esc en  = ", energy_esc_avg)
    println(outfile,
            " Esc. pxx flux FEB/UpS for i_iter = ", i_iter, ":  pxx esc = ",
            px_esc_flux_UpS[i_iter], "   Avg. esc pxx = ", px_esc_avg)

    if iszero(q_esc_cal_px_avg)
        println(outfile, " Avg q_px_MC/q_px_cal N/A, because q_px_cal = 0")
    else
        println(outfile, " Avg q_px_MC/q_px_cal. = ", px_esc_avg/q_esc_cal_px_avg)
    end
    if iszero(q_esc_cal_energy_avg)
        println(outfile, " Avg q_energy_MC/q_energy_cal N/A, because q_energy_cal = 0")
    else
        println(outfile, " Avg q_energy_MC/q_energy_cal. = ", energy_esc_avg/q_esc_cal_energy_avg)
    end
    println(outfile)


    # Compute various adiabatic indices and write them out
    n_avg     = min(i_iter, 4)
    γ_DwS_esc = mean(γ_adiab_DwS[i_iter-n_avg+1:i_iter])
    γ_UpS     = 5//3 # #assumecold

    println(outfile, " Iteration #", i_iter)
    println(outfile, "   r_comp = ", r_comp, "      r_RH = ",rRH)
    println(outfile, "   Adiab index for far UpS particles = ", γ_UpS)
    println(outfile, "   Adiab index for DwS PRP particles = ", γ_DwS_esc)
    println(outfile, "   Adiab index from R-H relations    = ", γ_adiab_2_RH)
    println(outfile)


    # If tcut tracking was enabled, print out the results here
    if do_tcuts
        tcut_print(
            wts_file, spectra_file, n_tcuts, tcuts, n_ions,
            num_psd_mom_bins, psd_mom_bounds,
            i_iter, wt_coupled, spectra_coupled,
        )
    end

end # loop_itr
#--------------------------------------------------------------------------
# Conclusion of iteration loop

close(wts_file)
close(spectra_file)

# End the run
t_end = now()

run_time = round(t_end - t_start, Second)

println()
@info(" Finished. Run time = $(run_time) sec, $(round(run_time, Minute)) min")

println(outfile)
println(outfile, " Finished. Run time = ", run_time, ", ", round(run_time, Minute))
println(outfile)

# vim: set textwidth=92:
