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
begin
    include("iter_init.jl")
    include("iter_finalize.jl")
    include("particle_loop.jl")
    include("pcut_finalize.jl")
    include("ion_init.jl")
    include("ion_finalize.jl")
end

# Start the wall clock for this run
t_start = now()

outfile = open("mc_out.dat", "w")

# Get input, control variables, etc.
include("data_input.jl") # FIXME

# Set up the computational grid
(
 x_grid_rg, x_grid_start, x_grid_stop
) = setup_grid(x_grid_start_rg, x_grid_stop_rg, use_prp, feb_DwS, rg₀)
const n_grid = length(x_grid_rg) - 2
const grid_axis = axes(x_grid_rg, 1)
x_grid_cm = x_grid_rg * rg₀ # Convert everything from rg₀ units to cgs units
Γ_grid = zeros(n_grid, 2)

# Set quantities related to the phase space distribution, including the bins
psd_cos_fine = 1 - 2 / (psd_lin_cos_bins+1)
psd_θ_fine = acos(psd_cos_fine)
psd_θ_min  = psd_θ_fine / 10^psd_log_θ_decs
@debug("Phase space angle parameters", psd_cos_fine, psd_θ_fine, psd_θ_min)

if inp_distr == 1
    # Set minimum PSD energy using thermal distribution for UpS plasma
    Emin_keV = ustrip(u"keV", kB_cgs * minimum(T₀_ion[1:n_ions])*u"erg")

    # Allow for a few extra zones below the thermal peak
    Emin_keV *= emin_therm_fac

elseif inp_distr == 2
    # Set minimum PSD energy using δ-function dist for UpS plasma;
    # allow for a few extra zones below the location of the distribution
    Emin_keV = 0.2 * energy_inj
end

# Determine minimum momentum associated with the given energy, which will occur for
# the lightest particle species. Use a cutoff of 0.1% of the rest-mass energy for
# relativistic/nonrelativistic calculation.
aa_min = minimum(aa_ion)
if ustrip(u"erg", Emin_keV*u"keV") < 1e-3*aa_min*E₀_proton
    psd_mom_min = √(2 * aa_min*mₚ_cgs * ustrip(u"erg", Emin_keV*u"keV"))
else
    γ           = 1 + ustrip(u"erg", Emin_keV*u"keV")/(aa_min*E₀_proton)
    psd_mom_min = aa_min * mₚ_cgs * c_cgs * √(γ^2 - 1)
end

# Now find the maximum momentum for the PSD (this will be adjusted due to SF->PF Lorentz
# transformation). How to actually calculate it depends on the user-specified maximum energy
# "ENMAX"
aa_max = maximum(aa_ion)
if Emax_keV > 0
    γ          = 1 + ustrip(u"erg", Emax_keV*u"keV")/(aa_max*E₀_proton)
    psd_mom_max = aa_max * mₚ_cgs * c_cgs * √(γ^2 - 1)
elseif Emax_keV_per_aa > 0
    γ          = 1 + ustrip(u"erg", Emax_keV_per_aa*u"keV")/E₀_proton
    psd_mom_max = aa_max * mₚ_cgs * c_cgs * √(γ^2 - 1)
elseif pmax_cgs > 0
    psd_mom_max = pmax_cgs
else
    # Something has gone very wrong.
    error("Max CR energy not set in data_input, so can not set PSD bins.")
end

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
    x_shell_midpoints, x_shell_endpoints = set_photon_shells(num_UpS_shells, num_DwS_shells, use_prp,
                                                             feb_UpS, feb_DwS, rg₀, x_grid_stop_rg)
end

begin # "module" iteration_vars
    pxx_flux = Vector{Float64}(undef, n_grid)
    pxz_flux = Vector{Float64}(undef, n_grid)
    energy_flux  = Vector{Float64}(undef, n_grid)
    esc_flux = zeros(na_ions)
    px_esc_feb = zeros(na_ions, na_itrs)
    energy_esc_feb = zeros(na_ions, na_itrs)
    esc_energy_eff = zeros(0:psd_max, na_ions)
    esc_num_eff = zeros(0:psd_max, na_ions)

    # Arrays for holding thermal distribution information; they're set at the start of
    # the run, but included here because of chance they could change due to fast push
    n_pts_MB   = zeros(Int, na_ions)
    ptot_inj   = zeros(na_particles, na_ions)
    weight_inj = zeros(na_particles, na_ions)
    # Arrays for holding information about particle counts and spectra at various tcuts
    weight_coupled  = Matrix{Float64}(undef, na_c, na_ions)
    spectra_coupled = zeros(0:psd_max, na_c, na_ions)
end # "module" iteration_vars



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
# the location of the UpS FEB. Add the photon shells to the output file
n_shell_endpoints = zeros(Int, num_UpS_shells+num_DwS_shells+1)
if do_photons
    let i_tmp = 1
        for i in 1:n_grid
            if x_grid_cm[i] ≤ x_shell_endpoints[i_tmp] && x_grid_cm[i+1] > x_shell_endpoints[i_tmp]
                n_shell_endpoints[i_tmp] = i
                i_tmp += 1
            end
        end
    end

    for i in 1:(num_UpS_shells+num_DwS_shells)
        println(outfile, "x boundaries[rg₀] for photon shell ", i, ":",
                x_shell_endpoints[i]/rg₀, " -> ", x_shell_endpoints[i+1]/rg₀, "  (",
                n_shell_endpoints[i], "->", n_shell_endpoints[i+1], ")")
    end

    println(outfile)
end

i_grid_feb = findfirst(>(feb_UpS), x_grid_cm) - 1

# Because redshift will be needed to compute radiative losses (it affects both the energy
# and density of CMB photons), calculate it here. If redshift was provided during
# data_input, cosmo_calc will return distance instead
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
# (b) they are true electrons, with aa = mₑ/mₚ. (electron_rm may in fact = E₀_proton,
# but in that case it won't be used because radiative losses will never be calculated.)
# To convert rad_loss_fac to dp/dt, multiply by p²B², both in cgs. Prefactor of rad_loss_fac
# comes from average over pitch and is Eq (16) of Sturner+ (1997) [1997ApJ...490..619S].
# Note the extra factor of c in the denominator, because code tracks dp/dt, not dE/dt as
# given in Sturner+ (1997).
σ_T = 8π/3 * (qₚ_cgs^2 / E₀_electron)^2 # Thomson scattering cross-section
rad_loss_fac = 4//3 * c_cgs * σ_T / (c_cgs^3 * mₑ_cgs^2 * 8π)
B_CMBz = B_CMB0 * (1 + redshift)^2


# Zero out total escaping fluxes and calculate the far UpS fluxes
#px_esc_flux_UpS_tot = 0.0
#energy_esc_flux_UpS_tot = 0.0
px_esc_flux_UpS     = zeros(na_itrs)
energy_esc_flux_UpS = zeros(na_itrs)
(
 flux_px_UpS, flux_pz_UpS, flux_energy_UpS
) = upstream_fluxes(oblique, n_ions, ρ_N₀_ion, T₀_ion, aa_ion,
                    bmag₀, θ_B₀, u₀, β₀, γ₀)

# Determine upstream Mach numbers (sonic & Alfvén)
mach_sonic, mach_alfven = upstream_machs(β₀, n_ions, ρ_N₀_ion, T₀_ion, aa_ion, bmag₀)


# Set up the initial shock profile, or read it in from a file
if ! do_old_prof
    (
     ux_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid,
     β_ef_grid, γ_ef_grid, btot_grid, θ_grid, εB_grid, bmag₂,
    ) = setup_profile(
                      u₀, β₀, γ₀, bmag₀, θ_B₀, r_comp, bturb_comp_frac, bfield_amp, use_custom_εB,
                      n_ions, aa_ion, ρ_N₀_ion, flux_px_UpS, flux_energy_UpS, grid_axis,
                      x_grid_cm, x_grid_rg,
                     )
else
    error("Reading old profiles not yet supported")
    (
     x_grid_rg, x_grid_cm, ux_sk_grid, uz_sk_grid, utot_grid, γ_sf_grid, β_ef_grid, γ_ef_grid,
     btot_grid, εB_grid, θ_grid, n_grid, u₀, γ₀, rg₀, r_comp, r_RH, β₀, bmag₀,
     u₂, β₂, γ₂, θᵤ₂, bmag₂, θ_B₀, θ_B₂, flux_px_UpS, flux_pz_UpS, flux_energy_UpS
    ) = read_old_prof(n_old_skip, n_old_profs, n_old_per_prof)

    # Must set far UpS and DwS limits manually, since they won't be read in from the file
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
electron_weight_fac = 1 / ρ_N₀_ion[end]


# Print a bunch of data about the run to screen/file
print_input(n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,
            NaN, NaN, n_xspec, n_pcuts, n_grid, r_RH, r_comp,
            u₀, β₀, γ₀, u₂, β₂, γ₂, ρ_N₀_ion, bmag₀, bmag₂, θ_B₀, θ_B₂, θᵤ₂,
            aa_ion, T₀_ion, mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine,
            feb_UpS, feb_DwS, rg₀, age_max, energy_pcut_hi, do_fast_push, bturb_comp_frac,
            outfile)

weights_file = jldopen("mc_coupled_weights.hdf5", "a+")
spectra_file = jldopen("mc_coupled_spectra.hdf5", "a+")

pressure_psd_par   = Vector{Float64}(undef, n_grid)
pressure_psd_perp  = Vector{Float64}(undef, n_grid)
energy_density_psd = Vector{Float64}(undef, n_grid)

energy_transfer_pool = Vector{Float64}(undef, n_grid)
energy_recv_pool     = Vector{Float64}(undef, n_grid)

energy_density       = Matrix{Float64}(undef, na_c, na_ions)
therm_energy_density = Matrix{Float64}(undef, na_c, na_ions)

psd = OffsetArray{Float64}(undef, (psd_mom_axis, psd_θ_axis, n_grid))

begin # "module" species_vars
    # Arrays will hold crossing data for thermal particles; px and pt are shock frame values
    #integer :: n_cr_count
    num_crossings = zeros(Int, n_grid)
    therm_grid = zeros(Int, na_cr)
    therm_px_sk = zeros(na_cr)
    therm_pt_sk = zeros(na_cr)
    therm_weight = zeros(na_cr)

    # Spectra at x_spec locations
    spectra_sf = zeros(0:psd_max, n_grid)
    spectra_pf = zeros(0:psd_max, n_grid)

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
    pion_photon_sum = zeros(na_photons, n_grid)
    ic_photon_sum = zeros(na_photons, n_grid)
end # "module" photons

weight_new = zeros(na_particles)
ptot_pf_new = zeros(na_particles)
pb_pf_new = zeros(na_particles)
x_PT_cm_new = zeros(na_particles)

pcuts_use = zeros(na_c)

begin # "module" pcut_vars
    l_save = zeros(Bool, na_particles) # Whether or not to save particle for next pcut
    grid_sav        = zeros(Int, na_particles)
    tcut_sav        = zeros(Int, na_particles)
    DwS_sav         = zeros(Bool, na_particles)
    inj_sav         = zeros(Bool, na_particles)
    weight_sav      = zeros(na_particles)
    ptot_pf_sav     = zeros(na_particles)
    pb_pf_sav       = zeros(na_particles)
    x_PT_cm_sav     = zeros(na_particles)
    xn_per_sav      = zeros(na_particles)
    #zz_sav         = zeros(na_particles)
    prp_x_cm_sav    = zeros(na_particles)
    acctime_sec_sav = zeros(na_particles)
    φ_rad_sav       = zeros(na_particles)

    grid_new        = zeros(Int, na_particles)
    tcut_new        = zeros(Int, na_particles)
    DwS_new         = zeros(Bool, na_particles)
    inj_new         = zeros(Bool, na_particles)
    xn_per_new      = zeros(na_particles)
    #zz_new         = zeros(na_particles)
    prp_x_cm_new    = zeros(na_particles)
    acctime_sec_new = zeros(na_particles)
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
    iter_init()

    #------------------------------------------------------------------------
    #  Start of loop over particle species
    #
    #  First species always protons. If electrons present they MUST be last species.
    #  Each species has mass number "aa" and charge number "zz" in units of
    #  proton mass and charge, respectively.
    #------------------------------------------------------------------------
    for i_ion in 1:n_ions # loop_ion

        ion_init()

        #----------------------------------------------------------------------
        #  Start of loop over pcuts
        #
        #  Particle splitting and resetting of n_pts_use is handled by the call
        #  to "new_pcut" at the end of each iteration of loop_pcut.
        #----------------------------------------------------------------------
        for i_cut in 1:n_pcuts # loop_pcut

            pcut_init()

            #--------------------------------------------------------------------
            #  Start of loop over particles
            #
            #  Quick explanation of particle properties. All units cgs where a unit exists.
            #
            #  weight: weighting value (used in momentum splitting)
            #  ptot_pf: total plasma frame momentum
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
            #  φ_rad: phase angle of gyration relative to z axis
            #  tcut_curr: next time at which particle tracking takes place
            #--------------------------------------------------------------------
            #$omp parallel for default(none), schedule(dynamic,1), num_threads(6)
            for i_prt in 1:n_pts_use # loop_pt

                particle_loop()

            end # loop_pt
            #$omp end parallel do
            #--------------------------------------------------------------------
            # Conclusion of particle loop

            break_pcut = pcut_finalize()
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

# vim: set textwidth=92:shiftwidth=2:
