using LinearAlgebra: dot
using Unitful: ustrip, MeV, erg, c, ħ
using UnitfulGaussian: qcgs
using QuadGK: quadgk
using SpecialFunctions: besselk
using .constants: E₀ₚ, E₀ₑ

"""
    synch_emission(...)

This subroutine takes an electron distribution and calculates the synchrotron emission.

### Arguments

- `i_grid`: current grid zone, used to find magnetic field strength
- `num_hist_bins`: number of momentum bins in the distribution of thermal particles
- `p_pf_therm`: momentum boundary values, cgs units, of thermal distribution histogram
- `dN_therm`: thermal particle distribution. It is a pure number of particles in the spectral region
- `num_psd_mom_bins`: number of momentum bins in the distribution of accelerated particles
- `p_pf_cr`: momentum boundary values, cgs units, of cosmic ray distribution histogram
- `dN_cr`: cosmic ray distribution. It is a pure number of particles in the spectral region
- `n_photon_synch`: number of energy bins to use for photon production
- `photon_synch_min_MeV`: minimum photon energy, in MeV, to use for synchrotron spectrum
- `bins_per_dec_photon`: number of energy bins per decade of photon spectrum

### Returns

- `energy_γ`: energy bin values of resultant photon distribution (in ergs)
- `synch_emis`: emitted synchrotron spectrum, units of erg/s
"""
function synch_emission(
        i_grid, num_hist_bins, p_pf_therm, dN_therm,
        num_psd_mom_bins, p_pf_cr, dN_cr, n_photon_synch,
        photon_synch_min_MeV, bins_per_dec_photon,
        n_ions, aa_ion, n₀_ion, γ₀, u₀, flux_px_upstream, flux_energy_upstream, u₂,
        n_grid, btot_grid,
        i_ion, mc,
    )

    # Minimum energy for synch photon spectrum.
    # Max set by value of n_γ_synch and bins_per_dec_photon.
    energy_γ_min_log = log10(ustrip(erg, photon_synch_min_MeV*MeV)) # Min photon energy in log(ergs)
    Δγ = 1 / bins_per_dec_photon

    # Get strength of magnetic field, which will be used to find p_fac below. Note that
    # i_grid may be greater than n_grid. This only happens when this subroutine is called
    # to find the synchrotron emission for SSC purposes, which requires special treatment:
    # assume εB = .001 and calculate necessary magnetic field.
    if i_grid ≤ n_grid
        bmag_curr = btot_grid[i_grid]
    else
        n₀ = dot(n₀_ion, aa_ion)
        energy_density = (flux_energy_upstream + γ₀*u₀*n₀*E₀ₚ) / u₂ - flux_px_upstream
        bmag_curr = √(8π * 1e-3 * energy_density)
    end

    # Rybicki & Lightman eq. 6.18 without F factor.
    # units are power per unit frequency per electron (cgs)
    # Above: Note there is no sin(α) factor
    p_fac = √3/2π * (qcgs^3 * bmag_curr/E₀ₑ)


    # Initialize emission array and set energy of output photons
    synch_emis = fill(1e-99, n_photon_synch)
    energy_γ = exp10.(range(start = energy_γ_min_log, step = Δγ, length = n_photon_synch)) # Photon energy in ergs

    nu = 5//3

    # Now loop over momenta (i.e. energy) of electrons and calculate the synchrotron emission.
    #-------------------------------------------------------------------------
    # First up, the thermal particles.
    synch_emission_thermal_particles!(synch_emis, num_hist_bins, dN_therm, p_pf_therm, bmag_curr, mc, n_photon_synch, energy_γ, nu, p_fac)
    # With the thermal particles finished, move on to the cosmic ray population
    synch_emission_cosmic_ray!(synch_emis, num_psd_mom_bins, dN_cr, p_pf_cr, mc, bmag_curr, n_photon_synch, energy_γ, nu, p_fac)
    #-------------------------------------------------------------------------
    # Finished


    # Additional lines necessary for synchrotron self-Compton calculations in
    # subroutine ssc_cooling_array. Open a scratch file to hold data (if it
    # isn't already open), and write d²N/dEdt to it.
    #-------------------------------------------------------------------------
    synch_unit = 62
    lopen = inquire(:isopen, unit=synch_unit)
    if !lopen
        synch_unit = open(status="scratch", form="unformatted")
    end


    # Only write to the scratch file if there actually was synchrotron emission
    # produced by this grid zone
    if maximum(synch_emis) > 1e-90

        # Write the current grid zone number and number of synchrotron bins to the scratch file
        write(synch_unit, i_grid, n_photon_synch)

        # Convert synch_emis from dP/d(lnE) to d²N/dEdt by dividing (twice) by
        # photon energy, then write to the scratch file
        for i in 1:n_photon_synch
            if synch_emis[i] > 1e-90
                p_fac = max(synch_emis[i] / energy_γ[i]^2, 1e-99)
                write(synch_unit, energy_γ[i], p_fac)
            else
                write(synch_unit, energy_γ[i], 1e-99)
            end
        end

    end
    #-------------------------------------------------------------------------
    # Scratch file created/updated as needed

    return energy_γ, synch_emis
end


function synch_emission_thermal_particles!(
    synch_emis, num_hist_bins, dN_therm,
    p_pf_therm, bmag_curr, mc, n_photon_synch, energy_γ, nu, p_fac
)
    for i in 0:num_hist_bins-1

        # Total number of electrons in Δp
        xnum_electron = dN_therm[i]
        xnum_electron ≤ 1e-60 && continue # skip empty bins

        # Assume electrons with E < 3 MeV contribute no synchrotron emission
        p1 = √(p_pf_therm[i] * p_pf_therm[i+1]) # Geometric mean
        p1*c < 3MeV && continue

        # Lorentz factor for electron
        γ_electron = hypot(p1/mc, 1)

        # Eq. 6.17c Rybicki & Lightman without sin(α)
        ω_c = 3*(γ_electron^2)*qcgs*bmag_curr / 2mc

        # Calculate F factor in eq. 6.18 Rybicki & Lightman (see Eq 6.31c)
        #     F(x) ≡ x ∫_x^∞ K_{5/3}(ξ) dξ              (6.31c)
        for j in 1:n_photon_synch
            if bmag_curr < 1e-20 || ω_c < 1e-55
                F = 0.0
            else
                ω_γ = energy_γ[j] / ħ    # from E = ħω
                x = ω_γ/ω_c

                xxx_max_set = 30.0
                if x ≥ xxx_max_set || x < 1e-15
                    F = 0.0
                else
                    F = x * quadgk(t -> besselk(nu, t), x, Inf)
                end
            end

            # In the MC code, the electron spectra passed to synch_emission contain the
            # total number of particles in each momentum shell. Therefore, the emission
            # returned by this subroutine has units of erg/s. (It would be erg/(s⋅cm³)
            # if we had passed density instead of number.)
            # Also, Eq. (6.18) from Rybicki & Lightman, [xnum_electron * p_fac * F], is
            # energy production rate per frequency, dP/dω. Since ω_γ is E/ħ, dω = dE/ħ, and
            # so ω_γ/dω = E/dE. Then [dP/dω * ω_γ] is equal to [dP/dE * E], or dP/d(lnE).
            # This is what `synch_emis()` expects as output. Has units [erg/s].
            tmp_add = xnum_electron * ω_γ * p_fac * F

            # Only include emission if it's sufficiently positive
            if tmp_add > 1e-55
                synch_emis[j] += tmp_add
            end

        end # loop over n_photon_synch
    end
end

function synch_emission_cosmic_ray!(
    synch_emis, num_psd_mom_bins, dN_cr, p_pf_cr, mc, bmag_curr, n_photon_synch, energy_γ, nu, p_fac
)
    for i in 0:num_psd_mom_bins

        # Total number of electrons in Δp
        xnum_electron = dN_cr[i]
        xnum_electron ≤ 1e-60 && continue # skip empty bins

        # Assume electrons with E < 3 MeV contribute no synchrotron emission
        p1 = √(p_pf_cr[i] * p_pf_cr[i+1]) # Geometric mean
        p1*c < 3MeV && continue

        # Lorentz factor for electron
        γ_electron = hypot(p1/mc, 1)

        # Eq. 6.17c Rybicki & Lightman without sin(α)
        ω_c = 3*(γ_electron^2)*qcgs*bmag_curr / 2mc

        # Calculate F factor in eq. 6.18 Rybicki & Lightman (see Eq 6.31c)
        #     F(x) ≡ x ∫_x^∞ K_{5/3}(ξ) dξ              (6.31c)
        for j in 1:n_photon_synch
            if bmag_curr < 1e-20 || ω_c < 1e-55
                F = 0.0
            else
                ω_γ = energy_γ[j] / ħ    # from E = ħω
                x = ω_γ/ω_c

                xxx_max_set = 30.0
                if x ≥ xxx_max_set || x < 1e-15
                    F = 0.0
                else
                    F = x * quadgk(t -> besselk(nu, t), x, Inf)
                end
            end

            # In the MC code, the electron spectra passed to synch_emission contain the
            # total number of particles in each momentum shell. Therefore, the emission
            # returned by this subroutine has units of erg/s. (It would be erg/(s⋅cm³)
            # if we had passed density instead of number.)
            # Also, Eq. (6.18) from Rybicki & Lightman, [xnum_electron * p_fac * F], is
            # energy production rate per frequency, dP/dω. Since ω_γ is E/ħ, dω = dE/ħ, and
            # so ω_γ/dω = E/dE. Then [dP/dω * ω_γ] is equal to [dP/dE * E], or dP/d(lnE).
            # This is what `synch_emis()` expects as output. Has units [erg/s].
            tmp_add = xnum_electron * ω_γ * p_fac * F

            # Only include emission if it's sufficiently positive
            if tmp_add > 1e-55
                synch_emis[j] += tmp_add
            end

        end # loop over n_photon_synch
    end
end
