using Unitful, UnitfulAstro
using UnitfulAstro: GeV, erg
using .constants: mₚ_cgs, E₀_proton, T_th
using .KATV2014: get_σ_π, get_Ffunc, get_Amax

"""
Calculates pion decay emission from given spectrum of particles interacting
with a medium of density ds_proton_density.

See: Kafexhiu et al., PhRvD 2014, V.90, 123014

### Arguments

- num_hist_bins: number of momentum bins in the distribution of thermal particles
- p_pf_cgs_therm: momentum boundary values, cgs units, of thermal distribution histogram
- dN_therm: thermal particle distribution. It is a pure number of particles in the spectral region
- num_psd_mom_bins: number of momentum bins in the distribution of accelerated particles
- p_pf_cgs_cr: momentum boundary values, cgs units, of cosmic ray distribution histogram
- dN_cr: cosmic ray distribution. It is a pure number of particles in the spectral region
- n_photon_pion: number of energy bins to use for photon production. Must be an even number.
- target_density: density of "thermal" protons with which the particle spectrum will interact
- ID: which of the various decay products to calculate the spectrum of. Currently only photons (ID=1) enabled.
- photon_pion_min_MeV: minimum photon energy, in MeV, to use for pion decay spectrum.
- bins_per_dec_photon: number of energy bins per decade of photon spectrum
- aa: mass number of the nucleus being handled

### Returns

- energy_γ_cgs: energy bin values for pion decay emission spectrum
- pion_emis: pion decay emission spectrum
"""
function pion_kafexhiu(
        num_hist_bins, p_pf_cgs_therm, dN_therm,
        num_psd_mom_bins, p_pf_cgs_cr, dN_cr, n_photon_pion,
        target_density, ID, photon_pion_min_MeV, bins_per_dec_photon, aa,
        n_ions, aa_ion, ρ_N₀_ion, mc
    )


    # ID != 1 is yet to be implemented
    ID != 1 && error("ID must be 1 (photons), ID = $ID is not yet implemented")


    # Calculate the scaling due to the presence of nuclei heavier than hydrogen.
    # This uses Eq. (26) from Baring et al. (1999) [1999ApJ...513..311B].
    #     σ ~ (A_{cr}^⅜ + A_{ISM}^⅜ - 1)² σ_{pp→π⁰X}........(Eq 26, Baring 1999)
    # Note that, per Orth & Buffington (1976), "the error incurred from
    # using [both aa's] greater than 1 is unknown".
    # The target density passed to the subroutine is number of *protons* per
    # cubic centimeter. This must be modified as well:
    #   n_targ*σ -->  n_p*σ_Xp + n_He*σ_X-He + ...
    #         = n_p * [n_p/n_p * σ_Xp + n_He/n_p * σ_X-He ...]
    #         = n_p * σ_pp * [ den1_in[1]*SF_Xp + den1_in[2]*SF_X-He ...]
    scaling_factor = 0.0
    for i in 1:n_ions
        if aa_ion[i] ≥ 1
            scaling_factor += (aa^0.375 + aa_ion[i]^0.375 - 1)^2 * ρ_N₀_ion[i]/ρ_N₀_ion[1]
        end
    end

    # Set parameters that affect energy_γ_cgs and "zero" out emission prior to the calculation
    γ_min_log = log10(photon_pion_min_MeV)
    energy_γ_cgs = ustrip(u"erg", exp10.(range(start = γ_min_log,
                                               step = 1/bins_per_dec_photon,
                                               length = n_photon_pion)) * u"MeV")
    pion_emis = fill(1e-99, n_photon_pion)

    # "1" to use GEANT 4 data
    # "2" to use PYTHIA 8 data
    # "3" to use SIBYLL 2.1 data
    # "4" to use QGSJET-I data
    i_data = 1

    if i_data < 1 || i_data > 4
        throw(ArgumentError("Invalid selection for cross-section data. i_data must be between 1 and 4, not $i_data"))
    end
    #-------------------------------------------------------------------------
    # Constants fixed


    # Loop over the number of bins in the distribution function, calculating
    # the total pion emission from each. First up, the thermal particles
    #-------------------------------------------------------------------------
    for i_fp in 0:num_hist_bins-1

        bin_count = dN_therm[i_fp]
        bin_count ≤ 1e-99 && continue # skip empty bins

        p_pf_sq = p_pf_cgs_therm[i_fp] * p_pf_cgs_therm[i_fp+1] # Geometric mean
        γ = √(p_pf_sq/mc^2 + 1)
        Tp = (γ - 1) * aa*ustrip(GeV, E₀_proton*erg) # particle kinetic energy in GeV
        Tp /= aa  # kinetic energy per nucleon
        vel = √p_pf_sq / (γ*aa*mₚ_cgs)

        # Tp must be at T_th; otherwise no possibility to produce pions/photons
        Tp < T_th && continue

        # Square of proton energy in center-of-mass frame will be used repeatedly
        s_ECM = 2rmp * (Tp + 2rmp)

        # Calculate inclusive pion production cross section using material from section 4
        σ_π = get_σ_π(Tp, i_data, s_ECM)

        # γ-ray production is parametrized as Amax(Tp) * F(Tp, Eγ). Since Amax
        # doesn't depend on photon energy, calculate it here; note that Eγ_max,
        # the maximum photon energy allowed for this value of Tp, is also an output
        Eγ_max, Amax = get_Amax(Tp, i_data, s_ECM, σ_π)


        # Now, loop over photon energies. Make sure that kinematic limits are
        # respected, i.e. that Xγ parameter of Equation (12) falls between 0 and 1.
        #-----------------------------------------------------------------------
        for i_γ in 1:n_photon_pion

            Eγ = ustrip(u"GeV", energy_γ_cgs[i_γ]*u"erg")  # in GeV

            # Calculate F function for current value of Tp and Eγ
            F_func = get_Ffunc(Tp, Eγ, i_data, Eγ_max)

            # Per Equation (8), differential cross section is Amax * F_func. Multiply by Eγ
            # (in consistent units!) to go from differential cross-section to differential
            # cross section per log energy bin, i.e.,
            #     dσ/dE ⋅ E = dσ/d(lnE).
            # This latter is expected by the rest of the subroutine. Note that dσ/d(lnE) has
            # a unitless denominator; it doesn't matter whether the energy used for
            # calculating it was GeV, ergs, or something else entirely
            σ_total = Amax * F_func * Eγ

            # σ_total is the differential cross section per logarithmic energy bin, i.e.
            # dσ/d(lnE). Convert to number production rate per log energy bin, d²N/d(lnE)dt,
            # by multiplying by target density, number of primaries, and velocity of
            # primaries. (Convert cross section from millibarns to cm² in the process.)
            # pion_emis_photon has units of photons/sec, since bin_count is a number of
            # particles rather than a density.
            pion_emis_photon = target_density * bin_count * vel * (σ_total*1e-27)

            # As calculated previously, pion_emis_photon is d²N/d(lnE)dt. Multiplying by energy
            # makes it an energy production rate (power) per logarithmic energy bin, dP/d(lnE).
            # This is what photon_pion_decay expects as output, with units of [erg/sec].
            pion_emis[i_γ] += pion_emis_photon*energy_γ_cgs[i_γ]

        end
        #-----------------------------------------------------------------------
        # Loop over photons

    end # loop over distribution function momentum bins
    #--------------------------------------------------------------------------


    # With the thermal particles finished, move on to the cosmic ray population
    #--------------------------------------------------------------------------
    for i_fp in 0:num_psd_mom_bins

        bin_count = dN_cr[i_fp]
        bin_count ≤ 1e-99 && continue # skip empty bins

        p_pf_sq = p_pf_cgs_cr[i_fp] * p_pf_cgs_cr[i_fp+1] # Geometric mean
        γ = √(p_pf_sq/mc^2 + 1)
        Tp  = (γ - 1) * aa*ustrip(GeV, E₀_proton*erg) # particle K.E. in GeV
        Tp  = Tp / aa  # kinetic energy per nucleon
        vel = √(p_pf_sq)/(γ*aa*mₚ_cgs)

        # Tp must be at T_th; otherwise no possibility to produce pions/photons
        Tp < T_th && continue

        # Square of proton energy in center-of-mass frame will be used repeatedly
        s_ECM = 2rmp * (Tp + 2rmp)

        # Calculate inclusive pion production cross section using material from section 4
        σ_π = get_σ_π(Tp, i_data, s_ECM)

        # γ-ray production is parametrized as Amax(Tp)⋅F(Tp, Eγ). Since Amax doesn't depend
        # on photon energy, calculate it here; note that Eγ_max, the maximum photon energy
        # allowed for this value of Tp, is also an output
        Eγ_max, Amax = get_Amax(Tp, i_data, s_ECM, σ_π)


        # Now, loop over photon energies. Make sure that kinematic limits are respected,
        # i.e. that Xγ parameter of Equation (12) falls between 0 and 1.
        #-----------------------------------------------------------------------
        for i_γ in 1:n_photon_pion

            Eγ = ustrip(GeV, energy_γ_cgs[i_γ]*erg)  # in GeV

            # Calculate F function for current value of Tp and Eγ
            F_func = get_Ffunc(Tp, Eγ, i_data, Eγ_max)

            # Per Equation (8), differential cross section is Amax * F_func. Multiply by Eγ
            # (in consistent units!) to go from differential cross-section to differential
            # cross section per log energy bin, i.e.,
            #     dσ/dE * E  =  dσ/d(lnE).
            # This latter is expected by the rest of the subroutine. Note that dσ/d(lnE) has
            # a unitless denominator; it doesn't matter whether the energy used for
            # calculating it was GeV, ergs, or something else entirely
            σ_total = Amax * F_func * Eγ

            # σ_total is the differential cross section per logarithmic energy bin, i.e.
            # dσ/d(lnE). Convert to number production rate per log energy bin, d²N/d(lnE)dt,
            # by multiplying by target density, number of primaries, and velocity of
            # primaries. (Convert cross section from millibarns to cm² in the process.)
            # pion_emis_photon has units of photons/sec, since bin_count is a number of
            # particles rather than a density.
            pion_emis_photon = target_density * bin_count * vel * (σ_total*1e-27)


            # As calculated previously, pion_emis_photon is d²N/d(lnE)dt.
            # Multiplying by energy makes it an energy production rate (power)
            # per logarithmic energy bin, dP/d(lnE). This is what
            # photon_pion_decay expects as output, with units of [erg/sec].
            pion_emis[i_γ] += pion_emis_photon*energy_γ_cgs[i_γ]

        end
        #-----------------------------------------------------------------------
        # Loop over photons


    end # loop over distribution function momentum bins
    #---------------------------------------------------------------------------
    # Finished


    # Put floor on emission for plotting purposes
    for i_γ in 1:n_photon_pion
        if pion_emis[i_γ] < 1e-99
            pion_emis[i_γ] = 1e-99
        else # Incorporate the scaling factor calculated at the beginning of the subroutine
            pion_emis[i_γ] *= scaling_factor
        end
    end

    return energy_γ_cgs, pion_emis
end
