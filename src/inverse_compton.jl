using OffsetArrays: OffsetMatrix
using Unitful: ustrip, erg, MeV, c, k as kB, h
using UnitfulGaussian: qcgs
include("constants.jl")
using .constants: E₀ₚ, E₀ₑ, T_CMB0
include("parameters.jl")
using .parameters: na_photons, E_rel_pt
include("io.jl")

"""
    photon_IC(...)

Calculates photon production by inverse Compton emission for a provided distribution of electrons.

While the subroutine can handle any number of photon distributions, and can readily be
changed to read in an arbitrary distribution (e.g. the photon fields calculated in
photon_pion_decay and photon_synch), right now the subroutine is set up to only work with
one photon field: the CMB.

### Arguments
- `n_grid`: current grid zone. Used for tracking emission output
- `p_pf_cr`: momentum boundary values, cgs units, of cosmic ray distribution histogram
- `num_psd_mom_bins`: number of momentum bins in the distribution of accelerated particles
- `cos_bounds`: boundaries of the angular bins of PSD, but in cosine form rather than mixed
  cosine and θ form (as in psd_θ_bounds)
- `num_psd_θ_bins`: number of angular bins in the distribution array
- `d²Ndp_slice`: particle distribution array. It is number of particles per dp (NOT #/cm³/dp),
  split by pitch angle (NOT per pitch angle)
- `n_photon_IC`: number of energy bins to use for photon production
- `photon_ic_min_MeV`: minimum photon energy, in MeV, to use for inverse Compton spectrum.
  Passed to `IC_emission_FCJ`.
- `bins_per_dec_photon`: number of energy bins per decade of photon spectrum. Passed to `IC_emission_FCJ`.
- `dist_lum`: luminosity distance (i.e. including redshift correction) to source; passed to `IC_emission_FCJ`.
- `redshift`: redshift of source, used to adjust photon energies; passed to `IC_emission_FCJ`.
"""
function photon_IC(
        n_grid, p_pf_cr, num_psd_mom_bins, cos_bounds,
        num_psd_θ_bins, d²Ndp_slice, n_photon_IC, photon_ic_min_MeV,
        bins_per_dec_photon, dist_lum, redshift,
        n_IC_specs, energy_IC_MeV, ic_photon_sum,
        dp_pf_cr, jet_sph_frac, mc
    )

    energy_γ_MeV = zeros(na_photons)

    # Our distribution function has already been normalized to the total number of emitting
    # particles, so no additional scaling is needed. However, it must be converted from
    # particles per momentum into pure particle count. Note the ordering of the dimensions
    # in d²N_slice:
    #    1  -  angle          2  -  momentum
    d²N_slice = OffsetMatrix{Float64}(undef, 0:num_psd_mom_bins, 0:num_psd_θ_bins)
    dp = diff(p_pf_cr)
    for i in 0:num_psd_mom_bins, j in 0:num_psd_θ_bins
        if d²Ndp_slice[j,i] ≤ 1e-99
            d²N_slice[j,i] = 1e-99
        else
            d²N_slice[j,i] = d²Ndp_slice[j,i] * dp_pf_cr[i]
        end
    end

    local j_unit

    # Loop over the photon fields, creating an output file for each field's
    # contribution to IC emission.
    # Right now, there's just one photon field: the CMB.
    # TODO: have access to synchrotron photon field. Use it.
    # TODO: investigate interplay between j3 and n_IC_specs if more than
    #  one field is used. Not sure current code is right
    for j3 in 1:1# i_photon_fields

        # Set number of inverse Compton spectra here. Will be 1 unless more than one photon
        # field was used. Also, zero out the column of ic_photon_sum that would hold the
        # summed emission if more than one photon field is used.
        n_IC_specs = j3
        if j3 == 1
            ic_photon_sum[:,n_grid] = 1e-99
            energy_IC_MeV .= 1e-99
        end

        # Compute the spectra
        energy_γ, ic_emis = IC_emission_FCJ(
            num_psd_mom_bins, p_pf_cr,
            num_psd_θ_bins, cos_bounds, d²N_slice, n_photon_IC, j3,
            photon_ic_min_MeV, bins_per_dec_photon, dist_lum, redshift,
            jet_sph_frac, mc)

        # Convert units of energy_γ and ic_emis
        # NOTE: unlike the other two photon production subroutines, ic_emis comes out of
        # IC_emission_FCJ already in the form of observed energy flux at Earth per
        # log energy bin, i.e. dP/(d(lnE)-dA). Pion production and synchrotron both require
        # processing from dP/d(lnE). The units of ic_emis are [erg/s⋅cm²].
        for i in 1:n_photon_IC
            energy_γ_MeV[i] = ustrip(MeV, energy_γ[i]*erg)

            # Add the current photon flux to the running total over all fields
            if ic_emis[i] > 1e-99
                ic_photon_sum[i,n_grid] += ic_emis[i] / energy_γ[i]
            end
            energy_IC_MeV[i] = energy_γ_MeV[i]

        end   # units on emis_γ: erg/(cm²⋅s)

        # Don't write out anything if ic_emis is empty
        count(ic_emis .> 1e-99) < 1 && continue

        # Do necessary unit conversion and write out results to file
        iplot = 0
        for i in 1:n_photon_IC
            iplot += 1

            # Above is energy flux [MeV/(cm²⋅s)] per log energy bin d(lnE) = dE/E.
            # Energy in spectrum is area under curve when plotted with a
            # logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if ic_emis[i] > 1e-99
                emis_γ_MeV = ustrip(MeV, ic_emis[i]*erg)  # MeV/(cm²⋅s) at earth
            else
                emis_γ_MeV = 1e-99               # "zero" emission
            end

            # This is photon flux [#/(cm²⋅sec)] per log energy bin d(lnE) = dE/E.
            # Number of photons in spectrum is area under curve when plotted
            # with a logarithmic energy axis, i.e. [dΦ/d(lnE) * d(lnE)].
            if emis_γ_MeV ≤ 1e-99
                photon_flux = 1e-99                  # "zero" emission
            else
                photon_flux = emis_γ_MeV/energy_γ_MeV[i]
            end

            # Open the file to which we will write the spectral data.
            if j3 == 1
                lopen = inquire(:isopen, file="./photon_IC_grid.dat")
                if !lopen
                    j_unit = open(status="unknown", file="./photon_IC_grid.dat")
                end
            end

            # Don't write out the last data point
            i == n_photon_IC && continue

            ##TODO: incorporate redshift into this write-out;
            # right now, everything is handled in time_seq_photons so this section is irrelevant
            write(j_unit, n_grid, iplot, # photon_IC_grid.dat
                  j3,                                   # 1 photon source
                  log10(photon_flux),                   # 2 log10(photons/(cm²⋅s))
                  log10(energy_γ_MeV[i]),               # 3 log10(MeV)
                  log10(emis_γ_MeV),                    # 4 log[MeV/(cm²⋅s)] at earth
                  log10(photon_flux/energy_γ_MeV[i]))   # 5 log10[photons/(cm²⋅s⋅MeV)]

        end # loop over n_photon_IC

        print_plot_vals(j_unit)

    end # Loop over photon fields

    close(j_unit)
end

"""
This subroutine uses the derivation of Jones (1968) to determine the inverse Compton emission
from a given electron distribution upscattering a specified photon field.

### Arguments

- `num_psd_mom_bins`: number of momentum bins in the distribution of particles
- `p_pf_cr`: momentum boundary values, cgs units, of cosmic ray distribution histogram
- `num_psd_θ_bins`: number of angular bins in the distribution array
- `cos_bounds`: boundaries of the angular cosine bins for determining CR pitch angle
- `d²N_slice`: particle distribution array. It is pure number of particles
  in each dp (NOT per dp), split by pitch angle (NOT per pitch angle)
- `n_photon_IC`: number of energy bins to use for photon production
- `j3`: which photon field to use (e.g. CMB, synchrotron photons, etc.)
- `photon_ic_min_MeV`: minimum photon energy, in MeV, to use for IC spectrum
- `bins_per_dec_photo`n: number of energy bins per decade of photon spectrum
- `dist_lum`: luminosity distance (i.e. including redshift correction) to source
- `redshift`: redshift of source, used to adjust photon energies and CMB

### Returns

- `energy_γ`: energy bin values of resultant photon distribution
- `ic_emis`: observed IC spectrum at Earth, units of erg/(sec⋅cm²). Note that this is NOT the
  same units as produced by the other two photon production subroutines. They do not
  convert to flux, while IC_emission_FCJ does.
"""
function IC_emission_FCJ(
        num_psd_mom_bins, p_pf_cr,
        num_psd_θ_bins, cos_bounds, d²N_slice, n_photon_IC, j3,
        photon_ic_min_MeV, bins_per_dec_photon, dist_lum, redshift,
        jet_sph_frac, mc,
    )

    # Constants that will be used during the calculation
    #----------------------------------------------------------------------------
    photon_out_erg_min = ustrip(erg, photon_ic_min_MeV*MeV) # Minimum photon energy in erg

    photon_out_min_rm  = photon_out_erg_min/E₀ₑ  # minimum en/(mₑc²)
    photon_out_min_log = log10(photon_out_min_rm)
    Δphoton_out        = 1 / bins_per_dec_photon

    # α_out here is outgoing photon energy in units of electron rest mass;
    # calculated here to save time during loops to follow
    α_out = exp10.(range(start = photon_out_min_log, step = Δphoton_out, length = n_photon_IC))

    # Calculate the pitch angle bin that corresponds to the opening angle of the jet; only
    # electrons whose pitch angle is *less* than that can emit photons that are detected at
    # Earth. Remember that cos_bounds is defined so that negative pitch angles point upstream,
    # so need to negate it when counting from jet axis.
    jθ_max = findfirst(>(2jet_sph_frac - 1), cos_bounds)

    # for Eq.(9), Frank Jones, PhysRev. 1968, V.167, p.1159
    r_z = qcgs^2/E₀ₑ
    #-------------------------------------------------------------------------
    # End constants section

    # Set various quantities related to the photon field:
    # - photon_temp_IC is the effective temperature of the field
    # - energy_density_IC is the energy density of the field, RELATIVE to the CMB
    # Note that if you are reading in the photon spectrum from a file (not an
    # option at present) these settings are irrelevant
    #----------------------------------------------------------------------------
    if j3 == 1 # the CMB
        local photon_temp_IC = T_CMB0 * (1 + redshift)
        local energy_density_IC = 1.0
    end

    # Wien displacement law in freq (Hz) (Wikipedia)
    freq_peak_Hz = 5.879e10 * photon_temp_IC

    # More constants related to the photon field
    n_freq       = 60   # 100 maximum
    freq_min     = freq_peak_Hz/30
    freq_max     = freq_peak_Hz*20   ##*10
    freq_min_log = log10(freq_min)
    Δfreq        = (log10(freq_max) - freq_min_log)/n_freq

    ∑xnum_photons          = 0.0
    ∑energy_photons        = 0.0
    ∑photon_energy_density = 0.0

    photon_energy_rm   = Vector{Float64}(undef, na_photons)
    xnum_photons_p_vol = Vector{Float64}(undef, na_photons)
    # Now actually populate the photon distribution. If j3 = -1, we are reading
    # in a distribution from a file that should be specified here. This
    # functionality has not yet been enabled, though, so don't have j3 = -1!
    if j3 == -1
        #for i in 1:n_freq
        #    # read in the variables photon_seed_energy and photon_seed_density from file,
        #    # and then use those to set the five variables below
        #
        #    photon_energy_rm[i] = photon_seed_energy[i]*eV/E₀ₑ
        #    xnum_photons_p_vol[i] = photon_seed_density[i]/photon_seed_energy[i]
        #
        #    ∑energy_photons += photon_seed_density[i]*eV # erg/cm³
        #    ∑xnum_photons += xnum_photons_p_vol[i]
        #    ∑photon_energy_density += photon_seed_density[i]*eV #erg/cm³
        #end
    else
        # Below are constants for spectral energy density of CMB in units of
        # energy per volume per Hz taken from Wikipedia
        con_f1 = 8π*h/(c^3)
        con_f2 = h/(kB*photon_temp_IC)

        for j_in in 1:n_freq     # loop over incoming photon energy density
            f1 = exp10(freq_min_log + (j_in-1)*Δfreq) # = freq_min (freq_max/freq_min)^[(j-1)/n_freq]
            f2 = exp10(freq_min_log + (j_in)*Δfreq)   # = freq_min (freq_max/freq_min)^[ j   /n_freq]
            f_avg = √(f1*f2) # = freq_min ⋅ (freq_max/freq_min)*[(2j-1)/2n_freq]
            exp_fac = exp(min(con_f2*f_avg, 200.0))

            # Below is incoming photon energy density derived from CMB energy density
            photon_energy_density   = (f2 - f1) * con_f1 * f_avg^3 / (exp_fac - 1)
            ∑photon_energy_density += photon_energy_density

            photon_energy_erg      = h*f_avg      # incoming photon energy (erg)
            photon_energy_rm[j_in] = photon_energy_erg/E₀ₑ

            # Below is number density of incoming photons [/cm³] in frequency bin
            xnum_photons_p_vol[j_in] = photon_energy_density/photon_energy_erg
            ∑xnum_photons   += xnum_photons_p_vol[j_in]
            ∑energy_photons += photon_energy_erg*xnum_photons_p_vol[j_in]
        end
    end
    #-------------------------------------------------------------------------
    # End photon field section


    # Now loop over momenta (i.e. energy) of electrons and calculate the inverse Compton emission.
    #-------------------------------------------------------------------------
    d²N_o_dtda = fill(1e-99, na_photons) # initialize output spectrum

    for i_el in 0:num_psd_mom_bins

        # Skip empty rows of PSD, as each one requires significant computation
        maximum(d²N_slice[begin:jθ_max,i_el]) ≤ 1e-99 && continue

        # Otherwise, sum number of electrons at this energy, since that's all equation (9) needs
        xnum_electron = sum(d²N_slice[begin:jθ_max,i_el])

        p1 = √(p_pf_cr[i_el] * p_pf_cr[i_el+1])
        γ = p1/mc < E_rel_pt ? 1.0 : hypot(p1/mc, 1)

        # Loop over incoming photons, then over outgoing photons, then over angle.
        # Fill d²N_o_dtdα as defined by equation (9) in the process.
        # Before collision, α₁ is incoming photon energy in units of the electron rest mass.
        # NOTE: right now, this loop is for a *single* electron. The normalization will be handled later.
        for j_in in 1:n_freq

            # for Eq.(9), Jones, PhysRev. 1968, V.167, p.1159
            #   d²N/dtdα ≈ 2πr₀²c/α₁γ² [2q″ ln q″ + (1+2q″)(1−q″) + ½(1-q″)(4α₁γq″)²/(1+4α₁γq″)]
            # where q″ = α/4α₁γ²(1−α/γ) and 1/4γ² < q″ ≤ 1
            α₁       = photon_energy_rm[j_in]
            norm_fac = xnum_photons_p_vol[j_in] * 2π* r_z^2 * c / (α₁ * γ^2)

            # Loop over outgoing photons. After collision, "α_out" is α in Eq. (9):
            # outgoing photon energy in units of electron rest mass
            for k in eachindex(α_out)

                # If the outgoing photon would have more energy than
                # the pre-scatter electron has energy, skip this zone.
                α_out[k] ≥ γ && continue

                # Determine q″ in equation (9)
                q″ = α_out[k] / (4α₁*γ^2 * (1 - α_out[k] / γ))

                # Equation (9) in Jones (1968)
                d²N_cur = norm_fac * xnum_electron * (2q″*log(q″) + (1+2q″)*(1-q″)
                                                      + 8*(α₁*γ*q″)^2*(1-q″)/(1+4α₁*γ*q″))

                # Only include photon production if it's sufficiently positive;
                # otherwise, assume it's zero
                if d²N_cur > 1e-60
                    d²N_o_dtda[k] += d²N_cur
                end
            end
            #-----------------------------------------------------------------
            # loop over outgoing photon energies

        end
        #---------------------------------------------------------------------
        # loop over incoming photon energies

    end # loop over electron energies
    #-------------------------------------------------------------------------
    # Finished!


    # Now have d²N/(dt dα), but need to convert that into energy flux observed at Earth.
    # Assume that the electrons are distributed homogeneously over the surface of the jet,
    # so that the photons are isotropic within that cone. NOTE: in the MC code, the electron
    # spectra passed to the emission subroutines contain the total number of particles in
    # each momentum bin, rather than number per unit momentum. Additionally, we're
    # converting from photons released per second to observed flux. Therefore, the emission
    # returned by this subroutine has units of erg/s⋅cm²
    beam_area = 4π * dist_lum^2 * jet_sph_frac
    # d²N/(dt dα) comes out of Eq.(9) in Jones 1968. This is (number of photons emitted)/(s⋅dα),
    # but needs to be converted to a flux so that it's actually (photons observed) / (s⋅dα⋅cm²).
    d²N_o_dtda ./= beam_area

    energy_γ = α_out * E₀ₑ   # photon energy in ergs
    # To get d²N/(dt dE) we multiply d²N/(dt dα) by dα/dE = 1 / mₑc²
    # photon_IC() expects energy production rate per logarithmic energy bin, dP/d(lnE).
    # Multiply d²N/(dt dE) by E once to make it dP/dE, then again to make it dP/d(lnE).
    ic_emis = @. d²N_o_dtda/E₀ₑ * energy_γ^2  # erg/(s⋅cm²)
    for k in eachindex(ic_emis)
        if ic_emis[k] ≤ 1e-55
            ic_emis[k] = 1e-99
        end
    end

    return energy_γ, ic_emis
end
