using .constants: MeV2erg, E₀_proton, qp_cgs, T_CMB0, h_cgs, c_cgs, kB_cgs, E₀_electron
using .parameters: psd_max, na_photons, energy_rel_pt, na_ions

"""
This subroutine uses the derivation of Jones (1968) to determine the inverse Compton emission
from a given electron distribution upscattering a specified photon field.

### Arguments

- num_psd_mom_bins: number of momentum bins in the distribution of particles
- p_pf_cgs_cr: momentum boundary values, cgs units, of cosmic ray distribution histogram
- num_psd_θ_bins: number of angular bins in the distribution array
- cos_bounds: boundaries of the angular cosine bins for determining CR pitch angle
- d2N_slice: particle distribution array. It is pure number of particles
  in each dp (NOT per dp), split by pitch angle (NOT per pitch angle)
- n_photon_IC: number of energy bins to use for photon production
- j3: which photon field to use (e.g. CMB, synchrotron photons, etc.)
- photon_ic_min_MeV: minimum photon energy, in MeV, to use for IC spectrum
- bins_per_dec_photon: number of energy bins per decade of photon spectrum
- dist_lum: luminosity distance (i.e. including redshift correction) to source
- redshift: redshift of source, used to adjust photon energies and CMB

### Returns

- energy_γ_cgs: energy bin values of resultant photon distribution
- ic_emis: observed IC spectrum at Earth, units of erg/sec-cm^2. Note that this is NOT the
  same units as produced by the other two photon production subroutines. They do not
  convert to flux, while IC_emission_FCJ does.
"""
function IC_emission_FCJ(
        num_psd_mom_bins, p_pf_cgs_cr,
        num_psd_θ_bins, cos_bounds, d2N_slice, n_photon_IC, j3,
        photon_ic_min_MeV, bins_per_dec_photon, dist_lum, redshift,
        jet_sph_frac, o_o_mc,
    )

    # Constants that will be used during the calculation
    #----------------------------------------------------------------------------
    photon_out_erg_min = photon_ic_min_MeV*MeV2erg # Minimum photon energy in erg

    photon_out_min_rm  = photon_out_erg_min/E₀_electron  # minimum en/(m_e c^2)
    photon_out_min_log = log10(photon_out_min_rm)
    Δphoton_out        = 1 / bins_per_dec_photon

    # α_out here is outgoing photon energy in units of electron rest mass;
    # calculated here to save time during loops to follow
    α_out = exp10.(range(start = photon_out_min_log, step = Δphoton_out, length = n_photon_IC))

    # Calculate the pitch angle bin that corresponds to the opening angle of the jet; only
    # electrons whose pitch angle is *less* than that can emit photons that are detected at
    # Earth. Remember that cos_bounds is defined so that negative pitch angles point upstream,
    # so need to negate it when counting from jet axis.
    cos_θ_max = 1 - 2jet_sph_frac
    jθ_max = findfirst(>(-cos_θ_max), cos_bounds)

    # for Eq.(9), Frank Jones, PhysRev. 1968, V.167, p.1159
    r_z = qp_cgs^2/E₀_electron
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

    xnum_photons_sum          = 0.0
    energy_photons_sum        = 0.0
    photon_energy_density_sum = 0.0

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
        #    photon_energy_rm[i] = (photon_seed_energy[i]/erg2eV)/E₀_electron
        #    xnum_photons_p_vol[i] = photon_seed_density[i]/photon_seed_energy[i]
        #
        #    energy_photons_sum += photon_seed_density[i]/erg2eV # erg/cm^3
        #    xnum_photons_sum += xnum_photons_p_vol[i]
        #    photon_energy_density_sum += photon_seed_density[i]/erg2eV #erg/cm^3
        #end
    else
        # Below are constants for spectral energy density of CMB in units of
        # energy per volume per Hz taken from Wikipedia
        con_f1 = 8π*h_cgs/(c_cgs^3)
        con_f2 = h_cgs/(kB_cgs*photon_temp_IC)

        for j_in in 1:n_freq     # loop over incoming photon energy density
            f1 = exp10(freq_min_log + (j_in-1)*Δfreq) # = freq_min (freq_max/freq_min)^[(j-1)/n_freq]
            f2 = exp10(freq_min_log + (j_in)*Δfreq)   # = freq_min (freq_max/freq_min)^[ j   /n_freq]
            f_avg = √(f1*f2) # = freq_min ⋅ (freq_max/freq_min)*[(2j-1)/2n_freq]
            exp_dum = min(con_f2*f_avg, 200.0)

            # Below is incoming photon energy density derived from CMB en. density
            photon_energy_density = (f2 - f1) * con_f1 * f_avg^3  /  (exp(exp_dum) - 1)
            photon_energy_density_sum = photon_energy_density_sum + photon_energy_density

            photon_energy_erg      = h_cgs*f_avg      # incoming photon energy (erg)
            photon_energy_rm[j_in] = photon_energy_erg/E₀_electron

            # Below is number density of incoming photons [cm^{-3}] in freq. bin
            xnum_photons_p_vol[j_in] = photon_energy_density/photon_energy_erg
            xnum_photons_sum        += xnum_photons_p_vol[j_in]
            energy_photons_sum      += photon_energy_erg*xnum_photons_p_vol[j_in]
        end
    end
    #-------------------------------------------------------------------------
    # End photon field section


    # Now loop over momenta (i.e. energy) of electrons and calculate the inverse Compton emission.
    #-------------------------------------------------------------------------
    d2N_o_dtda = fill(1e-99, na_photons) # initialize output spectrum

    for i_el in 0:num_psd_mom_bins

        # Skip empty rows of PSD, as each one requires significant computation
        maximum(d2N_slice[begin:jθ_max,i_el]) ≤ 1e-99 && continue

        # Otherwise, sum number of electrons at this energy, since that's all equation (9) needs
        xnum_electron = sum(d2N_slice[begin:jθ_max,i_el])


        p1_cgs = √(p_pf_cgs_cr[i_el] * p_pf_cgs_cr[i_el+1])
        γ = p1_cgs*o_o_mc < energy_rel_pt ? 1.0 : √( (p1_cgs * o_o_mc)^2 + 1)

        # Loop over incoming photons, then over outoing photons, then over angle.
        # Fill d2N_o_dtda as defined by equation (9) in the process.
        # Before collision, α_1 is incoming photon energy in units of the electron rest mass.
        # NOTE: right now, this loop is for a *single* electron. The normalization will be handled later.
        for j_in in 1:n_freq

            # for Eq.(9), Jones, PhysRev. 1968, V.167, p.1159
            #   d²N/dtdα ≈ 2πr₀²c/α₁γ² [2q″\ln q″ + (1+2q″)(1−q″) + ½(1-q″)(4α₁γq″)²/(1+4α₁γq″)]
            # where q″ = α/4α₁γ²(1−α/γ) and 1/4γ² < q″ ≤ 1
            α_1      = photon_energy_rm[j_in]
            norm_fac = xnum_photons_p_vol[j_in] * 2π*(r_z^2)*c_cgs / (α_1 * γ^2)

            # Loop over outgoing photons. After collision, "α_out" is α in eq. (9):
            # outgoing photon energy in units of electron rest mass
            for k in eachindex(α_out)

                # If the outgoing photon would have more energy than
                # the pre-scatter electron has energy, skip this zone.
                α_out[k] ≥ γ && continue

                # Determine q″ in equation (9)
                q″ = α_out[k] / (4α_1*γ^2 * (1 - α_out[k] / γ))

                # Equation (9) in Jones (1968)
                d2N_tmp = norm_fac * xnum_electron * ( 2q″*log(q″) + (1 + 2q″)*(1 - q″) + 8 * (α_1*γ*q″)^2 * (1-q″) / (1 + 4*α_1*γ*q″)  )

                # Only include photon production if it's sufficiently positive;
                # otherwise, assume it's zero
                if d2N_tmp > 1e-60
                    d2N_o_dtda[k] += d2N_tmp
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


    # Now have d2N_o_dtda, but need to convert that into energy flux observed at Earth.
    # Assume that the electrons are distributed homogeneously over the surface of the jet,
    # so that the photons are isotropic within that cone. NOTE: in the MC code, the electron
    # spectra passed to the emission subroutines contain the total number of particles in
    # each momentum bin, rather than number per unit momentum. Additionally, we're
    # converting from photons released per second to observed flux. Therefore, the emission
    # returned by this subroutine has units of erg/sec-cm^2
    beam_area = 4π * dist_lum^2 * jet_sph_frac

    energy_γ_cgs = Vector{Float64}(undef, n_photon_IC)
    ic_emis = Vector{Float64}(undef, n_photon_IC)
    for k in 1:n_photon_IC

        # "d2N_o_dtda" comes out of Eq.(9) in Jones 1968.
        # This is number of photons emitted per sec per dα, but needs to be converted
        # to a flux so that it's actually photons observed per sec per dα per cm^2.
        d2N_o_dtda[k] /= beam_area


        # To get "dN_o_dtdE" we multiply dN_o_dtda by da/dE = 1/(m_e c^2)
        # photon_IC() expects energy production rate per logarithmic energy bin, dP/d(lnE).
        # Multiply dN/dtdE by E once to make it dP/dE, then again to make it dP/d(lnE).
        energy_γ_cgs[k] = α_out[k] * E₀_electron   # photon energy in ergs
        ic_emis[k] = d2N_o_dtda[k]/E₀_electron * energy_γ_cgs[k]^2  # erg/(sec-cm^2)

        if ic_emis[k] ≤ 1e-55
            ic_emis[k] = 1e-99
        end
    end

    return energy_γ_cgs, ic_emis
end
