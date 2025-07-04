using .constants: E₀_proton
using .parameters: energy_rel_pt
using .transformers: transform_p_PS

const _pf_spike_away = 1000.0

"""
    particle_finish!(...)

Handles particles that leave the system during loop_helix for any reason.

### Arguments
- aa: particle atomic mass
- pb_pf: component of ptot_pf parallel to magnetic field
- p_perp_b_pf: component of ptot_pf perpendicular to magnetic field
- γₚ_pf: Lorentz factor associated with ptot_pf
- φ_rad: phase angle of gyration; looking UpS, counts clockwise from +z axis
- uₓ_sk: bulk flow speed along x axis
- uz_sk: bulk flow speed along z axis
- utot: total bulk flow speed
- γᵤ_sf: Lorentz factor associated with utot
- b_cosθ: component of magnetic field along x axis
- b_sinθ: component of magnetic field along z axis
- i_reason: integer reason for why particle left system
- weight: particle's weight

### Modifies
TODO

### Returns
Nothing; modifies input argument arrays as needed.
"""
function particle_finish!(
        pₓ_esc_feb, energy_esc_feb, esc_energy_eff, esc_num_eff,
        esc_flux, esc_psd_feb_DwS, esc_psd_feb_UpS,
        i_reason, i_iter, i_ion,
        num_psd_θ_bins,
        aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad, uₓ_sk, uz_sk, utot, γᵤ_sf,
        b_cosθ, b_sinθ, weight, mc,
    )

    @debug("Input arguments:", aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad, uₓ_sk, uz_sk, utot, γᵤ_sf,
           b_cosθ, b_sinθ, i_reason, weight, esc_psd_feb_DwS, esc_psd_feb_UpS, esc_flux,
           pₓ_esc_feb, energy_esc_feb, esc_energy_eff, esc_num_eff,
           i_iter, i_ion, mc)


    # Transform plasma frame momentum into shock frame for binning
    ptot_sk, p_sk, γₚ_sk = transform_p_PS(
        aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad, uₓ_sk, uz_sk, utot, γᵤ_sf,
        b_cosθ, b_sinθ, mc)
    @debug("Values from transform_p_PS:", ptot_sk, p_sk, γₚ_sk)

    # Get PSD bins for this particle
    ip = get_psd_bin_momentum(ptot_sk, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
    jθ = get_psd_bin_angle(p_sk.x, ptot_sk, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)

    if ptot_sk > abs(_pf_spike_away*p_sk.x)
        weight_factor = γₚ_sk * aa*mp * _pf_spike_away / ptot_sk
    else
        weight_factor = γₚ_sk * aa*mp / abs(p_sk.x)
    end

    # Now take additional action based on *how* the particle left the grid
    if i_reason == 1     # Particle escape: DwS, with or without scattering enabled

        @info "About to increase this boy" weight weight_factor
        esc_psd_feb_DwS[ip, jθ] += weight * weight_factor

    elseif i_reason == 2 # Particle escape: pmax, UpS FEB, transverse distance

        esc_flux[i_ion] += weight
        esc_psd_feb_UpS[ip, jθ] += weight * weight_factor

        if (γₚ_sk - 1) < (energy_rel_pt/(aa*E₀_proton))
            energy_flux_add = ptot_sk^2 / (2 * aa*mp) * weight
        else
            energy_flux_add = (γₚ_sk - 1) * aa*E₀_proton * weight
        end

        # Update escape arrays that will be averaged over consecutive iterations
        pₓ_esc_feb[i_ion, i_iter] += abs(p_sk.x) * weight
        energy_esc_feb[i_ion, i_iter] += energy_flux_add

        # Update escape arrays to be printed out with spectral information
        esc_energy_eff[ip, i_ion] += energy_flux_add
        esc_num_eff[ip, i_ion] += weight

    elseif i_reason == 3 # Particle escape: age_max

        # TODO: write out particles to be read in during a later run, i.e. a
        # pre-existing population of CRs. Also include new keyword for this purpose

    elseif i_reason == 4 # Zero energy after radiative losses

        # Do nothing

    else
        error("Unknown i_reason passed: $i_reason. Can only handle 1-4")
    end
end
