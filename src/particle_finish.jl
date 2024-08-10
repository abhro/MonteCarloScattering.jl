using .constants: mp_cgs, E₀_proton
using .parameters: na_ions, na_itrs, psd_max, energy_rel_pt
using .transformers: transform_p_PS

const _pf_spike_away = 1000.0

"""
Handles particles that leave the system during loop_helix for any reason.

### Arguments
- aa: particle atomic mass
- pb_pf: component of pt_pf parallel to magnetic field
- p_perp_b_pf: component of pt_pf perpendicular to magnetic field
- γ_pt_pf: Lorentz factor associated with pt_pf
- φ_rad: phase angle of gyration; looking UpS, counts clockwise from +z axis
- ux_sk: bulk flow speed along x axis
- uz_sk: bulk flow speed along z axis
- utot: total bulk flow speed
- γ_usf: Lorentz factor associated with utot
- b_cosθ: component of magnetic field along x axis
- b_sinθ: component of magnetic field along z axis
- i_reason: integer reason for why particle left system
- xwt: particle's weight

### Modifies
TODO

### Returns
Nothing; modifies input argument arrays as needed.
"""
function particle_finish!(
        aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad, ux_sk, uz_sk, utot, γ_usf,
        b_cosθ, b_sinθ, i_reason, xwt, esc_psd_feb_DwS, esc_psd_feb_UpS, esc_flux,
        px_esc_feb, energy_esc_feb, esc_energy_eff, esc_num_eff,
        i_iter, i_ion, oblique, o_o_mc)

    @debug("Input arguments:", aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad, ux_sk, uz_sk, utot, γ_usf,
           b_cosθ, b_sinθ, i_reason, xwt, esc_psd_feb_DwS, esc_psd_feb_UpS, esc_flux,
           px_esc_feb, energy_esc_feb, esc_energy_eff, esc_num_eff,
           i_iter, i_ion, oblique, o_o_mc)


    # Transform plasma frame momentum into shock frame for binning
    pt_sk, px_sk, pz_sk, γ_pt_sk = transform_p_PS(
        aa, pb_pf, p_perp_b_pf, γ_pt_pf, φ_rad, ux_sk, uz_sk, utot, γ_usf,
        b_cosθ, b_sinθ, oblique, o_o_mc)
    @debug "Values from transform_p_PS:" pt_sk px_sk pz_sk γ_pt_sk

    # Get PSD bins for this particle
    i_particle = get_psd_bin_momentum(pt_sk, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
    jθ = get_psd_bin_angle(px_sk, pt_sk, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)

    if pt_sk > abs(_pf_spike_away*px_sk)
        wtfac = γ_pt_sk * aa*mp_cgs * _pf_spike_away / pt_sk
    else
        wtfac = γ_pt_sk * aa*mp_cgs / abs(px_sk)
    end

    # Now take additional action based on *how* the particle left the grid
    if i_reason == 1     # Particle escape: DwS, with or without scattering enabled

        esc_psd_feb_DwS[i_particle, jθ] += xwt * wtfac

    elseif i_reason == 2 # Particle escape: pmax, UpS FEB, transverse distance

        esc_flux[i_ion] += xwt
        esc_psd_feb_UpS[i_particle, jθ] += xwt * wtfac

        if (γ_pt_sk - 1) < (energy_rel_pt/(aa*E₀_proton))
            energy_flux_add = pt_sk^2 / (2 * aa*mp_cgs) * xwt
        else
            energy_flux_add = (γ_pt_sk - 1) * aa*E₀_proton * xwt
        end

        # Update escape arrays that will be averaged over consecutive iterations
        px_esc_feb[i_ion, i_iter] += abs(px_sk) * xwt
        energy_esc_feb[i_ion, i_iter] += energy_flux_add

        # Update escape arrays to be printed out with spectral information
        esc_energy_eff[i_particle, i_ion]  += energy_flux_add
        esc_num_eff[i_particle, i_ion] += xwt

    elseif i_reason == 3 # Particle escape: age_max

        # TODO: write out particles to be read in during a later run, i.e. a
        # pre-existing population of CRs. Also include new keyword for this purpose

    elseif i_reason == 4 # Zero energy after radiative losses

        # Do nothing

    else
        error("Unknown i_reason passed: $i_reason. Can only handle 1-4")
    end
end
