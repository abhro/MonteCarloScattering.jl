function clear_psd!(num_crossings, therm_grid, therm_pₓ_sk, therm_ptot_sk,
                    therm_weight, psd, esc_psd_feb_upstream, esc_psd_feb_downstream)
    # zero out each array in a type stable manner
    zero!(num_crossings)
    zero!(therm_grid)
    zero!(therm_pₓ_sk)
    zero!(therm_ptot_sk)
    zero!(therm_weight)
    fill!(psd,                      1e-99)
    fill!(esc_psd_feb_upstream,     1e-99)
    fill!(esc_psd_feb_downstream,   1e-99)
    n_cr_count = 0
    return n_cr_count
end

function assign_slice!(recv_arrays::NTuple{N,AbstractArray}, send_arrays::NTuple{N,AbstractArray}, indices) where N
    for (receiver, sender) in zip(recv_arrays, send_arrays)
        assign_slice!(receiver, sender, indices)
    end
end
function assign_slice!(recv_array::AbstractArray{T}, send_array::AbstractArray{T}, indices) where T
    recv_array[indices] .= view(send_array, indices)
end

function assign_particle_properties_to_population!(
        rng, n_pts_use, xn_per_fine, x_grid_stop,
        weight_new, weight_in, ptot_pf_new, ptot_pf_in,
        pb_pf_new, pb_pf_in, x_PT_cm_new, x_PT_cm_in, grid_new, i_grid_in,
        downstream_new, inj_new, xn_per_new, prp_x_cm_new, acctime_sec_new, tcut_new, φ_rad_new)

    slice = 1:n_pts_use

    assign_slice!((weight_new, ptot_pf_new, pb_pf_new, x_PT_cm_new, grid_new),
                  (weight_in, ptot_pf_in, pb_pf_in, x_PT_cm_in, i_grid_in),
                  slice)

    @inbounds downstream_new[slice]  .= false
    @inbounds inj_new[slice]         .= false
    @inbounds xn_per_new[slice]      .= xn_per_fine
    @inbounds prp_x_cm_new[slice]    .= x_grid_stop
    @inbounds acctime_sec_new[slice] .= 0.0s
    @inbounds tcut_new[slice]        .= 1

    @inbounds φ_rad_new[slice] .= 2π*Random.rand(rng, length(slice))
end

function get_pmax_cutoff(Emax, Emax_per_aa, pmax, aa)
    if Emax > 0keV
        γ = 1 + Emax/(aa*E₀ₚ)
        pmax_cutoff = aa*mp * c * √(γ^2 - 1)
    elseif Emax_per_aa > 0keV
        γ = 1 + Emax_per_aa/E₀ₚ
        pmax_cutoff = aa*mp * c * √(γ^2 - 1)
    elseif pmax > 0g*cm/s
        pmax_cutoff = pmax
    else
        # Something has gone very wrong.
        error("Max CR energy not set in data_input, so can't set pmax_cutoff.")
    end

    return pmax_cutoff
end

function pcut_hi(energy_pcut_hi, E_rel_pt, m)
    E_pcut_hi_rmproton = energy_pcut_hi*keV / E₀ₚ # FIXME pick better name
    if E_pcut_hi_rmproton < E_rel_pt
        p_pcut_hi = √(2 * E_pcut_hi_rmproton)
    else
        p_pcut_hi = m*c * √((E_pcut_hi_rmproton + 1)^2 - 1)
    end
    return p_pcut_hi
end
