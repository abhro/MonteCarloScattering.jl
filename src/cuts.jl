using .parameters: na_particles

"""
    new_pcut(...)

Takes the `**_sav` arrays filled over the course of loop_pt and splits the
saved particles to form the population of the next pcut

### Arguments
TODO
- `n_pts_target`
- `n_saved`
- `grid_sav`
- `tcut_sav`
- `l_save`
- `downstream_sav`
- `inj_sav`
- `weight_sav`
- `ptot_pf_sav`
- `pb_pf_sav`
- `x_PT_cm_sav`
- `xn_per_sav`
- `prp_x_cm_sav`
- `acctime_sec_sav`
- `φ_rad_sav`
- `n_pts_use`

### Returns
TODO

### Modifies
- `weight_running`: weight factor of each particle remaining after this pcut
"""
function new_pcut(
        n_pts_target, n_saved, l_save, grid_sav, downstream_sav,
        inj_sav, weight_sav, ptot_pf_sav, pb_pf_sav, x_PT_cm_sav, xn_per_sav,
        prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav,
        n_pts_use, weight_running)

    # Determine multiplicity of splitting; perhaps none needed
    i_mult = max(n_pts_target ÷ n_saved, 1) # In case n_pts_target drops between pcuts

    grid_new        = zeros(Int,        n_pts_use*i_mult)
    tcut_new        = zeros(Int,        n_pts_use*i_mult)
    downstream_new  = zeros(Bool,       n_pts_use*i_mult)
    inj_new         = zeros(Bool,       n_pts_use*i_mult)
    weight_new      = zeros(Float64,    n_pts_use*i_mult)
    ptot_pf_new     = zeros(MomentumCGS, n_pts_use*i_mult)
    pb_pf_new       = zeros(MomentumCGS, n_pts_use*i_mult)
    x_PT_cm_new     = zeros(LengthCGS,  n_pts_use*i_mult)
    xn_per_new      = zeros(Float64,    n_pts_use*i_mult)
    prp_x_cm_new    = zeros(LengthCGS,  n_pts_use*i_mult)
    acctime_sec_new = zeros(TimeCGS,    n_pts_use*i_mult)
    φ_rad_new       = zeros(Float64,    n_pts_use*i_mult)

    # Calculate effect on particle weights and the weighting factor of each
    # remaining particle in the simulation
    # CHECKTHIS: is that last claim still true if old particles are imported
    # into the simulation?
    weight_running /= i_mult

    # Perform the splitting
    n_pts_new = 0

    for j in 1:n_pts_use

        l_save[j] || continue # Don't multiply particles that weren't kept, obviously

        for i in 1:i_mult

            n_pts_new += 1

            weight_new[n_pts_new]      = weight_sav[j] / i_mult
            ptot_pf_new[n_pts_new]     = ptot_pf_sav[j]
            pb_pf_new[n_pts_new]       = pb_pf_sav[j]
            x_PT_cm_new[n_pts_new]     = x_PT_cm_sav[j]
            grid_new[n_pts_new]        = grid_sav[j]
            downstream_new[n_pts_new]  = downstream_sav[j]
            inj_new[n_pts_new]         = inj_sav[j]
            xn_per_new[n_pts_new]      = xn_per_sav[j]
            prp_x_cm_new[n_pts_new]    = prp_x_cm_sav[j]
            acctime_sec_new[n_pts_new] = acctime_sec_sav[j]
            φ_rad_new[n_pts_new]       = φ_rad_sav[j]
            tcut_new[n_pts_new]        = tcut_sav[j]

        end  # loop over splits
    end  # loop over saved particles

    return (grid_new, tcut_new, downstream_new, inj_new, weight_new, ptot_pf_new, pb_pf_new,
            x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new,
            n_pts_new, weight_running)
end

function pcut_finalize(
        i_iter, i_ion, i_cut, p_pcut_hi, n_pts_pcut, n_pts_pcut_hi, n_pts_use,
        weight_running, l_save, t_start, pcuts_in, pcuts_use, outfile)
    break_pcut = false
    n_saved = count(l_save)

    t_end = now()
    run_time = t_end - t_start

    @info("Finalizing pcut", i_iter, i_ion, i_cut,
          pcuts_in[i_cut], "pcuts_use[i_cut]"=NoUnits(pcuts_use[i_cut]/(mp*c)),
          n_saved, n_pts_use, weight_running, run_time)

    # If no particles saved, don't bother with remaining pcuts
    if n_saved == 0
        break_pcut = true
        return break_pcut
    end

    # Prepare population for next pcut
    n_pts_target = pcuts_use[i_cut] < p_pcut_hi ? n_pts_pcut : n_pts_pcut_hi
    return break_pcut, n_pts_target, n_saved
end

"""
    tcut_track!(...)

If directed in data_input, track particles' momenta -- and whether they are
still interacting with the shock -- at various times since acceleration began.

### Modifies

- `weight_coupled`
- `spectra_coupled`

### Arguments

- `tcut_curr`: current tcut for tracking
- `weight`: particle weight
- `ptot_pf`: plasma-frame total momentum
- `i_ion`
- `num_psd_mom_bins`

### Returns

Nothing. All adjustments made to input arrays
"""
function tcut_track!(
        weight_coupled::AbstractArray, spectra_coupled::AbstractArray,
        tcut_curr, weight, ptot_pf,
        i_ion, num_psd_mom_bins, psd_mom_min, psd_bins_per_dec_mom)
    # Since particle is still coupled to shock (i.e. being accelerated),
    # add its weight to appropriate bin of weight_coupled
    weight_coupled[tcut_curr,i_ion] += weight

    # For spectra, need to convert ptot_pf into a psd bin, then add it to array;
    #  note that we don't care about angular component
    i_pt = get_psd_bin_momentum(ptot_pf, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
    spectra_coupled[i_pt, tcut_curr, i_ion] += weight
end
