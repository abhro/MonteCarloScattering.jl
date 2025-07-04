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
        i_ion, num_psd_mom_bins)
    # Since particle is still coupled to shock (i.e. being accelerated),
    # add its weight to appropriate bin of weight_coupled
    weight_coupled[tcut_curr,i_ion] += weight

    # For spectra, need to convert ptot_pf into a psd bin, then add it to array;
    #  note that we don't care about angular component
    i_pt = get_psd_bin_momentum(ptot_pf, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
    spectra_coupled[i_pt, tcut_curr, i_ion] += weight
end
