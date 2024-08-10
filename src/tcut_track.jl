"""
If directed in data_input, track particles' momenta -- and whether they are
still interacting with the shock -- at various times since acceleration began.

### Arguments

- `tcut_curr`: current tcut for tracking
- `weight`: particle weight
- `pt_pf`: plasma-frame total momentum

### Returns

None. All adjustments made to input arrays
"""
function tcut_track!(tcut_curr, weight, pt_pf, wt_coupled::AbstractArray, spectra_coupled::AbstractArray, i_ion)
    # Since particle is still coupled to shock (i.e. being accelerated),
    # add its weight to appropriate bin of wt_coupled
    wt_coupled[tcut_curr,i_ion] += weight

    # For spectra, need to convert pt_pf into a psd bin, then add it to array;
    #  note that we don't care about angular component
    i_pt = get_psd_bin_momentum(pt_pf, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)

    spectra_coupled[i_pt, tcut_curr, i_ion] += weight
end
