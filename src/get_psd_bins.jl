using .constants: mp_cgs, me_cgs, c_cgs

"""
Given parallel component and total of a particle's momentum in the shock frame,
determine which bin of psd particle will fall into.

Binning is done by the value of pt_sk in code units. All binning is logarithmic.

### Arguments
- px_sk: component of momentum parallel to shock normal (to B-field?), in code units
- pt_sk: total particle momentum, in code units

### Returns
Bin in momentum into which particle falls
"""
function get_psd_bin_momentum(pt_sk, psd_bins_per_dec_mom, psd_mom_min, num_psd_mom_bins)
    # Bin in total momentum (bin)
    if pt_sk < psd_mom_min
        # Momentum close enough to 0
        bin = 0
    else
        # Particle falls into logarithmic spacing region
        @debug "" pt_sk psd_mom_min psd_bins_per_dec_mom
        bin = trunc(Int, log10(pt_sk/psd_mom_min) * psd_bins_per_dec_mom )  +  1
    end

    # Sanity check
    if bin > num_psd_mom_bins
        @warn("Particle momentum exceeded PSD's bounds!",
              bin, num_psd_mom_bins, pt_sk/(mp_cgs*c_cgs),
              psd_mom_min * exp10(num_psd_mom_bins) / (mp_cgs*c_cgs))
        bin = num_psd_mom_bins
    end

    return bin
end

"""
Given parallel component and total of a particle's momentum in the shock frame,
determine which bin of psd particle will fall into.

Binning is done by the value of pt_sk in code units.
Angular binning is done with cos(θ) for large angles or θ for small angles.

The value of θ_fine marks the division between linear bin spacing (above θ_fine)
and logarithmic spacing (from θ_fine down to θ_min). In both logarithmic regions
the parameter bins_per_decade_*** determines the fineness of the bins.


!!! warning
    binning in angle is actually done with the *NEGATIVE* of the particle's
    cosine. This lets the most finely spaced bins correspond to UpS-pointing
    particles rather than DwS-pointing ones. This also has the effect that
    angles are essentially measured from the -x axis rather than the +x axis.

### Arguments
- px_sk: component of momentum parallel to shock normal (to B-field?), in code units
- pt_sk: total particle momentum, in code units

### Returns
Bin in angle into which particle falls
"""
function get_psd_bin_angle(px_sk, pt_sk, psd_bins_per_dec_θ, num_psd_θ_bins, psd_cos_fine, Δcos, psd_θ_min)
    # Bin in angle (bin); note that we negate the pitch angle to provide the
    # finest resolution (i.e. the logarithimcally-spaced angle bins rather than
    # the linearly-spaced cosine bins) for particles that are directed upstream
    p_cos = -px_sk / pt_sk

    if p_cos < psd_cos_fine
        # Pitch angle falls within linear spacing
        bin = num_psd_θ_bins - trunc(Int, (p_cos + 1) / Δcos )
    else
        θ = acos(p_cos) # Particle falls into logarithmic spacing region
        bin = θ < psd_θ_min ? 0 : # θ is close enough to zero that it might as well be
            trunc(Int, log10(θ/psd_θ_min) * psd_bins_per_dec_θ)  +  1
    end

    # Check to avert floating point error
    bin = min(bin, num_psd_θ_bins)

    return bin
end
