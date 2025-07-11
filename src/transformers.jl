module transformers
using StaticArrays: SVector
using LinearAlgebra: norm
using Unitful, UnitfulAstro, UnitfulGaussian, UnitfulEquivalences
using Unitful: mp, c    # physical constants

using ..parameters: psd_max
using ..constants: E₀ₚ

export transform_p_PS, transform_p_PSP, transform_psd_corners, get_transform_dN

"""
    get_transform_dN(...)

Calculate dN(p) (a 1-D array) for the passed slice of PSD in the specified inertial frame.

### Arguments

- `psd`: the (2-D) slice of the larger phase space distribution
- `m`: integer specifying frame into which we're transforming
- `transform_corner_**`: array holding transformed corner values, both ptot and cos(θ)
- `γᵤ`: conversion factor from number density to flux
- `i_approx`: degree of approximation to use in computing the `dNdp_out`

### Returns
`dN_out`: dN(p) for the given slice of PSD once transformed into the specified frame
"""
function get_transform_dN(
        psd, m, transform_corner_pt, transform_corner_ct, γᵤ, i_approx::Integer,
        num_psd_mom_bins, psd_mom_bounds)

    dN_out = zeros(psd_mom_axis)

    # Loop over cos(θ) and ptot space to re-bin input PSD slice
    #--------------------------------------------------------------------------
    for j in 0:psd_max, i in 0:psd_max

        psd[i,j] < 1e-66 && continue # Skip empty cells in PSD

        # In the below block, "cell_weight" includes n₀⋅u₁ normalization and division by
        # |vₓ| to change #/(cm²⋅s) to #/cm³. Divide by γ of flow speed as part of
        # Lorentz transformation of phase-space density (Rybicki & Lightman, p.146)
        #
        # Obtain cell_weight, p_cell_lo and p_cell_hi
        #----------------------------------------------------------------------
        cell_weight = psd[i, j] / γᵤ
        (pt_lo_pt, pt_lo_ct, pt_hi_pt, pt_hi_ct,
         ct_lo_pt, ct_lo_ct, ct_hi_pt, ct_hi_ct,
         pt_lo_tied, pt_hi_tied) = identify_corners(i, j, transform_corner_pt, transform_corner_ct, m)
        # Use of m in argument list signals what frame we're in; not used except in case of error/printout

        p_cell_lo = pt_lo_pt
        p_cell_hi = pt_hi_pt
        #----------------------------------------------------------------------
        # cell_weight, p_cell_lo, p_cell_hi found


        # Find lower and upper boundaries of psd_mom_bounds that we'll be
        # dealing with based on p_cell_lo and p_cell_hi
        #----------------------------------------------------------------------
        l_lo = findfirst(>(p_cell_lo), psd_mom_bounds) - 1

        # Error check to make sure that Lorentz transformations give reasonable results for l_lo
        if p_cell_lo > psd_mom_bounds[end]
            @info(" In get_dNdp_cr, p_cell_lo > psd_mom_max!",
                  m, j, i, p_cell_lo, psd_mom_bounds[end])
            l_lo = num_psd_mom_bins
        end


        l_hi = findnext(≥(p_cell_hi), # find phase space break FIXME write better comments
                        psd_mom_bounds,
                        l_lo) # start searching from l_lo


        # Error check to make sure that Lorentz transformations give reasonable results
        # for l_hi
        if p_cell_hi > psd_mom_bounds[num_psd_mom_bins+1]
            @info(" In get_dNdp_cr, p_cell_hi > psd_mom_max!",
                  m, j, i, p_cell_hi, psd_mom_bounds[num_psd_mom_bins+1])
            l_hi = num_psd_mom_bins
        end
        #----------------------------------------------------------------------
        # l_lo and l_hi found


        #----------------------------------------------------------------------
        # Distribute the weight value of cell_weight into the appropriate bins of dN_out,
        # assigning fractional weights where the cell of PSD covers more than one bin
        # of psd_mom_bounds
        # NOTE: virtually every line referring to psd_mom_bounds uses the top
        # edge of the bin under consideration (i.e. the value of the next bin
        # of psd_mom_bounds) because this code block originally used xph
        #----------------------------------------------------------------------

        # If assuming uniform distribution of cell_weight between p_cell_lo and p_cell_hi
        #----------------------------------------------------------------------
        if i_approx == 0
            uniform_cell_distribution!(dN_out, l_lo, l_hi, p_cell_lo, p_cell_hi)
        #----------------------------------------------------------------------
        # i_approx = 0 finished


        # If assuming isosceles trianglular distribution of cell_weight, peak is located
        # above geometric mean of p_cell_lo and p_cell_hi. If assuming scalene
        # triangular distribution of cell_weight, peak of triangle is located on (geometric)
        # mean of ct_hi_pt and ct_lo_pt, noting that both are logarithms.
        # Only difference is location of p_peak, so handle both with same block of code.
        #----------------------------------------------------------------------
        elseif i_approx == 1 || i_approx == 2
            triangular_distribution!(dN_out)
        #----------------------------------------------------------------------
        # i_approx = 1, i_approx = 2 finished


        # If not assuming anything about shape, i.e. performing exact
        # calculation of fractional areas
        #----------------------------------------------------------------------
        elseif i_approx == 3

            error("i_approx = 3 not currently enabled")

            # Determine which of ct_lo and ct_hi is peak_left and peak_right
            # Find heights of both peaks, i.e. vertical distance between, e.g.,
            # ct_lo and line segment connecting ct_hi and pt_**
            # Run through same process as in original version of subroutine,
            # dividing shape into zones for determining partial areas. Will be
            # easier this time, though, since bottom line is flat

        #----------------------------------------------------------------------
        # i_approx = 3 finished

        else # Invalid selection of i_approx. Flag error and stop program
            throw(DomainError(i_approx, "i_approx must be 0, 1, 2, or 3"))
        end
        #----------------------------------------------------------------------
        # cell_weight distributed


        # This block tracks the distribution of pitch angles in the shock and plasma frames.
        # In the interest of speed, it uses only the equal-weights method (i_approx = 0 above).
        # Because ct_bounds DECREASES as index increases (since it's derived from increasing θ),
        # many things are *slightly* different from the i_approx = 0 case for total momentum.
        # WARNING: this must be re-checked and re-tested since it was brought into new version of code
        #----------------------------------------------------------------------
        #track_pitch_angles() # TODO figure out arguments
        # pitch angles tracked

    end # loop over phase-space momentum space cells

    return dN_out
end

"""
    uniform_cell_distribution!(...)

TODO
"""
function uniform_cell_distribution!(dN_out, l_lo, l_hi, p_cell_lo, p_cell_hi) # TODO fix arguments
    length_tot = 1 / (p_cell_hi - p_cell_lo)
    p_bottom   = p_cell_lo

    for l in l_lo:l_hi
        # Cell fits entirely within bin of psd_mom_bounds
        if p_cell_hi < psd_mom_bounds[l_lo+1]
            dN_out[l] += cell_weight
            break
        end

        # Top of bin in psd_mom_bounds less than p_cell_hi
        if psd_mom_bounds[l+1] < p_cell_hi
            dN_out[l] += cell_weight*(psd_mom_bounds[l+1] - p_bottom) / length_tot
            # Adjust p_bottom to mark counting of current bin
            p_bottom = psd_mom_bounds[l+1]
        end

        # Top of bin in psd_mom_bounds is equal to/greater than p_cell_hi
        if psd_mom_bounds[l+1] ≥ p_cell_hi
            dN_out[l] += cell_weight*(p_cell_hi - psd_mom_bounds[l]) / length_tot
            break
        end
    end
    return dN_out
end

"""
    triangular_distribution!(...)

TODO
"""
function triangular_distribution!(dN_out) # TODO fix arguments
    length_tot = 1 / (p_cell_hi - p_cell_lo)
    ct_height  = 2 * cell_weight / length_tot # A = 1/2*b*h

    p_bottom     = p_cell_lo
    if i_approx == 1
        p_peak   = (p_cell_lo + p_cell_hi) / 2
    else # i_approx == 2
        p_peak   = (ct_lo_pt + ct_hi_pt) / 2
    end
    p_denom_lo   = 1 / (p_peak - p_cell_lo)
    p_denom_hi   = 1 / (p_cell_hi - p_peak)

    fractional_area = 0 # Total amount of cell_weight accounted for

    for l in l_lo:l_hi

        # Cell fits entirely within bin of psd_mom_bounds
        #------------------------------------------------------------------
        if p_cell_hi < psd_mom_bounds[l_lo+1]
            dN_out[l] += cell_weight
            break
        end
        #------------------------------------------------------------------
        # cell fits entirely within bin of psd_mom_bounds


        # Top of bin in psd_mom_bounds less than or equal to p_peak.
        # Area is
        # - a triangle if p_bottom = p_cell_lo, or
        # - a trapezoid if p_bottom > p_cell_lo
        #------------------------------------------------------------------
        if psd_mom_bounds[l+1] ≤ p_peak
            # Calculate current base
            p_base = psd_mom_bounds[l+1] - p_bottom

            # Calculate right-hand height
            ct_rh_height = (psd_mom_bounds[l+1] - p_cell_lo) * p_denom_lo * ct_height

            # Calculate left-hand height, taking advantage of right triangle similarity
            if p_bottom == p_cell_lo
                ct_lh_height = 0.0
            else
                ct_lh_height = (p_bottom - p_cell_lo) * p_denom_lo * ct_height
            end

            # Calculate partial area and add it to dN_out
            partial_area = p_base/2 * (ct_lh_height + ct_rh_height)
            dN_out[l] += partial_area

            # Adjust p_bottom to mark counting of current bin, update fractional_area...
            p_bottom = psd_mom_bounds[l+1]
            fractional_area += partial_area

            # ...and move on to next bin
            continue
        end
        #------------------------------------------------------------------
        # top of bin in psd_mom_bounds between p_cell_lo and p_peak


        # Top of bin in psd_mom_bounds between p_peak and p_cell_hi; area
        # is remaining area minus missing triangle at right
        #------------------------------------------------------------------
        if psd_mom_bounds[l+1] < p_cell_hi
            # Calculate missing base
            p_base = p_cell_hi - psd_mom_bounds[l+1]

            # Calculate left-hand height of missing area, taking advantage
            # of right triangle similarity
            ct_lh_height = p_base * p_denom_hi * ct_height

            # Calculate missing area and partial area and add it to dN_out
            missing_area = p_base/2 * ct_lh_height
            partial_area = (cell_weight - fractional_area) - missing_area
            dN_out[l] += partial_area

            # Adjust p_bottom to mark counting of current bin, update fractional_area
            p_bottom = psd_mom_bounds[l+1]
            fractional_area += partial_area

            # ...and move on to next bin
            continue
        end
        #------------------------------------------------------------------
        # top of bin in psd_mom_bounds between p_peak and p_cell_hi


        # Top of bin in psd_mom_bounds above p_cell_hi; area is remaining
        # fraction of total
        #------------------------------------------------------------------
        if psd_mom_bounds[l+1] ≥ p_cell_hi
            # Calculate partial area and add it to dN_out
            dN_out[l] += cell_weight - fractional_area
            break
        end
        #------------------------------------------------------------------
        # top of bin in psd_mom_bounds above p_cell_hi
    end
    return dN_out
end

"""
    track_pitch_angles(...)

TODO
"""
function track_pitch_angles() # TODO figure out arguments
    # Adjust cell_weight to remove the velocity-weighting applied in PSD
    if i == 0
        pt_sk = 0.0
        i_pt_sk  = i_ct_pt_sk_min
    else
        pt_sk = exp10(psd_mom_bounds[i] + psd_mom_bounds[i+1])
        pt_sk = √pt_sk * mp * utsk_cm
        i_pt_sk  = floor(Int, psd_mom_bounds[i])
    end
    γₚ_sk = hypot(1, pt_sk/(rest_mass*c))

    cell_weight = cell_weight * pt_sk / (γₚ_sk * rest_mass) / proton_num_density_upstream

    # Binning shock frame values very easy; just add directly to correct bin of histogram
    if m == 1
        cθ_sk_xw[j,i_pt_sk] += cell_weight
    end


    # Identify the min and max extent of the plasma frame cell,
    # as well as the correct decade of momentum for binning
    ct_cell_lo = min(pt_lo_ct, pt_hi_ct, ct_lo_ct, ct_hi_ct)
    ct_cell_hi = max(pt_lo_ct, pt_hi_ct, ct_lo_ct, ct_hi_ct)
    i_pt_pf = floor(Int, (pt_lo_pt + pt_hi_pt + ct_lo_pt + ct_hi_pt)/4)


    # Determine the spread in cos(θ)
    # Remember that ct_bounds counts DOWN from +1 to -1 !!!!
    l_hi = findfirst(<(ct_cell_hi), ct_bounds) - 1
    l_lo = findnext(≤(ct_cell_lo), ct_bounds, l_hi)


    # Distribute cell weight among all bins crossed by cell
    ct_length_tot = ct_cell_hi - ct_cell_lo
    ct_bottom     = ct_cell_lo

    for l in l_lo:-1:l_hi
        # Cell fits entirely within bin of ct_bounds
        if ct_cell_hi < ct_bounds[l_lo-1]
            if m == 2
                cθ_pf_xw[l-1, i_pt_pf] += cell_weight
            elseif m == 3
                cθ_ef_xw[l-1, i_pt_pf] += cell_weight
            end
            break
        end

        # Top of bin in ct_bounds less than ct_cell_hi
        if ct_bounds[l-1] < ct_cell_hi
            frac_ct_length = (ct_bounds[l-1] - ct_bottom) / ct_length_tot

            if m == 2
                cθ_pf_xw[l-1, i_pt_pf] += cell_weight*frac_ct_length
            elseif m == 3
                cθ_ef_xw[l-1, i_pt_pf] += cell_weight*frac_ct_length
            end

            # Adjust ct_bottom to mark counting of current bin
            ct_bottom = ct_bounds[l-1]
        end

        # Top of bin in ct_bounds is ≥ ct_cell_hi
        if ct_bounds[l-1] ≥ ct_cell_hi
            frac_ct_length = (ct_cell_hi - ct_bounds[l]) / ct_length_tot

            if m == 2
                cθ_pf_xw[l-1, i_pt_pf] += cell_weight*frac_ct_length
            elseif m == 3
                cθ_ef_xw[l-1, i_pt_pf] += cell_weight*frac_ct_length
            end

            break
        end
    end # loop over bins of ct_bounds
end

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

"""
    transform_p_PS(...)

Takes particle's Plasma frame momentum components and transforms them to the shock frame.
Transformation formulas are correct for all obliquities.

Notes from Glen Double (31 Oct 2002):
Method: define flow velocity and particle momentum in component vectors with respect to
the plasma xyz frame. Find p_para (with respect to u) component of particle momentum
using scalar product of p and u. Then `p_perp = p - p_para` (using vectors).
Next, make relativistic transformation on p_para while `p_perp` remains constant.
Finally, recreate `p_rel = p_rel_para + p_perp`, and find `p_rel` components in the
original xyz frame by taking scalar products along the xyz axes.

### Arguments

- `aa`: particle atomic mass
- `pb_pf`: component of ptot_pf parallel to magnetic field
- `p_perp_b_pf`: component of ptot_pf perpendicular to magnetic field
- `γₚ_pf`: Lorentz factor associated with ptot_pf
- `φ_rad`: phase angle of gyration; looking upstream, counts clockwise from +z axis
- `uₓ_sk`: bulk flow speed along x axis
- `uz_sk`: bulk flow speed along z axis
- `utot`: total bulk flow speed
- `γᵤ_sf`: Lorentz factor associated with `utot`
- `b_cosθ`: component of magnetic field along x axis
- `b_sinθ`: component of magnetic field along z axis

### Returns

- `ptot_sk`: total shock frame momentum in new grid zone
- `p_sk`: components of shock frame momentum
- `γₚ_sk`: Lorentz factor associated with ptot_sk
"""
function transform_p_PS(
        aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad, uₓ_sk, uz_sk, utot,
        γᵤ_sf, b_cosθ, b_sinθ,
        mc)

    φ_p = φ_rad + π/2

    p_p_cos = p_perp_b_pf * cos(φ_p)


    # xyz plasma frame components
    #CHECKTHIS: θ in shock frame may be different from θ in plasma frame because of lorentz
    # transformation between the two. Are the next lines correct in light of this?
    #@debug "" pb_pf p_perp_b_pf p_p_cos φ_p b_cosθ b_sinθ
    p_pf = SVector(pb_pf*b_cosθ - p_p_cos*b_sinθ,
                   p_perp_b_pf * sin(φ_p),
                   pb_pf*b_sinθ + p_p_cos*b_cosθ)

    # xyz shock frame components
    p_sk = SVector(γᵤ_sf * (p_pf.x + γₚ_pf * aa*mp * uₓ_sk),
                   p_pf.y,
                   p_pf.z)

    # Parallel/perpendicular (new) shock frame components
    ptot_sk = norm(p_sk)
    pb_sk = p_sk.x*b_cosθ + p_sk.z*b_sinθ
    #@debug "" ptot_sk pb_sk

    if ptot_sk < abs(pb_sk)
        p_perp_b_sk = 1e-6 * ptot_sk
        pb_sk       = copysign(√(ptot_sk^2 - p_perp_b_sk^2), pb_sk)
        #CHECKTHIS: does this *ever* happen?!
        @warn("ptot_sk < pb_sk in transform_p_PS")
    else
        p_perp_b_sk = √(ptot_sk^2 - pb_sk^2)
    end

    γₚ_sk = hypot(ptot_sk/mc, 1)

    return ptot_sk, p_sk, γₚ_sk
end

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

"""
    transform_p_PSP(...)

Takes particle's old plasma frame momentum components and transforms them twice:
from old Plasma frame to new Shock frame, then from new Shock frame to new Plasma frame.
Only ever called if particle scattered between zones with different bulk flow velocities.
Transformation formulas are correct for all obliquities.

Notes from Glen Double (31 Oct 2002):
Method: define flow velocity and particle momentum in component vectors with respect to the
plasma xyz frame. Find `p_para` (with respect to u) component of particle momentum using
scalar product of p and u. Then `p_perp = p - p_para` (using vectors). Next, make relativistic
transformation on `p_para` while `p_perp` remains constant. Finally, recreate
`p_rel = p_rel_para + p_perp`, and find `p_rel` components in the original xyz frame by taking
scalar products along the xyz axes.

### Arguments

- `aa`: particle atomic mass
- `uₓ_sk`/`old`: current and old bulk flow speed along x axis
- `uz_sk`/`old`: current and old bulk flow speed along z axis
- `utot`/`old`: current and old total bulk flow speed
- `γᵤ_sf`/`old`: Lorentz factor associated with utot/old
- `b_cosθ`/`old`: current and old component of magnetic field along x axis
- `b_sinθ`/`old`: current and old component of magnetic field along z axis

### Returns

- `ptot_pf`: total plasma frame momentum in new grid zone
- `ptot_sk`: total shock frame momentum in new grid zone
- `p_sk`: xyz components of ptot_sk
- `pb_sk`: component of ptot_sk parallel to magnetic field
- `p_perp_b_sk`: component of ptot_sk perpendicular to magnetic field
- `γₚ_sk`: Lorentz factor associated with ptot_sk

### Modifies

- `pb_pf`: component of ptot_pf parallel to magnetic field
- `p_perp_b_pf`: component of ptot_pf perpendicular to magnetic field
- `γₚ_pf`: Lorentz factor associated with ptot_pf
- `φ_rad`: phase angle of gyration; looking upstream, counts clockwise from +z axis
"""
function transform_p_PSP(
        aa, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad,
        uₓ_sk_old, uz_sk_old, utot_old, γᵤ_sf_old, b_cos_old, b_sin_old,
        uₓ_sk, uz_sk, utot, γᵤ_sf, b_cosθ, b_sinθ,
        mc)

    φ_p = φ_rad + π/2

    p_p_cos = p_perp_b_pf * cos(φ_p)

    # xyz (old) plasma frame components
    p_pf = SVector(pb_pf*b_cos_old - p_p_cos*b_sin_old,
                   p_perp_b_pf * sin(φ_p),
                   pb_pf*b_sin_old + p_p_cos*b_cos_old)

    # xyz (new) shock frame components
    p_sk = SVector(
        (((γᵤ_sf_old-1) * (uₓ_sk_old/utot_old)^2 + 1) * p_pf.x +
         (γᵤ_sf_old-1) * (uₓ_sk_old*uz_sk_old/utot_old^2) * p_pf.z +
         γᵤ_sf_old * γₚ_pf * aa*mp * uₓ_sk_old),
        p_pf.y,
        ((γᵤ_sf_old-1) * (uₓ_sk_old*uz_sk_old/utot_old^2) * p_pf.x +
         ((γᵤ_sf_old-1) * (uz_sk_old/utot_old)^2 + 1) * p_pf.z +
         γᵤ_sf_old * γₚ_pf * aa*mp * uz_sk_old)
    )

    # Parallel/perpendicular (new) shock frame components
    ptot_sk = norm(p_sk)
    pb_sk = p_sk.x*b_cosθ + p_sk.z*b_sinθ

    if ptot_sk < abs(pb_sk)
        p_perp_b_sk = 1e-6 * ptot_sk
        pb_sk = copysign(√(ptot_sk^2 - p_perp_b_sk^2), pb_sk)
        @warn("ptot_sk < pb_sk in transform_p_PSP")
    else
        p_perp_b_sk = √(ptot_sk^2 - pb_sk^2)
    end

    γₚ_sk = hypot(ptot_sk/mc, 1)


    # xyz (new) plasma frame components
    p_pf = SVector(( # x-component
                    ((γᵤ_sf - 1) * (uₓ_sk/utot)^2 + 1) * p_sk.x + (γᵤ_sf - 1) *
                    (uₓ_sk*uz_sk/utot^2) * p_sk.z - γᵤ_sf * γₚ_sk * aa*mp * uₓ_sk
                   ),
                   p_sk.y,
                   ( # z-component
                    (γᵤ_sf - 1) * (uₓ_sk*uz_sk/utot^2) * p_sk.x +
                    ((γᵤ_sf - 1) * (uz_sk/utot)^2 + 1) * p_sk.z - γᵤ_sf * γₚ_sk * aa*mp * uz_sk
                   ))


    # Parallel/perpendicular (new) plasma frame components, including new phase angle
    ptot_pf = norm(p_pf)
    pb_pf = p_pf.x*b_cosθ + p_pf.z*b_sinθ
    if ptot_pf < abs(pb_pf)
        p_perp_b_pf = 1e-6 * ptot_pf
        pb_pf       = copysign(√(ptot_pf^2 - p_perp_b_pf^2), pb_pf)
        @warn("ptot_pf < pb_pf in transform_p_PSP")
    else
        p_perp_b_pf = √(ptot_pf^2 - pb_pf^2)
    end

    γₚ_pf = hypot(ptot_pf/mc, 1)

    # See Figure 14 of Ellison, Baring, Jones (1996) [1996ApJ...473.1029E] for more details on φ_p
    φ_p = atan(p_pf.y, -p_pf.x*b_sinθ + p_pf.z*b_cosθ)
    φ_rad = φ_p - π/2

    return ptot_pf, ptot_sk, p_sk, pb_sk, p_perp_b_sk, γₚ_sk, pb_pf, p_perp_b_pf, γₚ_pf, φ_rad
end

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

"""
    transform_psd_corners(...)

Given a relative Lorentz factor between two frames, transform the corners of the PSD into
the new frame. Outputs total momentum as log of cgs values.

TODO: remove conversion to/from log space here and in `get_dndp_cr`

### Arguments
FIXME with actual argument list
- `γ_in`: relative Lorentz factor between shock and resultant frame

### Returns

- `transform_corner_pt`: total momenta at the corners
- `transform_corner_ct`: cos(θ) (NOT θ!!!) values at the corners
"""
function transform_psd_corners(
        γ_in,
        aa_ion, psd_lin_cos_bins, num_psd_θ_bins, psd_θ_bounds,
        num_psd_mom_bins, psd_mom_bounds, i_ion)

    # Administrative constants
    rest_mass_energy = aa_ion[i_ion] * E₀ₚ
    βᵤ = γ_in ≥ 1.000001 ? √(1 - 1/γ_in^2) : 0.0 # Prevent floating point issues

    # Fill transform_corner_** arrays, looping over angle outermost
    transform_corner_pt = OffsetMatrix{Float64}(undef, 0:num_psd_mom_bins+1, 0:num_psd_θ_bins+1)
    transform_corner_ct = OffsetMatrix{Float64}(undef, 0:num_psd_mom_bins+1, 0:num_psd_θ_bins+1)
    for j in eachindex(psd_θ_bounds)

        # Determine current cosine, remembering that psd_θ_bounds has both a linearly-spaced
        # region in cosine and logarithmically-spaced region in θ. Also need to remember that
        # the most finely spaced bins should occur in the upstream-pointing direction, so need to
        # negate psd_θ_bounds to get true cosine value.
        if j > (num_psd_θ_bins - psd_lin_cos_bins)
            cosθ = -psd_θ_bounds[j]
        else
            cosθ = -cos(psd_θ_bounds[j])
        end

        # With angle fixed, loop over total momenta
        for i in eachindex(psd_mom_bounds)

            # psd_mom_bounds uses logarithmic spacing for its bins, so undo that before
            # continuing the calculation
            pt_sk = exp10(psd_mom_bounds[i])

            #if i == 0
            #    pt_sk = 0.0 # Edge case when i = 0
            #end

            pₓ_sk   = pt_sk * cosθ
            etot_sk = hypot(pt_sk*c, rest_mass_energy)

            pₓ_Xf  = γ_in * (pₓ_sk - βᵤ*etot_sk/c)
            pt_Xf  = √(pt_sk^2 + pₓ_Xf^2 - pₓ_sk^2)

            # Transform to log space because get_dNdp_cr expects it
            transform_corner_pt[i,j] = log10(pt_Xf)
            transform_corner_ct[i,j] = pₓ_Xf / pt_Xf


        end # loop over momentum
    end # loop over θ

    return transform_corner_pt, transform_corner_ct
end
end # module transformers
