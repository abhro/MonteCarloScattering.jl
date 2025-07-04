using .parameters: psd_max

"""
    identify_corners(...)

Given particular cell in PSD, identify the corners with lowest and highest
total momenta and the lower/higher of the two remaining cos(θ) values.

### Arguments

- `i`: marker for location on the total momentum axis
- `j`: marker for location on the cos(θ) axis
- `transform_corner_pt`: values of total momentum for the corners of PSD, transformed into desired frame
- `transform_corner_cθ`: values of cos(θ) for the corners of PSD, transformed into desired frame
- `ii_sk_pf`: integer giving frame into which corners were transformed

### Returns

- `pt_lo_pt`: momentum of the corner with lowest total momentum
- `pt_lo_cθ`: cos(θ) of the corner with lowest total momentum
- `pt_hi_pt`: momentum of the corner with highest total momentum
- `pt_hi_cθ`: cos(θ) of the corner with highest total momentum
- `cθ_lo_pt`: momentum of the remaining corner with lower cos(θ)
- `cθ_lo_cθ`: cos(θ) of the remaining corner with lower cos(θ)
- `cθ_hi_pt`: momentum of the remaining corner with higher cos(θ)
- `cθ_hi_cθ`: cos(θ) of the remaining corner with higher cos(θ)
- `pt_lo_tied`: flag used if multiple cells are tied for lowest ptot
- `pt_hi_tied`: flag used if multiple cells are tied for highest ptot
"""
function identify_corners(i, j, transform_corner_pt, transform_corner_cθ, ii_sk_pf)

    # Get momentum and cos(θ) values from the grid
    corner_pts = [transform_corner_pt[i,   j],   transform_corner_pt[i+1, j],
                  transform_corner_pt[i,   j+1], transform_corner_pt[i+1, j+1]]
    corner_cθs = [transform_corner_cθ[i,   j],   transform_corner_cθ[i+1, j],
                  transform_corner_cθ[i,   j+1], transform_corner_cθ[i+1, j+1]]

    # Pre-set logical array
    mask = fill(true, 4)

    # Pre-set flags for ties
    pt_lo_tied = 0
    pt_hi_tied = 0

    # Identify *a* corner with lowest total momentum;
    # if there are no ties, this corner is pt_lo
    #----------------------------------------------------------------------------
    pt_lo_pt, i_pt_lo = findmin(corner_pts)
    pt_lo_cθ = corner_cθs[i_pt_lo]
    mask[i_pt_lo] = false

    # Count number of corners with identified lowest momentum; if this number
    # isn't 1, flag as needing additional attention
    if count(corner_pts .== pt_lo_pt) > 1
        # There was a tie, so set flag for special handling later
        @error("Multiple corners with lowest momentum", ii_sk_pf, i, j, corner_pts)
        pt_lo_tied = 1
        #exit()
    elseif count(corner_pts .== pt_lo_pt) < 1
        # Who knows what happened
        error("No corner with lowest momentum", ii_sk_pf, i, j, pt_lo_pt, corner_pts)
    end
    #-------------------------------------------------------------------------
    # pt_lo identified


    # Identify *a* corner with highest total momentum; if there are no ties, this corner is pt_hi
    #----------------------------------------------------------------------------
    pt_hi_pt, i_pt_hi = findmax(corner_pts)
    pt_hi_cθ = corner_cθs[i_pt_hi]
    mask[i_pt_hi] = false

    # Count number of corners with identified highest momentum;
    # if this number isn't 1, flag as needing additional attention
    if count(corner_pts .== pt_hi_pt) > 1
        # There was a tie, so set flag for special handling later
        @error("multiple corners with highest momentum", ii_sk_pf, i, j, corner_pts)
        pt_hi_tied = 1
        #exit()
    elseif count(corner_pts .== pt_hi_pt) < 1
        # Who knows what happened
        error("No corner with highest momentum", ii_sk_pf, i, j, pt_hi_pt, corner_pts)
    end
    #-------------------------------------------------------------------------
    # pt_hi identified


    # Between two remaining corners, identify higher/lower values of cos(θ).
    # If the two values aren't identical, one corner will be cθ_hi and other will be cθ_lo.
    # NOTE: in event that cθ_hi_cθ and cθ_lo_cθ are equal, assign corner with lower total momentum as cθ_lo
    #----------------------------------------------------------------------------
    j_cθ_hi  = maxloc(corner_cθs, 1, mask) # FIXME
    cθ_hi_pt = corner_pts[j_cθ_hi]
    cθ_hi_cθ = corner_cθs[j_cθ_hi]
    mask[j_cθ_hi] = false

    j_cθ_lo  = minloc(corner_cθs, 1, mask) # FIXME
    cθ_lo_pt = corner_pts[j_cθ_lo]
    cθ_lo_cθ = corner_cθs[j_cθ_lo]

    # Make sure the two corners aren't somehow identical
    if cθ_hi_cθ == cθ_lo_cθ

        # cθ_hi and cθ_lo have same value of cos(θ); use total momentum as the tiebreaker
        if cθ_hi_pt > cθ_lo_pt
            # cθ_hi remains cθ_hi by virtue of its higher momentum

        elseif cθ_hi_pt < cθ_lo_pt
            # Switch cθ_hi and cθ_lo because of momentum ordering
            cθ_hi_pt = corner_pts[j_cθ_lo]
            cθ_hi_cθ = corner_cθs[j_cθ_lo]
            cθ_lo_pt = corner_pts[j_cθ_hi]
            cθ_lo_cθ = corner_cθs[j_cθ_hi]

        else
            # cθ_hi and cθ_lo have exactly the same values in momentum and in cos(θ)
            error("cθ_hi and cθ_lo identical", i, j, cθ_hi_pt, cθ_hi_cθ, cθ_lo_pt, cθ_lo_cθ)
        end
    end
    #-------------------------------------------------------------------------
    # cθ_hi and cθ_lo identified


    # In the case of a tie for lowest/highest momentum, determine which of cθ_lo and cθ_hi is
    # the source of the tie. Also make sure the pt_** involved in the tie has the lower cos(θ).
    #-------------------------------------------------------------------------
    if pt_lo_tied == 1

        # Determine which of cθ_lo and cθ_hi is tied with pt_lo
        if pt_lo_pt == cθ_lo_pt
            pt_lo_tied = 2
        elseif pt_lo_pt == cθ_hi_pt
            pt_lo_tied = 3
        else
            error("Something has gone horribly wrong in identify_corners. Code 1")
        end

        # Of the two tied corners, pt_lo should have the lower cos(θ)
        if pt_lo_tied == 2

            # Corner tied with pt_lo is cθ_lo
            if pt_lo_cθ > cθ_lo_cθ
                # Switch pt_lo and cθ_lo
                pt_lo_pt = corner_pts[j_cθ_lo]
                pt_lo_cθ = corner_cθs[j_cθ_lo]
                cθ_lo_pt = corner_pts[i_pt_lo]
                cθ_lo_cθ = corner_cθs[i_pt_lo]
            elseif pt_lo_cθ < cθ_lo_cθ
                # Nothing to do here. Assignment is correct
            else
                error("pt_lo and cθ_lo identical\n",
                      i, j, pt_lo_pt, pt_lo_cθ, cθ_lo_pt, cθ_lo_cθ,
                      "\nIf $i = 0, reduce EMNFC in mc_in.dat")
            end

        elseif pt_lo_tied == 3

            # Corner tied with pt_lo is cθ_hi
            if pt_lo_cθ > cθ_hi_cθ
                # Switch pt_lo and cθ_hi
                pt_lo_pt = corner_pts[j_cθ_hi]
                pt_lo_cθ = corner_cθs[j_cθ_hi]
                cθ_hi_pt = corner_pts[i_pt_lo]
                cθ_hi_cθ = corner_cθs[i_pt_lo]
            elseif pt_lo_cθ < cθ_hi_cθ
                # Nothing to do here. Assignment is correct
            else
                error("pt_lo and cθ_hi identical", i, j, pt_lo_pt, pt_lo_cθ, cθ_hi_pt, cθ_hi_cθ,
                      "If $i = 0, reduce EMNFC in mc_in.dat")
            end

        end

    end

    if pt_hi_tied == 1

        # Determine which of cθ_lo and cθ_hi is tied with pt_hi
        if pt_hi_pt == cθ_lo_pt
            pt_hi_tied = 2
        elseif pt_hi_pt == cθ_hi_pt
            pt_hi_tied = 3
        else
            error("Something has gone horribly wrong in identify_corners. Code 2")
        end

        # Of the two tied corners, pt_lo should have the lower cos(θ)
        if pt_hi_tied == 2

            # Corner tied with pt_hi is cθ_lo
            if pt_hi_cθ > cθ_lo_cθ
                # Switch pt_hi and cθ_lo
                pt_hi_pt = corner_pts[j_cθ_lo]
                pt_hi_cθ = corner_cθs[j_cθ_lo]
                cθ_lo_pt = corner_pts[i_pt_hi]
                cθ_lo_cθ = corner_cθs[i_pt_hi]
            elseif pt_hi_cθ < cθ_lo_cθ
                # Nothing to do here. Assignment is correct
            else
                error("pt_hi and cθ_lo identical",
                      i, j, pt_hi_pt, pt_hi_cθ, cθ_lo_pt, cθ_lo_cθ,
                      "\nIf $i = 0, reduce EMNFC in mc_in.dat")
            end

        elseif pt_hi_tied == 3

            # Corner tied with pt_hi is cθ_hi
            if pt_hi_cθ > cθ_hi_cθ
                # Switch pt_hi and cθ_hi
                pt_hi_pt = corner_pts[j_cθ_hi]
                pt_hi_cθ = corner_cθs[j_cθ_hi]
                cθ_hi_pt = corner_pts[i_pt_hi]
                cθ_hi_cθ = corner_cθs[i_pt_hi]
            elseif pt_hi_cθ < cθ_hi_cθ
                # Nothing to do here. Assignment is correct
            else
                error("pt_hi and cθ_lo identical",
                      i, j, pt_hi_pt, pt_hi_cθ, cθ_hi_pt, cθ_hi_cθ,
                      "If $i = 0, reduce EMNFC in mc_in.dat")
            end

        end # check on pt_hi having lower momentum

    end
    #-------------------------------------------------------------------------
    # Ties handled

    return (pt_lo_pt, pt_lo_cθ, pt_hi_pt, pt_hi_cθ,
            cθ_lo_pt, cθ_lo_cθ, cθ_hi_pt, cθ_hi_cθ,
            pt_lo_tied, pt_hi_tied)
end
