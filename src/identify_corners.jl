using .parameters: psd_max

"""
Given particular cell in PSD, identify the corners with lowest and highest
total momenta and the lower/higher of the two remaining cos(θ) values.

### Arguments

- i: marker for location on the total momentum axis
- j: marker for location on the cos(θ) axis
- transform_corn_pt: values of total momentum for the corners of PSD, transformed into desired frame
- transform_corn_ct: values of cos(θ) for the corners of PSD, transformed into desired frame
- ii_sk_pf: integer giving frame into which corners were transformed

### Returns

- pt_lo_pt: momentum of the corner with lowest total momentum
- pt_lo_ct: cos(θ) of the corner with lowest total momentum
- pt_hi_pt: momentum of the corner with highest total momentum
- pt_hi_ct: cos(θ) of the corner with highest total momentum
- ct_lo_pt: momentum of the remaining corner with lower cos(θ)
- ct_lo_ct: cos(θ) of the remaining corner with lower cos(θ)
- ct_hi_pt: momentum of the remaining corner with higher cos(θ)
- ct_hi_ct: cos(θ) of the remaining corner with higher cos(θ)
- pt_lo_tied: flag used if multiple cells are tied for lowest ptot
- pt_hi_tied: flag used if multiple cells are tied for highest ptot
"""
function identify_corners(i, j, transform_corn_pt, transform_corn_ct, ii_sk_pf)

    # Get momentum and cos(θ) values from the grid
    corn_pts = [transform_corn_pt[i,   j],
                transform_corn_pt[i+1, j],
                transform_corn_pt[i,   j+1],
                transform_corn_pt[i+1, j+1]]
    corn_cts = [transform_corn_ct[i,   j],
                transform_corn_ct[i+1, j],
                transform_corn_ct[i,   j+1],
                transform_corn_ct[i+1, j+1]]

    # Pre-set logical array
    mask = fill(true, 4)

    # Pre-set flags for ties
    pt_lo_tied = 0
    pt_hi_tied = 0

    # Identify *a* corner with lowest total momentum;
    # if there are no ties, this corner is pt_lo
    #----------------------------------------------------------------------------
    pt_lo_pt, i_pt_lo = findmin(corn_pts)
    pt_lo_ct = corn_cts[i_pt_lo]
    mask[i_pt_lo] = false

    # Count number of corners with identified lowest momentum; if this number
    # isn't 1, flag as needing additional attention
    if count(corn_pts .== pt_lo_pt) > 1
        # There was a tie, so set flag for special handling later
        @error("Multiple corners with lowest momentum", ii_sk_pf, i, j, corn_pts)
        pt_lo_tied = 1
        #exit()
    elseif count(corn_pts .== pt_lo_pt) < 1
        # Who knows what happened
        error("No corner with lowest momentum", ii_sk_pf, i, j, pt_lo_pt, corn_pts)
    end
    #-------------------------------------------------------------------------
    # pt_lo identified


    # Identify *a* corner with highest total momentum; if there are no ties, this corner is pt_hi
    #----------------------------------------------------------------------------
    pt_hi_pt, i_pt_hi = findmax(corn_pts)
    pt_hi_ct = corn_cts[i_pt_hi]
    mask[i_pt_hi] = false

    # Count number of corners with identified highest momentum;
    # if this number isn't 1, flag as needing additional attention
    if count(corn_pts .== pt_hi_pt) > 1
        # There was a tie, so set flag for special handling later
        @error("multiple corners with highest momentum", ii_sk_pf, i, j, corn_pts)
        pt_hi_tied = 1
        #exit()
    elseif count(corn_pts .== pt_hi_pt) < 1
        # Who knows what happened
        error("No corner with highest momentum", ii_sk_pf, i, j, pt_hi_pt, corn_pts)
    end
    #-------------------------------------------------------------------------
    # pt_hi identified


    # Between two remaining corners, identify higher/lower values of cos(θ).
    # If the two values aren't identical, one corner will be ct_hi and other will be ct_lo.
    # NOTE: in event that ct_hi_ct and ct_lo_ct are equal, assign corner with lower total momentum as ct_lo
    #----------------------------------------------------------------------------
    j_ct_hi  = maxloc(corn_cts, 1, mask) # FIXME
    ct_hi_pt = corn_pts[j_ct_hi]
    ct_hi_ct = corn_cts[j_ct_hi]
    mask[j_ct_hi] = false

    j_ct_lo  = minloc(corn_cts, 1, mask) # FIXME
    ct_lo_pt = corn_pts[j_ct_lo]
    ct_lo_ct = corn_cts[j_ct_lo]

    # Make sure the two corners aren't somehow identical
    if ct_hi_ct == ct_lo_ct

        # ct_hi and ct_lo have same value of cos(θ); use total momentum as the tiebreaker
        if ct_hi_pt > ct_lo_pt
            # ct_hi remains ct_hi by virtue of its higher momentum

        elseif ct_hi_pt < ct_lo_pt
            # Switch ct_hi and ct_lo because of momentum ordering
            ct_hi_pt = corn_pts[j_ct_lo]
            ct_hi_ct = corn_cts[j_ct_lo]
            ct_lo_pt = corn_pts[j_ct_hi]
            ct_lo_ct = corn_cts[j_ct_hi]

        else
            # ct_hi and ct_lo have exactly the same values in momentum and in cos(θ)
            error("ct_hi and ct_lo identical", i, j, ct_hi_pt, ct_hi_ct, ct_lo_pt, ct_lo_ct)
        end
    end
    #-------------------------------------------------------------------------
    # ct_hi and ct_lo identified


    # In the case of a tie for lowest/highest momentum, determine which of ct_lo and ct_hi is
    # the source of the tie. Also make sure the pt_** involved in the tie has the lower cos(θ).
    #-------------------------------------------------------------------------
    if pt_lo_tied == 1

        # Determine which of ct_lo and ct_hi is tied with pt_lo
        if pt_lo_pt == ct_lo_pt
            pt_lo_tied = 2
        elseif pt_lo_pt == ct_hi_pt
            pt_lo_tied = 3
        else
            error("Something has gone horribly wrong in identify_corners. Code 1")
        end

        # Of the two tied corners, pt_lo should have the lower cos(θ)
        if pt_lo_tied == 2

            # Corner tied with pt_lo is ct_lo
            if pt_lo_ct > ct_lo_ct
                # Switch pt_lo and ct_lo
                pt_lo_pt = corn_pts[j_ct_lo]
                pt_lo_ct = corn_cts[j_ct_lo]
                ct_lo_pt = corn_pts[i_pt_lo]
                ct_lo_ct = corn_cts[i_pt_lo]
            elseif pt_lo_ct < ct_lo_ct
                # Nothing to do here. Assignment is correct
            else
                error("pt_lo and ct_lo identical\n",
                      i, j, pt_lo_pt, pt_lo_ct, ct_lo_pt, ct_lo_ct,
                      "\nIf $i = 0, reduce EMNFC in mc_in.dat")
            end

        elseif pt_lo_tied == 3

            # Corner tied with pt_lo is ct_hi
            if pt_lo_ct > ct_hi_ct
                # Switch pt_lo and ct_hi
                pt_lo_pt = corn_pts[j_ct_hi]
                pt_lo_ct = corn_cts[j_ct_hi]
                ct_hi_pt = corn_pts[i_pt_lo]
                ct_hi_ct = corn_cts[i_pt_lo]
            elseif pt_lo_ct < ct_hi_ct
                # Nothing to do here. Assignment is correct
            else
                error("pt_lo and ct_hi identical", i, j, pt_lo_pt, pt_lo_ct, ct_hi_pt, ct_hi_ct,
                      "If $i = 0, reduce EMNFC in mc_in.dat")
            end

        end

    end

    if pt_hi_tied == 1

        # Determine which of ct_lo and ct_hi is tied with pt_hi
        if pt_hi_pt == ct_lo_pt
            pt_hi_tied = 2
        elseif pt_hi_pt == ct_hi_pt
            pt_hi_tied = 3
        else
            error("Something has gone horribly wrong in identify_corners. Code 2")
        end

        # Of the two tied corners, pt_lo should have the lower cos(θ)
        if pt_hi_tied == 2

            # Corner tied with pt_hi is ct_lo
            if pt_hi_ct > ct_lo_ct
                # Switch pt_hi and ct_lo
                pt_hi_pt = corn_pts[j_ct_lo]
                pt_hi_ct = corn_cts[j_ct_lo]
                ct_lo_pt = corn_pts[i_pt_hi]
                ct_lo_ct = corn_cts[i_pt_hi]
            elseif pt_hi_ct < ct_lo_ct
                # Nothing to do here. Assignment is correct
            else
                error("pt_hi and ct_lo identical",
                      i, j, pt_hi_pt, pt_hi_ct, ct_lo_pt, ct_lo_ct,
                      "\nIf $i = 0, reduce EMNFC in mc_in.dat")
            end

        elseif pt_hi_tied == 3

            # Corner tied with pt_hi is ct_hi
            if pt_hi_ct > ct_hi_ct
                # Switch pt_hi and ct_hi
                pt_hi_pt = corn_pts[j_ct_hi]
                pt_hi_ct = corn_cts[j_ct_hi]
                ct_hi_pt = corn_pts[i_pt_hi]
                ct_hi_ct = corn_cts[i_pt_hi]
            elseif pt_hi_ct < ct_hi_ct
                # Nothing to do here. Assignment is correct
            else
                error("pt_hi and ct_lo identical",
                      i, j, pt_hi_pt, pt_hi_ct, ct_hi_pt, ct_hi_ct,
                      "If $i = 0, reduce EMNFC in mc_in.dat")
            end

        end # check on pt_hi having lower momentum

    end
    #-------------------------------------------------------------------------
    # Ties handled

    return (pt_lo_pt, pt_lo_ct, pt_hi_pt, pt_hi_ct,
            ct_lo_pt, ct_lo_ct, ct_hi_pt, ct_hi_ct,
            pt_lo_tied, pt_hi_tied)
end
