using .parameters: na_particles

"""
Takes the **_sav arrays filled over the course of loop_pt and splits the
saved particles to form the population of the next pcut

### Arguments
TODO
- n_pts_target
- n_saved
- grid_sav
- tcut_sav
- l_save
- DwS_sav
- inj_sav
- xwt_sav
- pt_pf_sav
- pb_pf_sav
- x_PT_cm_sav
- xn_per_sav
- prp_x_cm_sav
- acctime_sec_sav
- φ_rad_sav

### Returns
TODO

### Modifies
- wt_running: weight factor of each particle remaining after this pcut
"""
function new_pcut(
        n_pts_target, n_saved, l_save, grid_sav, DwS_sav,
        inj_sav, xwt_sav, pt_pf_sav, pb_pf_sav, x_PT_cm_sav, xn_per_sav,
        prp_x_cm_sav, acctime_sec_sav, φ_rad_sav, tcut_sav,
        n_pts_use, wt_running)

    grid_new        = zeros(Int, na_particles)
    tcut_new        = zeros(Int, na_particles)
    DwS_new         = zeros(Bool, na_particles)
    inj_new         = zeros(Bool, na_particles)
    xwt_new         = zeros(na_particles)
    pt_pf_new        = zeros(na_particles)
    pb_pf_new        = zeros(na_particles)
    x_PT_cm_new     = zeros(na_particles)
    xn_per_new      = zeros(na_particles)
    prp_x_cm_new    = zeros(na_particles)
    acctime_sec_new = zeros(na_particles)
    φ_rad_new       = zeros(na_particles)

    # Determine multiplicity of splitting; perhaps none needed
    i_mult = max(n_pts_target ÷ n_saved, 1) # In case n_pts_target drops btwn pcuts

    # Calculate effect on particle weights and the weighting factor of each
    # remaining particle in the simulation
    # CHECKTHIS: is that last claim still true if old particles are imported
    # into the simulation?
    wt_running /= i_mult


    # Perform the splitting
    n_pts_new = 0

    for j in 1:n_pts_use

        l_save[j] || continue # Don't multiply particles that weren't kept, obviously

        for i in 1:i_mult

            n_pts_new += 1

            xwt_new[n_pts_new]         = xwt_sav[j] / i_mult
            pt_pf_new[n_pts_new]        = pt_pf_sav[j]
            pb_pf_new[n_pts_new]        = pb_pf_sav[j]
            x_PT_cm_new[n_pts_new]     = x_PT_cm_sav[j]
            grid_new[n_pts_new]        = grid_sav[j]
            DwS_new[n_pts_new]         = DwS_sav[j]
            inj_new[n_pts_new]         = inj_sav[j]
            xn_per_new[n_pts_new]      = xn_per_sav[j]
            prp_x_cm_new[n_pts_new]    = prp_x_cm_sav[j]
            acctime_sec_new[n_pts_new] = acctime_sec_sav[j]
            φ_rad_new[n_pts_new]       = φ_rad_sav[j]
            tcut_new[n_pts_new]        = tcut_sav[j]

        end  # loop over splits
    end  # loop over saved particles


    # Zero out the unused portions of the *_new arrays, just as a precaution
    xwt_new[n_pts_new+1:end]         .= 0.0
    pt_pf_new[n_pts_new+1:end]        .= 0.0
    pb_pf_new[n_pts_new+1:end]        .= 0.0
    x_PT_cm_new[n_pts_new+1:end]     .= 0.0
    grid_new[n_pts_new+1:end]        .= 0
    DwS_new[n_pts_new+1:end]         .= false
    inj_new[n_pts_new+1:end]         .= false
    xn_per_new[n_pts_new+1:end]      .= 0.0
    prp_x_cm_new[n_pts_new+1:end]    .= 0.0
    acctime_sec_new[n_pts_new+1:end] .= 0.0
    φ_rad_new[n_pts_new+1:end]       .= 0.0
    tcut_new[n_pts_new+1:end]        .= 0

    return (grid_new, tcut_new, DwS_new, inj_new, xwt_new, pt_pf_new, pb_pf_new,
            x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, φ_rad_new,
            n_pts_new, wt_running)
end
