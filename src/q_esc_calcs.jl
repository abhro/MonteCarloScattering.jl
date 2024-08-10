using .constants: kB_cgs, mp_cgs, me_cgs, c_cgs
using .parameters: β_rel_fl, na_ions

"""
Use the Rankine-Hugoniot relations to calculate the escaping momentum & energy flux.
"""
function q_esc_calcs(
        γ_adiab,
        r_comp, rRH, u_Z, β_Z, γ_Z, oblique, n_ions,
        aa_ion, zz_ion, tZ_ion, denZ_ion, sc_electron, tZ_electron, γ_2, β_2, u_2)


    # Quick test first. If r_comp = rRH, we expect no escaping flux.
    r_comp == rRH && return 0.0, 0.0


    #--------------------------------------------------------------------------
    #  Four possibilities for R-H relations: rel/nonrel and parallel/oblique.
    #  Determine which of the four to use. Cutoff for rel/nonrel is set in
    #  module 'controls'
    #--------------------------------------------------------------------------
    relativistic = (β_Z ≥ β_rel_fl)

    γβ_2  = γ_2 * β_2
    γ_fac = γ_adiab / (γ_adiab - 1)

    if !oblique
        if relativistic     # Possibility 1: Relativistic, parallel
            q_esc_calcs_relativistic_parallel()
        else                # Possibility 2: Nonrelativistic, parallel
            q_esc_calcs_nonrelativistic_parallel()
        end
    else
        if relativistic     # Possibility 3: Relativistic, oblique
            q_esc_calcs_relativistic_oblique()
        else                # Possibility 4: Nonrelativistic, oblique
            q_esc_calcs_nonrelativistic_oblique()
        end
    end

    return q_esc_cal_px, q_esc_cal_energy
end


# Solution comes from Ellison (1985) [1985JGR....90...29E] (Eqs 8-10).
# Note assumption of zero escaping momentum flux, which is good to
# within a couple percent for strong nonrelativistic shocks.
# #TODO: check how much of a difference this assumption makes
function q_esc_calcs_nonrelativistic_parallel()
    # Calculate thermal pressure of far upstream gas
    pressure_Z   = dot(denZ_ion, tZ_ion) * kB_cgs
    ρ_Z          = dot(denZ_ion, aa_ion) * mp_cgs
    mask = (aa_ion .≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])

    # If electrons were not a separate species, add them in here
    if !sc_electron
        pressure_Z +=  density_electron * kB_cgs*tZ_electron
        ρ_Z        +=  density_electron * me_cgs
    end

    # Calculate UpS incoming energy flux   #assumecold
    F_px_UpS_fl = ρ_Z * u_Z^2  +  pressure_Z
    F_energy_UpS_fl = ρ_Z*u_Z^3/2  +  5//2 * pressure_Z * u_Z

    # Calculate far DwS density (Eq 8) and pressure (Eq 9)
    ρ_2 = ρ_Z * γ_Z*β_Z / γβ_2
    pressure_2 = F_px_UpS_fl  -  ρ_2*u_2^2

    # Calculate escaping energy flux using Eq (10)
    Q_en = F_energy_UpS_fl  -  ρ_Z*u_Z*u_2^2/2  -  pressure_2 * u_2 * γ_fac

    # Finally, put in units of F_enZ
    q_esc_cal_energy = Q_en / F_energy_UpS_fl
    q_esc_cal_px = 0.0
end

#  Solution comes from Ellison+ (1990) [1991ApJ...378..214E]. Uses
#  relativistic Rankine-Hugoniot relations. See that paper for details
#  of equations and associated quantities. Briefly,
#     R-H1:             g0 * n0 * b0  =  g2 * n2 * b2
#     R-H2:  g0^2 * w0 * b0^2  +  P0  =  g2^2 * w2 * b2^2  +  P2  + Q_px
#     R-H3:       g0^2 * w0 * b0 * c  =  g2^2 * w2 * b2 * c       + Q_en
#  where
#     w    = E_rm  +  E_ke  +  P,  <--- enthalpy as total energy density + pressure
#     E_rm = n * m * c^2           <--- rest mass energy density
#     E_ke = n * m * c^2 * e(p)    <--- kinetic energy density, with e(p) =  √( 1 + (p/mc)^2 )  -  1
#     P    = 1/3 * n * p * v       <--- pressure
#
#  For closure, it is assumed that the two escaping fluxes are related by
#      Q_en = √[0.5*(1+β_Z)*c_cgs^2] * Q_px ,
#  i.e. the geometric mean of the arithmetic mean of u_Z and c. This
#  allows the solution to smoothly join with the non-rel version.
#  Use only fluid component of fluxes, not fluid+EM, for now.
function q_esc_calcs_relativistic_parallel()
    # Factor relating Q_en and Q_px
    q_fac = √( 0.5 * (1 + β_Z) * c_cgs^2 )

    # Calculate thermal pressure of far upstream gas
    pressure_Z  = dot(denZ_ion, tZ_ion) * kB_cgs
    ρ_Z         = dot(denZ_ion, aa_ion) * mp_cgs
    mask = (aa_ion .≥ 1)
    density_electron = dot(denZ_ion[mask], zz_ion[mask])

    # If electrons were not a separate species, add them in here
    if ! sc_electron
        pressure_Z  +=  density_electron * kB_cgs*tZ_electron
        ρ_Z         +=  density_electron * me_cgs
    end

    # Two terms to simplify the calculation of pressure_2.   #assumecold
    F_px_UpS_fl     = γ_Z^2 * β_Z^2 * (ρ_Z*c_cgs^2 + 5//2*pressure_Z) + pressure_Z
    F_energy_UpS_fl = γ_Z^2 * u_Z   * (ρ_Z*c_cgs^2 + 5//2*pressure_Z)
    term_aux = γ_2^2 * (q_fac * β_2^2  -  u_2)

    # Calculate far DwS density and pressure
    ρ_2 = ρ_Z * γ_Z*β_Z / γβ_2
    pressure_2 = (q_fac * F_px_UpS_fl - F_energy_UpS_fl - term_aux*ρ_2*c_cgs^2) / (q_fac + γ_fac*term_aux)

    # Calculate Q_px & Q_en (the physical escaping momentum & energy fluxes)
    Q_px = F_px_UpS_fl  -  (γ_2*β_2)^2 * ( ρ_2 * c_cgs^2  +  γ_fac * pressure_2 ) -  pressure_2
    Q_en = Q_px * q_fac

    # Convert Q_en & Q_px into q_esc_cal_**, which involves dividing by the
    # far UpS values. Subtract off mass-energy flux from F_energy_UpS_fl to
    # bring results in line with non-rel calculation.
    q_esc_cal_px = Q_px / F_px_UpS_fl
    q_esc_cal_energy = Q_en / ( F_energy_UpS_fl - γ_Z * u_Z * ρ_Z*c_cgs^2 )

end

function q_esc_calcs_nonrelativistic_oblique()
    error("Not implemented for oblique shocks yet.")
end

function q_esc_calcs_relativistic_oblique()
    error("Not implemented for oblique shocks yet.")
end
