module cosmo_calc
using Cosmology
using Roots
using QuadGK

export cosmo, get_redshift

# Parameters; cosmology values pulled from Planck 2013 results
const h = 0.678                 # dimensionless Hubble constant
const Ωᵣ = 0.4165 / (h * 100)^2 # Fraction of density in neutrinos?
const Ω_vac = 0.683 - 0.5 * Ωᵣ  # Fraction of density in dark energy and in matter
const Ωₘ = 0.317 - 0.5 * Ωᵣ
const Ωₖ = 0.0                  # Assume flat Universe, i.e. Ωₖ = 1 - ∑(Ω_*)
const cosmo = cosmology(; h, OmegaK = Ωₖ, OmegaR = Ωᵣ, OmegaM = Ωₘ)
const d_H = hubble_dist(cosmo, 0)       # Hubble distance at z=0, Mpc

D_C′(z) = d_H / Cosmology.E(cosmo, z)

"""
    get_redshift(d_CM)

Calculator to get redshift from comoving distance.
Subroutine adapted from Hogg (1999). For more detail see Hogg (1999)
[DOI 10.48550/arXiv.astro-ph/9905116]

### Arguments
- `d_CM`: Comoving distance in Mpc

### Returns
- `z`: Redshift
"""
function get_redshift(d_CM)
    d_CM ≤ 0 && throw(DomainError(d_CM, "d_CM must be positive"))

    # d_CM less than critical value. Skip stepping through z integral and
    # integration to find t_look
    if d_CM < 0.443
        z = 0.0
        return z
    end

    return find_zero(
        (
            (z -> comoving_radial_dist(cosmo, z) - d_CM),   # function
            D_C′                                            # derivative
        ),
        0,                              # initial guess
        Roots.Newton(),                 # method to use
    )
end
end # module
