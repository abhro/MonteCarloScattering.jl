module cosmo_calc
using Cosmology
using Roots
using QuadGK

#export H0, h, Ω_r, Ω_v, Ω_m, Ω_k, c, d_H, t_H
export cosmo

export get_redshift


# Parameters; cosmology values pulled from Planck 2013 results
const H0    = 67.8            # Hubble constant, km/s Mpc^-1
const h     = H0/100          # dimensionless Hubble constant
const Ω_r   = 0.4165 / H0^2   # Fraction of density in neutrinos?
const Ω_vac = 0.683 - 0.5*Ω_r # Fraction of density in dark energy and in matter
const Ω_m   = 0.317 - 0.5*Ω_r
const Ω_k   = 0.0             # Assume flat Universe, i.e. Ω_k = 1 - ∑(Ω_*)
const c     = 2.99792458e5    # Speed of light, km/s
const d_H   = c / H0          # Hubble distance, Mpc
const t_H   = 9.778e11 / H0   # Hubble time, years (9.778e11 = 1 s⋅Mpc/yr⋅km)

const cosmo = cosmology(h = h, OmegaK = Ω_k, OmegaR = Ω_r, OmegaM = Ω_m)

function E(c, z)
    Ω_k = 1 - (c.Ω_Λ + c.Ω_m + c.Ω_r)
    return √(c.Ω_r * (1 + z)^4 + c.Ω_m * (1 + z)^3 + Ω_k * (1 + z)^2 + c.Ω_Λ)
end
# use constants from within module
E(z) = √(Ω_r * (1 + z)^4 + Ω_m * (1 + z)^3 + Ω_k * (1 + z)^2 + Ω_vac)
D_C(z) = d_H * quadgk(t -> 1/E(t), 0, z)[begin] # Hogg (1999) Equation 15
D_C′(z) = d_H / E(z)

#using PolynomialRoots: roots

"""
    get_redshift(d_CM)

Calculator to get redshift from comoving distance.
Subroutine adapted from Hogg (1999). For more detail see Hogg (1999)
[DOI 10.48550/arXiv.astro-ph/9905116]

Note one major difference: Equation (13) in Hogg uses ``Ω_r`` where this
program uses ``Ω_k``, and does not include a term for what this program
calls ``Ω_r``. For more information, see Wright (2006): [DOI 10.1086/510102]

The program calculates lookback time (in years) as well, but since this result
is not necessary for GRB redshift calculation it is not passed back to the
calling subroutine.

### Arguments
- `d_CM`: Comoving distance in Mpc
### Returns
- `z`: Redshift
"""
function get_redshift(d_CM)
    d_CM ≤ 0 && throw(DomainError(d_CM, "d_CM must be positive"))

    # Step through integral in z (Equation (15) in Hogg 1999) until next step exceeds d_CM
    if d_CM < 0.443
        z      = 0.0
        #t_look = d_CM * 3.2616e6 # convert Mpc to years
        # d_CM less than critical value. Skip stepping through z integral and
        # integration to find t_look
        code_stat = 2
        return z
    end

    return find_zero(
        ((z -> D_C(z) - d_CM), D_C′),   # function-derivative pair
        0,                              # initial guess
        Roots.Newton(),                 # method to use
    )
end
end # module
