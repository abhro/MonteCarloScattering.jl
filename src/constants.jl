module constants

export E₀_proton, E₀_electron
export B_CMB0, T_CMB0
export T_th, M_res, Γ_res, rmpi
export Tₜₕ
export rad_loss_fac

using PhysicalConstants.CODATA2018: σ_e as σ_T
using Unitful, UnitfulAstro, UnitfulEquivalences
using Unitful: g, K, cm, s, erg, keV, GeV
using Unitful: mp, me, c, q, k as kB, h, ħ    # physical constants
using UnitfulGaussian: Fr, G

# Physical or arithmetic constants

"Proton rest energy"
const E₀_proton   = mp * c^2 |> erg
"Electron rest energy"
const E₀_electron = me * c^2 |> erg

"Equivalent B field to CMB energy density at a redshift of 0"
const B_CMB0 = 3.27e-6G
"Temperature of CMB at a redshift of 0"
const T_CMB0 = 2.725K

"Threshold kinetic energy, for pion production"
const T_th  = 0.2797GeV
"Resonance mass"
const M_res = 1.1883GeV
"Resonance width"
const Γ_res = 0.2264GeV
"Neutral pion rest energy"
const rmpi  = 0.134976GeV

const Tₜₕ = T_th # alias

# The numerator 4/3 c σ_T comes from Rybicki and Lightman Eq. (6.7b)
# The denominator is included because β² U_B = (v²/c²) B²/8π, so we can use v
# and B instead of of β and U_B (I think). Not sure where the extra factor of c
# and the electron mass comes from.
"Factor related to radiative losses"
const rad_loss_fac = 4//3 * c * σ_T / (c^3 * me^2 * 8π) |> s^2/g^2
end
