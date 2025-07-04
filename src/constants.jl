module constants

export qₚ_cgs
export E₀_proton, E₀_electron
export B_CMB0, T_CMB0
export T_th, M_res, Γ_res, rmpi
export Tₜₕ

using PhysicalConstants: CODATA2018
using Unitful, UnitfulAstro, UnitfulGaussian, UnitfulEquivalences
using Unitful: g, K, cm, s, erg, keV, GeV
using Unitful: mp, me, c, q, k as kB, h, ħ    # physical constants
using UnitfulGaussian: Fr, G

# Physical or arithmetic constants

const qₚ_cgs = uconvert(Fr, q, ChargeEquivalence())   # Proton charge in cgs ESU units

"Proton rest mass energy"
const E₀_proton   = mp * c^2 |> erg
"Electron rest mass energy"
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

end
