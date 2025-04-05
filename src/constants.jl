module constants

export mₚ_cgs, mₑ_cgs, qₚ_cgs, kB_cgs, c_cgs, h_cgs, ħ_cgs
export E₀_proton, E₀_electron
export B_CMB0, T_CMB0
export T_th, M_res, Γ_res, rmp, rmpi
export Tₜₕ

using PhysicalConstants: CODATA2018
using Unitful, UnitfulGaussian
using Unitful: g, K, cm, s, erg, GeV
using UnitfulGaussian: Fr, G

# Physical or arithmetic constants

const mₚ_cgs = uconvert(g,     Unitful.mp)      # Proton mass in grams
const mₑ_cgs = uconvert(g,     Unitful.me)      # Electron mass in grams
const qₚ_cgs = 4.803205e-10Fr                   # Proton charge in cgs ESU units
const kB_cgs = uconvert(erg/K, Unitful.k)       # Boltzmann's constant
const c_cgs  = uconvert(cm/s,  convert(Quantity{Int128}, Unitful.c0)) # speed of light
# ^ Int128 used because usually it's c^2, and in cgs it overflows Int64
const h_cgs  = uconvert(erg*s, Unitful.h)       # Planck constant
const ħ_cgs  = uconvert(erg*s, Unitful.ħ)

const E₀_proton   = mₚ_cgs * float(c_cgs)^2 # proton rest mass energy [erg]
const E₀_electron = mₑ_cgs * float(c_cgs)^2 # electron rest mass energy [erg]

const B_CMB0 = 3.27e-6G  # Equivalent B field to CMB energy density at a redshift of 0
const T_CMB0 = 2.725K    # Temperature of CMB at a redshift of 0

const T_th  = 0.2797GeV                 # Threshold kinetic energy, for pion production
const M_res = 1.1883GeV                 # Resonance mass
const Γ_res = 0.2264GeV                 # Resonance width
const rmp   = uconvert(GeV, E₀_proton)  # Proton rest energy
const rmpi  = 0.134976GeV               # Neutral pion rest energy

const Tₜₕ = T_th # alias

end
