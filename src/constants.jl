module constants

export mₚ_cgs, mₑ_cgs, qₚ_cgs, kB_cgs, c_cgs, h_cgs, ħ_cgs
export E₀_proton, E₀_electron
export B_CMB0, T_CMB0
export T_th, M_res, Γ_res, rmp, rmpi
export Tₜₕ

using PhysicalConstants: CODATA2018
using Unitful
using Unitful: g, K, cm, s, erg, GeV

# Physical or arithmetic constants

const mₚ_cgs = ustrip(g,     CODATA2018.ProtonMass)           # Proton mass in grams
const mₑ_cgs = ustrip(g,     CODATA2018.ElectronMass)         # Electron mass in grams
const qₚ_cgs = 4.803205e-10 # Proton charge in cgs ESU units: statC or Fr
const kB_cgs = ustrip(erg/K, CODATA2018.BoltzmannConstant)    # Boltzmann's constant [erg/K]
const c_cgs  = ustrip(cm/s,  CODATA2018.SpeedOfLightInVacuum) # speed of light [cm/s]
const h_cgs  = ustrip(erg*s, CODATA2018.PlanckConstant)       # Planck constant [erg⋅s]
const ħ_cgs  = ustrip(erg*s, CODATA2018.ReducedPlanckConstant)

const E₀_proton   = mₚ_cgs * c_cgs^2 # proton rest mass energy [erg]
const E₀_electron = mₑ_cgs * c_cgs^2 # electron rest mass energy [erg]

const B_CMB0 = 3.27e-6  # Equivalent B field[G] to CMB energy density at a redshift of 0
const T_CMB0 = 2.725    # Temperature[K] of CMB at a redshift of 0

const T_th  = 0.2797                    # Threshold kinetic energy, GeV, for pion production
const M_res = 1.1883                    # Resonance mass, GeV
const Γ_res = 0.2264                    # Resonance width, GeV
const rmp   = ustrip(GeV, E₀_proton*erg)# Proton rest mass, GeV
const rmpi  = 0.134976                  # Neutral pion rest mass, GeV

const Tₜₕ = T_th # alias

end
