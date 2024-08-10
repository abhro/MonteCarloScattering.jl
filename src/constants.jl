module constants

export mp_cgs, amu, me_cgs, qp_cgs, kB_cgs, c_cgs, h_cgs, ħ_cgs
export E₀_proton, E₀_electron
export erg2eV, erg2keV, erg2MeV, eV2erg, keV2erg, MeV2erg, Jy2cgs, pc2cm
export min_θ8
export B_CMB0, T_CMB0
export T_th, M_res, Γ_res, rmp, rmpi
export Tₜₕ

using PhysicalConstants: CODATA2018
using Unitful
using Unitful: g, K, cm, s, yr, erg, eV, keV, MeV
using UnitfulAstro: pc, Jy

# Physical or arithmetic constants

const mp_cgs = ustrip(g,     CODATA2018.ProtonMass)           # Proton mass in grams
const amu    = ustrip(g,     CODATA2018.AtomicMassConstant)   # Atomic mass unit in grams
const me_cgs = ustrip(g,     CODATA2018.ElectronMass)         # Electron mass in grams
const qp_cgs = 4.803205e-10 # Proton charge in cgs ESU units: statC or Fr
const kB_cgs = ustrip(erg/K, CODATA2018.BoltzmannConstant)    # Boltzmann's constant (cgs)
const c_cgs  = ustrip(cm/s,  CODATA2018.SpeedOfLightInVacuum) # speed of light (cm/s)
const h_cgs  = ustrip(erg*s, CODATA2018.PlanckConstant)       # Planck constant (erg-sec)
const ħ_cgs  = ustrip(erg*s, CODATA2018.ReducedPlanckConstant)

const E₀_proton    = mp_cgs*(c_cgs^2) # proton rest mass energy [erg]
const E₀_electron  = me_cgs*(c_cgs^2) # electron rest mass energy [erg]

const B_CMB0 = 3.27e-6  # Equivalent B field[G] to CMB energy density at a redshift of 0
const T_CMB0 = 2.725    # Temperature[K] of CMB at a redshift of 0

const erg2eV  = ustrip(eV,  1erg)       # conversion ergs to eV
const erg2keV = ustrip(keV, 1erg)       # conversion ergs to keV
const erg2MeV = ustrip(MeV, 1erg)       # conversion ergs to MeV
const eV2erg  = ustrip(erg, 1eV)        # conversion eV to ergs
const keV2erg = ustrip(erg, 1keV)       # conversion keV to ergs
const MeV2erg = ustrip(erg, 1MeV)       # conversion MeV to ergs
const Jy2cgs  = ustrip(erg/cm^2, 1Jy)   # Jansky in erg cm^(-2)
const pc2cm   = ustrip(cm, 1pc)         # cm in a parsec

const T_th  = 0.2797                    # Threshold kinetic energy, GeV, for pion production
const M_res = 1.1883                    # Resonance mass, GeV
const Γ_res = 0.2264                    # Resonance width, GeV
const rmp   = E₀_proton * erg2eV * 1e-9 # Proton rest mass, GeV
const rmpi  = 0.134976                  # Neutral pion rest mass, GeV

const Tₜₕ = T_th # alias

# Set the minimum value of θ such that cos(θ) can be distinguished from 1.0 in whichever precision
const min_θ8 = 10 * √(2eps(Float64))
end
