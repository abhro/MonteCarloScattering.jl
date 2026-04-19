module CGSTypes
using Unitful: g, K, cm, s, dyn, erg, keV, GeV
using UnitfulGaussian: Fr, G, qcgs

export LengthCGS, MassCGS, TimeCGS, SpeedCGS, MomentumCGS, BFieldCGS, EnergyCGS, PressureCGS,
    MomentumFluxCGS, MomentumDensityFluxCGS, EnergyFluxCGS, EnergyDensityFluxCGS

const LengthCGS = typeof(1.0 * cm)
const MassCGS = typeof(1.0 * g)
const TimeCGS = typeof(1.0 * s)
const AreaCGS = typeof(1.0 * cm^2)
const SpeedCGS = typeof(1.0 * cm / s)
const MomentumCGS = typeof(1.0 * g * cm / s)
const BFieldCGS = typeof(1.0 * G)
const EnergyCGS = typeof(1.0 * erg)
const EnergyDensityCGS = typeof(1.0 * erg / cm^3)
const PressureCGS = typeof(1.0 * erg / cm^3)
const MomentumFluxCGS = typeof(1.0 * g * cm^2 / s^2) # (g*cm/s) * cm/s
const MomentumDensityFluxCGS = typeof(1.0 * erg / cm^3)
const EnergyFluxCGS = typeof(1.0 * erg * cm / s)
const EnergyDensityFluxCGS = typeof(1.0 * erg / (cm^2 * s))
end
