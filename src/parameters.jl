module parameters

# Parameters governing array sizes. All of these should be interpreted as maximum values for
# their respective quantity; the code will quite happily run if not all the array is used.

export na_particles, na_c
const na_particles = 100_000    # Max number of particles at each pcut
const na_c = 100                # Max number elements in pcut array

export psd_max, num_therm_bins
# Maximum number of bins usable for phase space distribution and associated calculations.
# Applies to both momentum and angular dimensions
const psd_max = 200
# number of bins in thermal dist
const num_therm_bins = 150

export na_cr, na_photons
# Maximum number of thermal particle crossings to hold before writing to scratch file
const na_cr = 10 * na_particles
# Maximum size of photon array
const na_photons = 300

# Cutoffs for nonrelativistic vs relativistic equations, for fluid & for particles
export β_rel_fl, energy_rel_pt
const β_rel_fl = 0.02
const energy_rel_pt = 0.005
end
