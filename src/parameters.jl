module parameters

export na_particles, na_ions, na_grid, na_c, na_itrs
export psd_max, num_therm_bins
export na_cr, na_photons
export β_rel_fl, energy_rel_pt

# Parameters governing array sizes. All of these should be interpreted as
# maximum values for their respective quantity; the code will quite
# happily run if not all the array is used.

const na_particles = 100_000    # Max number of particles at each pcut
const na_ions = 5               # Max number of different ion species
const na_grid = 110             # Max number of elements in grid arrays
const na_c = 100                # Max number elements in pcut array
const na_itrs = 50              # Max number of iterations

# Max # of bins usable for phase space distribution and associated calculations.
# Applies to both momentum and angular dimensions
const psd_max = 200

const num_therm_bins = 150 # number of bins in thermal dist

# Max # of thermal particle crossings to hold before writing to scratch file
const na_cr = 10 * na_particles

const na_photons = 300  # Max size of photon array

# Cutoffs for nonrelativistic vs relativistic equations, for fluid & for particles
const β_rel_fl = 0.02
const energy_rel_pt = 0.005
end
