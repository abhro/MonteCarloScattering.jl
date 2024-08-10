module parameters

export na_particles, na_ions, na_grid, na_c, na_itrs
export psd_max, num_therm_bins
export na_cr, na_photons
export β_rel_fl, energy_rel_pt

# Parameters governing array sizes. All of these should be interpreted as
# maximum values for their respective quantity; the code will quite
# happily run if not all the array is used.

const na_particles = 100_000    # Max no. of particles at each pcut
const na_ions = 5               # Max no. of different ion species
const na_grid = 110                # Max no. of elements in grid arrays
const na_c = 100                # Max no. elements in pcut array
const na_itrs = 50              # Max no. of iterations

# Max # of bins usable for phase space distribution and associated calcs.
# Applies to both momentm and angular dimensions
const psd_max = 200

const num_therm_bins = 150 # number of bins in thermal dist

# Max # of thermal particle crossings to hold before writing to scratch file
const na_cr = 10 * na_particles

const na_photons = 300  # Max size of photon array

# Cutoffs for non-rel vs rel eqs, for fluid & for particles
const β_rel_fl = 0.02
const energy_rel_pt = 0.005
end
