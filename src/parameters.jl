"""
Parameters governing array sizes. All of these should be interpreted as maximum values for
their respective quantity; the code will quite happily run if not all the array is used.
"""
module parameters

export na_particles, na_c
"Max number of particles at each pcut"
const na_particles = 100_000
"Max number elements in pcut array"
const na_c = 100

export psd_max, num_therm_bins
"""
Maximum number of bins usable for phase space distribution and associated calculations.
Applies to both momentum and angular dimensions
"""
const psd_max = 200
"Number of bins in thermal distribution"
const num_therm_bins = 150

export na_cr, na_photons
"Maximum number of thermal particle crossings to hold before writing to scratch file"
const na_cr = 10 * na_particles
"Maximum size of photon array"
const na_photons = 300

export β_rel_fl, E_rel_pt
"Cutoff for nonrelativistic vs relativistic equations, for fluid"
const β_rel_fl = 0.02
"Cutoff for nonrelativistic vs relativistic equations, for particles"
const E_rel_pt = 0.005
end
