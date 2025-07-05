module debug
export energy_density, therm_energy_density, zone_vol, energy_pool

import ..n_grid, ..n_ions

const energy_density = zeros(n_grid, n_ions)
const therm_energy_density = zeros(n_grid, n_ions)

const zone_vol = zeros(n_grid)
const energy_pool = zeros(n_grid)
end
