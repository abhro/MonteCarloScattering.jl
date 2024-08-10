module debug
using ..parameters: na_grid, na_ions

#export energy_density, therm_energy_density, zone_vol, energy_pool

energy_density = zeros(na_grid, na_ions)
therm_energy_density = zeros(na_grid, na_ions)

zone_vol = zeros(na_grid)
energy_pool = zeros(na_grid)
end
