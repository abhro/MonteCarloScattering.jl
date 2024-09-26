module debug
#export energy_density, therm_energy_density, zone_vol, energy_pool

energy_density = zeros(n_grid, n_ions)
therm_energy_density = zeros(n_grid, n_ions)

zone_vol = zeros(n_grid)
energy_pool = zeros(n_grid)
end
