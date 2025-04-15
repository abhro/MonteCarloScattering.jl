function populate_ε_target!(ε_target, z_max, γ_sf_grid, uₓ_sk_grid, u₀, γ₀, energy_transfer_frac)
    for i in eachindex(ε_target)
        if uₓ_sk_grid[i] != u₀
            z_curr = γ₀ * u₀ / (γ_sf_grid[i] * uₓ_sk_grid[i])
            ε_target[i] = energy_transfer_frac * (z_curr - 1) / (z_max - 1)
        end
    end
    return ε_target
end

# A little bit of type piracy here to make things go smoother
function Base.fill!(arrs::NTuple{N,<:AbstractArray{T}}, x::T) where {N,T}
    for arr in arrs
        fill!(arr, x)
    end
end
