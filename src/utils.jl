"""
Functions that might be nice to have when the code base is refactored
"""

using Unitful: g, K, cm
using UnitfulGaussian: Fr

"""
    adjacent_apply(f, x)

Return `f(x[i], x[i+1])` for each `i` in `x`'s index, excluding the last `i`.
"""
function adjacent_apply(f, x::AbstractVector{T}) where {T}
    y = Vector{T}(undef, length(x) - 1)
    adjacent_apply!(f, y, x)
    return y
end

"""
    adjacent_apply!(f, y, x)

In-place version of `adjacent_apply`.
"""
function adjacent_apply!(f, y, x)
    length(y) == length(x) - 1 || throw(DimensionMismatch("y must have length one less than x"))
    for i in eachindex(@view(x[begin:(end - 1)]), y)
        y[i] = f(x[i], x[i + 1])
    end

    return y
end

#geometric_center!(x, y) = adjacent_apply!((a, b) -> √(a*b), x, y)
function geometric_center!(x, y)
    n = length(x)
    n == length(y) - 1 || throw(DimensionMismatch())

    for i in eachindex(y[begin:(end - 1)])
        x[i] = √(y[i] * y[i + 1])
    end
    return
end

#geometric_center(y) = adjacent_apply((a, b) -> √(a*b), y)
function geometric_center(y)
    x = zeros(eltype(y), length(y) - 1)

    for i in eachindex(@view(y[begin:(end - 1)]))
        x[i] = √(y[i] * y[i + 1])
    end
    return
end

using Unitful
lorentz(v::Unitful.Velocity) = lorentz(v / c)
"""
    γ = lorentz(v)
    γ = lorentz(β)

Get Lorentz factor γ from velocity v or β (in units of c).
"""
lorentz(β::Real) = 1 / √(1 - β^2)

"""
    β(γ)

Get velocity β (in units of c) from Lorentz factor γ.
"""
β(γ) = √(1 - 1 / γ^2)


@kwdef struct Species
    mass::typeof(1.0g)
    charge::typeof(1.0Fr)
    temperature::typeof(1.0K)
    number_density::typeof(1.0/cm^3)
end
function Base.getproperty(s::Species, sym::Symbol)
    if sym == :m
        sym = :mass
    elseif sym == :q || sym == :Z
        sym = :charge
    elseif sym == :T || sym == :temp
        sym = :temperature
    elseif sym == :n || sym == :ρ_N || sym == :density
        sym = :number_density
    end
    return getfield(s, sym)
end
Base.show(io::IO, s::Species) =
    print(io, "Species(m = ", s.m, ", q = ", s.q, ", T = ", s.T, ", n = ", s.n, ")")
mass(s::Species) = s.mass
charge(s::Species) = s.charge
temperature(s::Species) = s.temperature
density(s::Species) = s.number_density
number_density(s::Species) = s.number_density
