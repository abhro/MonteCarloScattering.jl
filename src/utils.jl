"""
Functions that might be nice to have when the code base is refactored
"""

include("constants.jl")
using .constants: c_cgs


"""
Return f(x[i], x[i+1]) for each i in x's index, excluding the last i.
"""
function adjacent_apply(f, x)
    y = Vector{eltype(x)}(undef, length(x) - 1)
    adjacent_apply!(f, y, x)
    return y
end

function adjacent_apply!(f, y, x)
    length(y) == length(x) - 1 || throw(DimensionMismatch("y must have length one less than x"))
    for i in eachindex(x[begin:end-1], y)
        y[i] = f(x[i], x[i+1])
    end

    return y
end

#geometric_center!(x, y) = adjacent_apply!((a, b) -> √(a*b), x, y)
function geometric_center!(x, y)
    n = length(x)
    n == length(y) - 1 || throw(DimensionMismatch())

    for i in eachindex(y[begin:end-1])
        x[i] = √(y[i] * y[i+1])
    end
end

#geometric_center(y) = adjacent_apply((a, b) -> √(a*b), y)
function geometric_center(y)
    x = zeros(eltype(y), length(y) - 1)

    for i in eachindex(y[begin:end-1])
        x[i] = √(y[i] * y[i+1])
    end
end

"""
Get Lorentz factor γ from velocity β (in units of c)
"""
γ(β) = 1 / √(1 - β^2)
"""
Get velocity β (in units of c) from Lorentz factor γ
"""
β(γ) = √(1 - 1/γ^2)

struct RelativisticVelocity{T}
    u::T # TODO parametrize u differently because it can have units
    β::T
    γ::T
end
function RelativisticVelocity(u::T) where T
    u < c_cgs || throw(DomainError(u, "speed is greater than speed of light"))

    β = u/c_cgs
    return RelativisticVelocity(u, β, γ(β))
end
function Base.show(io::IO, v::RelativisticVelocity{T}) where T
    print(io, "RelativisticVelocity{", T, "}")
    print(io, "(")
    print(io, "u=", v.u, ", ")
    print(io, "β=", v.β, ", ")
    print(io, "γ=", v.γ)
    print(io, ")")
end

function velocity_from_β(β::T) where T
    return RelativisticVelocity{T}(β*c_cgs, β, γ(β))
end

function velocity_from_γ(γ::T) where T
    β_from_γ = β(γ)
    return RelativisticVelocity{T}(β_from_γ*c_cgs, β_from_γ, γ)
end
