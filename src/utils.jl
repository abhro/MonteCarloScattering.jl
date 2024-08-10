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

function geometric_center!(x, y)
    n = length(x)
    n == length(y) - 1 || throw(DimensionMismatch())

    for i in eachindex(y[begin:end-1])
        x[i] = √(y[i] * y[i+1])
    end
end

function geometric_center(y)
    x = zeros(eltype(y), length(y) - 1)

    for i in eachindex(y[begin:end-1])
        x[i] = √(y[i] * y[i+1])
    end
end
