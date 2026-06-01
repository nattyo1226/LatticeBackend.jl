function build_matrix(::Type, pr::AbstractOperatorPrimitive)
    throw(ArgumentError("Unsupported operator type: $(typeof(pr))"))
end

function build_matrix(pr::AbstractOperatorPrimitive)
    return build_matrix(ComplexF64, pr)
end

function build_matrix(::Type{T}, ::Identity) where T<:Number
    return reshape(T[1, 0, 0, 1], (2, 2))
end

function build_matrix(::Type{T}, ::PauliX) where T<:Number
    return reshape(T[0, 1, 1, 0], (2, 2))
end

function build_matrix(::Type{T}, ::PauliY) where T<:Number
    return reshape(T[0, -im, im, 0], (2, 2))
end

function build_matrix(::Type{T}, ::PauliZ) where T<:Number
    return reshape(T[1, 0, 0, -1], (2, 2))
end
