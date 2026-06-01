function build_matrix(::Type, pr::AbstractOperatorPrimitive)
    throw(ArgumentError("Unsupported operator type: $(typeof(pr))"))
end

function build_matrix(pr::AbstractOperatorPrimitive)
    return build_matrix(ComplexF64, pr)
end

function build_matrix(::Type{T}, ::Identity) where T<:Number
    return T[1 0; 0 1]
end

function build_matrix(::Type{T}, ::PauliX) where T<:Number
    return T[0 1; 1 0]
end

function build_matrix(::Type{T}, ::PauliY) where T<:Number
    return T[0 -1im; 1im 0]
end

function build_matrix(::Type{T}, ::PauliZ) where T<:Number
    return T[1 0; 0 -1]
end
