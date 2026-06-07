function build_matrix(::Type{T}, ::Type{P}) where {T<:Number, P<:AbstractOperatorPrimitive}
    throw(ArgumentError("Unsupported operator type: $(typeof(pr))"))
end

function build_matrix(pr::AbstractOperatorPrimitive)
    return build_matrix(ComplexF64, pr)
end

function build_matrix(::Type{T}, ::Type{Identity}) where T<:Number
    return reshape(T[1, 0, 0, 1], (2, 2))
end

function build_matrix(::Type{T}, ::Type{PauliX}) where T<:Number
    return reshape(T[0, 1, 1, 0], (2, 2))
end

function build_matrix(::Type{T}, ::Type{PauliY}) where T<:Number
    return reshape(T[0, -im, im, 0], (2, 2))
end

function build_matrix(::Type{T}, ::Type{PauliZ}) where T<:Number
    return reshape(T[1, 0, 0, -1], (2, 2))
end
