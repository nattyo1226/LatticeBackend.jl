function build_matrix(::Type{T}, ::AbstractOperatorPrimitive) where {T<:Number}
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

function build_matrix(::Type{T}, pr::SummedOperatorPrimitive) where T<:Number
    return sum(build_matrix(T, pr_sub) for pr_sub in pr.prs)
end

function build_matrix(::Type{T}, pr::ProductedOperatorPrimitive) where T<:Number
    return foldl(*, [build_matrix(T, pr_sub) for pr_sub in pr.prs]; init=I)
end
