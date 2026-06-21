function build_matrix(
    ::Type{N},
    ids_all::AbstractVector{<:AbstractIndex{T}},
    ids::AbstractVector{<:AbstractIndex{T}},
    matrices::AbstractVector{<:AbstractMatrix},
) where {N<:Number,T<:AbstractSystemTag}
    if length(ids) != length(matrices)
        throw(ArgumentError("Length of ids and matrices must be the same"))
    end

    matrix_result = one(N)
    identity_matrix = build_matrix(N, Identity{T}())

    for id in ids_all
        if id in ids
            matrix = matrices[findfirst(x -> x == id, ids)]
            matrix_result = kron(matrix_result, matrix)
        else
            matrix_result = kron(matrix_result, identity_matrix)
        end
    end

    return matrix_result
end

function build_matrix(
    ::Type{N},
    space::Space{T},
    ids::AbstractVector{<:AbstractIndex{T}},
    matrices::AbstractVector{<:AbstractMatrix},
) where {N<:Number,T<:AbstractSystemTag}
    return build_matrix(N, indices(space), ids, matrices)
end

function build_matrix(
    ::Type{N},
    space::Space{T},
    id::AbstractIndex{T},
    matrix::AbstractMatrix,
) where {N<:Number,T<:AbstractSystemTag}
    return build_matrix(N, indices(space), [id], [matrix])
end
