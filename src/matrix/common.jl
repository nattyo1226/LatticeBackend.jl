function build_matrix(
    ::Type{T},
    num_sites::Int,
    ids::Vector{Int},
    matrices::Vector{<:AbstractMatrix},
) where {T<:Number}
    if length(ids) != length(matrices)
        throw(ArgumentError("Length of ids and matrices must be the same"))
    end

    matrix_result = one(T)
    identity_matrix = build_matrix(T, Identity())

    for id in 1:num_sites
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
    ::Type{T},
    num_sites::Int,
    id::Int,
    matrix::AbstractMatrix,
) where T<:Number
    return build_matrix(T, num_sites, [id], [matrix])
end
