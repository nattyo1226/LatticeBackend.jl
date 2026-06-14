function build_matrix(
    ::Type{T},
    num_sites::Int,
    ids_site::NTuple{K,Int},
    matrices::NTuple{K,<:AbstractMatrix},
) where {T<:Number,K}
    matrix_result = one(T)
    identity_matrix = build_matrix(T, Identity)

    for id in 1:num_sites
        if id in ids_site
            matrix = matrices[findfirst(x -> x == id, ids_site)]
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
    id_site::Int,
    matrix::AbstractMatrix,
) where T<:Number
    return build_matrix(T, num_sites, (id_site,), (matrix,))
end
