function build_zero_matrix(::Type{T}, num_sites::Int) where T<:Number
    dim = 2 ^ num_sites
    return spzeros(T, (dim, dim))
end

function build_zero_matrix(::Type{T}, lattice::Lattice) where T<:Number
    num_sites = Int(nsites(lattice))
    return build_zero_matrix(T, num_sites)
end

function build_matrix(::Type{T}, num_sites::Int, ids_site::NTuple{K,Int}, matrices::NTuple{K,<:AbstractMatrix}) where {T<:Number,K}
    dim = 2 ^ num_sites
    dim_matrix = 2 ^ K

    rows = Int[]
    cols = Int[]
    vals = T[]

    for col in 0:(dim-1)
        col_local = 0
        for (n, id) in enumerate(ids_site)
            bit = (col >> (id - 1)) & 1
            col_local |= bit << (n - 1)
        end

        for row_local in 0:(dim_matrix-1)
            val = one(T)
            for (n, matrix) in enumerate(matrices)
                row_bit = (row_local >> (n - 1)) & 1
                col_bit = (col_local >> (n - 1)) & 1
                val *= matrix[row_bit+1, col_bit+1]
            end

            iszero(val) && continue

            row = col
            for (n, id) in enumerate(ids_site)
                row_bit = (row_local >> (n - 1)) & 1
                mask = 1 << (id - 1)
                row = (row & ~mask) | (row_bit << (id - 1))
            end

            push!(rows, row + 1)
            push!(cols, col + 1)
            push!(vals, val)
        end
    end

    return sparse(rows, cols, vals)
end

function build_matrix(::Type{T}, num_sites::Int, id_site::Int, matrix::AbstractMatrix) where T<:Number
    return build_matrix(T, num_sites, (id_site,), (matrix,))
end
