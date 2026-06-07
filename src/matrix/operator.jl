function build_matrix(::Type, ::Lattice, op::AbstractOperator)
    throw(ArgumentError("Unsupported operator type: $(typeof(op))"))
end

function build_matrix(lattice::Lattice, op::AbstractOperator)
    return build_matrix(ComplexF64, lattice, op)
end

function build_matrix(::Type{T}, lattice::Lattice, op::OnsiteOperator{P}) where {T<:Number, P<:AbstractOperatorPrimitive}
    num_sites = Int(nsites(lattice))

    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    id = op.id
    matrix = build_matrix(T, P)
    return T(op.coeff) * build_matrix(T, num_sites, id, matrix)
end

function build_matrix(::Type{T}, lattice::Lattice, op::UniformOnsiteOperator{P}) where {T<:Number, P<:AbstractOperatorPrimitive}
    num_sites = Int(nsites(lattice))
    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    dim_system = 2 ^ num_sites
    matrix = build_matrix(T, P)
    matrix_result = zeros(T, dim_system, dim_system)

    for id in 1:num_sites
        matrix_result += build_matrix(T, num_sites, id, matrix)
    end

    return T(op.coeff) * matrix_result
end

function build_matrix(::Type{T}, lattice::Lattice, op::PairOperator{P1, P2}) where {T<:Number, P1<:AbstractOperatorPrimitive, P2<:AbstractOperatorPrimitive}
    num_sites = Int(nsites(lattice))

    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    id1 = op.id1
    id2 = op.id2
    matrix1 = build_matrix(T, P1)
    matrix2 = build_matrix(T, P2)
    return T(op.coeff) * build_matrix(T, num_sites, (id1, id2), (matrix1, matrix2))
end

function build_matrix(::Type{T}, lattice::Lattice, op::UniformPairOperator{P1, P2}) where {T<:Number, P1<:AbstractOperatorPrimitive, P2<:AbstractOperatorPrimitive}
    num_sites = Int(nsites(lattice))

    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    dim_system = 2 ^ num_sites
    matrix1 = build_matrix(T, P1)
    matrix2 = build_matrix(T, P2)
    matrix_result = zeros(T, dim_system, dim_system)

    for (id1, id2) in neighbor_pairs(lattice)
        matrix_result += build_matrix(T, num_sites, (id1, id2), (matrix1, matrix2))
    end

    return T(op.coeff) * matrix_result
end

function build_matrix(::Type{T}, lattice::Lattice, op::SummedOperator) where T<:Number
    num_sites = Int(nsites(lattice))
    dim_system = 2 ^ num_sites
    matrix_result = zeros(T, dim_system, dim_system)

    for term in op.ops
        if iszero(term.coeff)
            continue
        end
        matrix_result += build_matrix(T, lattice, term)
    end

    return matrix_result
end
