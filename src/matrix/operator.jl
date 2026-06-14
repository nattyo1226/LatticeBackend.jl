function build_matrix(::Type, ::Int, op::AbstractOperator)
    throw(ArgumentError("Unsupported operator type: $(typeof(op))"))
end

function build_matrix(num_sites::Int, op::AbstractOperator)
    return build_matrix(ComplexF64, num_sites, op)
end

function build_matrix(::Type{T}, lattice::Lattice, op::AbstractOperator) where T<:Number
    num_sites = Int(nsites(lattice))
    return build_matrix(T, num_sites, op)
end

function build_matrix(lattice::Lattice, op::AbstractOperator)
    return build_matrix(ComplexF64, lattice, op)
end

function build_matrix(::Type{T}, num_sites::Int, op::TensoredOperator) where T<:Number
    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    mats_local = Dict{Int,AbstractMatrix}(
        pr.id => build_matrix(T, pr.pr)
        for pr in op.prs
    )

    mats = [
        get(mats_local, id, build_matrix(T, Identity()))
        for id in 1:num_sites
    ]

    return T(op.coeff) * reduce(kron, mats)
end

function build_matrix(::Type{T}, num_sites::Int, op::SummedOperator) where T<:Number
    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    ops_filtered = filter(op_sub -> !iszero(op_sub.coeff), op.ops)
    if isempty(ops_filtered)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    return sum(build_matrix(T, num_sites, op_sub) for op_sub in ops_filtered)
end
