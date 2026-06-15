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

function build_matrix(
    ::Type{T},
    num_sites::Int,
    ids::Vector{Int},
    prs::Vector{<:AbstractOperatorPrimitive},
) where {T<:Number}
    matrices = [build_matrix(T, pr) for pr in prs]
    return build_matrix(T, num_sites, ids, matrices)
end

function build_matrix(
    num_sites::Int,
    ids::Vector{Int},
    prs::Vector{<:AbstractOperatorPrimitive},
)
    return build_matrix(ComplexF64, num_sites, ids, prs)
end

function build_matrix(
    ::Type{T},
    num_sites::Int,
    id::Int,
    pr::AbstractOperatorPrimitive,
) where T<:Number
    return build_matrix(T, num_sites, [id], [pr])
end

function build_matrix(
    num_sites::Int,
    id::Int,
    pr::AbstractOperatorPrimitive,
)
    return build_matrix(ComplexF64, num_sites, id, pr)
end

function build_matrix(::Type{T}, num_sites::Int, op::TensoredOperator) where T<:Number
    if iszero(op.coeff)
        dim_system = 2 ^ num_sites
        return zeros(T, dim_system, dim_system)
    end

    identity_matrix = build_matrix(T, Identity())
    mats = map(1:num_sites) do id
        idx = findfirst(pr -> pr.id == id, op.prs)
        idx === nothing ? identity_matrix : build_matrix(T, op.prs[idx].pr)
    end

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
