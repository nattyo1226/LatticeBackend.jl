function build_matrix(::Type, ::Space{T}, op::AbstractOperator{T}) where {T<:AbstractSystemTag}
    throw(ArgumentError("Unsupported operator type: $(typeof(op))"))
end

function build_matrix(space::Space{T}, op::AbstractOperator{T}) where {T<:AbstractSystemTag}
    return build_matrix(ComplexF64, space, op)
end

function build_matrix(
    ::Type{N},
    space::Space{T},
    ids::AbstractVector{<:AbstractIndex{T}},
    prs::AbstractVector{<:AbstractOperatorPrimitive{T}},
) where {N<:Number,T<:AbstractSystemTag}
    matrices = [build_matrix(N, pr) for pr in prs]
    return build_matrix(N, space, ids, matrices)
end

function build_matrix(
    space::Space{T},
    ids::AbstractVector{<:AbstractIndex{T}},
    prs::AbstractVector{<:AbstractOperatorPrimitive{T}},
) where {T<:AbstractSystemTag}
    return build_matrix(ComplexF64, space, ids, prs)
end

function build_matrix(
    ::Type{N},
    space::Space{T},
    id::AbstractIndex{T},
    pr::AbstractOperatorPrimitive{T},
) where {N<:Number,T<:AbstractSystemTag}
    return build_matrix(N, space, [id], [pr])
end

function build_matrix(
    space::Space{T},
    id::AbstractIndex{T},
    pr::AbstractOperatorPrimitive{T},
) where {T<:AbstractSystemTag}
    return build_matrix(ComplexF64, space, [id], [pr])
end

function build_matrix(
    ::Type{N},
    ids_all::AbstractVector{<:AbstractIndex{T}},
    op::TensoredOperator{T},
) where {N<:Number,T<:AbstractSystemTag}
    if iszero(op.coeff)
        dim_system = 2 ^ length(ids_all)
        return zeros(N, dim_system, dim_system)
    end

    identity_matrix = build_matrix(N, Identity{T}())
    z_matrix = build_matrix(N, PauliZ{T}())

    mats = map(ids_all) do id
        mats_tmp = map(op.prs) do pr
            grade = majorana_grade(pr.pr)

            if id < pr.id
                isodd(grade) ? z_matrix : identity_matrix
            elseif id == pr.id
                build_matrix(N, pr.pr)
            else
                identity_matrix
            end
        end

        foldl(*, mats_tmp; init=identity_matrix)
    end

    return N(op.coeff) * reduce(kron, mats)
end

function build_matrix(
    ::Type{N},
    space::Space{T},
    op::TensoredOperator{T},
) where {N<:Number,T<:AbstractSystemTag}
    return build_matrix(N, indices(space), op)
end

function build_matrix(
    ::Type{N},
    ids_all::AbstractVector{<:AbstractIndex{T}},
    op::SummedOperator{T},
) where {N<:Number,T<:AbstractSystemTag}
    ops_filtered = filter(op_sub -> !iszero(op_sub.coeff), op.ops)

    if isempty(ops_filtered)
        dim_system = 2 ^ length(ids_all)
        return zeros(N, dim_system, dim_system)
    end

    return sum(build_matrix(N, ids_all, op_sub) for op_sub in ops_filtered)
end

function build_matrix(
    ::Type{N},
    space::Space{T},
    op::SummedOperator{T},
) where {N<:Number,T<:AbstractSystemTag}
    return build_matrix(N, indices(space), op)
end
