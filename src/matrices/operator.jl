function apply(pr::IndexedOperatorPrimitive{T}, state::State) where {T<:AbstractSystemTag}
    id_bit = to_bit(pr.id)
    return apply(pr.pr, state, id_bit)
end

function apply(pr::IndexedOperatorPrimitive{T}, states::Vector{State}) where {T<:AbstractSystemTag}
    id_bit = to_bit(pr.id)
    return apply(pr.pr, states, id_bit)
end

function apply(pr::IndexedOperatorPrimitive{T}, bits::Int) where {T<:AbstractSystemTag}
    id_bit = to_bit(pr.id)
    return apply(pr.pr, bits, id_bit)
end

function apply(
    op::TensoredOperator,
    state::State,
)
    states = Vector{State}([op.coeff * state])

    for pr in reverse(op.prs)
        states = apply(pr, states)
    end

    return states
end

function apply(
    op::TensoredOperator,
    states::Vector{State},
)
    return reduce(vcat, [apply(op, state) for state in states])
end

function apply(
    op::TensoredOperator,
    bits::Int,
)
    return apply(op, State(bits))
end

function reverse_bits(bits::Int, num_bits::Int)
    bits_new = 0
    for i in 0:(num_bits-1)
        bits_new |= ((bits >> i) & 1) << (num_bits - 1 - i)
    end
    return bits_new
end

function build_matrix(::Type{N}, ::Space{T}, op::AbstractOperator{T}) where {N<:Number,T<:AbstractSystemTag}
    throw(ArgumentError("Unsupported operator type: $(typeof(op))"))
end

function build_matrix(space::Space{T}, op::AbstractOperator{T}) where {T<:AbstractSystemTag}
    return build_matrix(ComplexF64, space, op)
end

function build_matrix(::Type{N}, space::Space{T}, op::TensoredOperator{T}) where {N<:Number,T<:AbstractSystemTag}
    num_bits = nindices(space)
    dim_system = dim(space)
    basis_sector = basis(space)
    mat = zeros(N, dim_system, dim_system)

    for ket in basis_sector
        ket_reversed = reverse_bits(ket, num_bits)
        states = apply(op, ket_reversed)

        for state in states
            bra = reverse_bits(state.bits, num_bits)

            if bra in basis_sector
                bra_sector = searchsortedfirst(basis_sector, bra)
                ket_sector = searchsortedfirst(basis_sector, ket)
                mat[bra_sector, ket_sector] += state.coeff
            end
        end
    end

    return mat
end

function build_matrix(::Type{N}, space::Space{T}, op::SummedOperator{T}) where {N<:Number,T<:AbstractSystemTag}
    dim_system = dim(space)
    mat = zeros(N, dim_system, dim_system)

    for op_i in op.ops
        mat += build_matrix(N, space, op_i)
    end

    return mat
end
