function apply(::IndexedOperatorPrimitive, ::Nothing)
    return nothing, 0.0
end

function apply(pr::IndexedOperatorPrimitive{T}, state::Int) where {T<:AbstractSystemTag}
    id_bit = to_bit(pr.id)
    return apply(pr.pr, state, id_bit)
end

function apply(
    op::TensoredOperator,
    state::Int,
)
    state_new = state
    coeff_new = op.coeff

    for pr in reverse(op.prs)
        state_new, coeff_tmp = apply(pr, state_new)
        if state_new === nothing
            return nothing, 0.0
        end
        coeff_new *= coeff_tmp
    end

    return state_new, coeff_new
end

function apply(
    ::TensoredOperator,
    ::Nothing
)
    return nothing, 0.0
end

function reverse_bits(state::Int, num_bits::Int)
    state_new = 0
    for i in 0:(num_bits-1)
        state_new |= ((state >> i) & 1) << (num_bits - 1 - i)
    end
    return state_new
end

function reverse_bits(::Nothing, ::Int)
    return nothing
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
    mat = zeros(N, dim_system, dim_system)

    for ket in 0:(dim_system-1)
        ket_reversed = reverse_bits(ket, num_bits)
        bra_reversed, coeff = apply(op, ket_reversed)

        if bra_reversed === nothing
            continue
        end

        bra = reverse_bits(bra_reversed, num_bits)
        mat[bra+1, ket+1] += coeff
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
