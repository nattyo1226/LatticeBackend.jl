struct State
    bits::Int
    coeff::ComplexF64

    function State(bits::Int, coeff::Number=1.0)
        return new(bits, ComplexF64(coeff))
    end
end

function Base.:(==)(s1::State, s2::State)
    return s1.bits == s2.bits && s1.coeff == s2.coeff
end

function Base.:(*)(c::Number, s::State)
    return State(s.bits, c * s.coeff)
end

function apply!(::Vector{State}, ::Int, ap::AbstractPrimitive, ::State)
    throw(ArgumentError("Unsupported operator type: $(typeof(ap))"))
end

function apply!(out::Vector{State}, id::Int, ap::AbstractPrimitive, states::Vector{State})
    for state in states
        apply!(out, id, ap, state)
    end
    return out
end

function apply!(out::Vector{State}, id::Int, ap::AbstractPrimitive, bits::Int)
    return apply!(out, id, ap, State(bits))
end

function apply(id::Int, ap::AbstractPrimitive, state::State)
    out = Vector{State}()
    apply!(out, id, ap, state)
    return out
end

function apply!(::Vector{State}, ao::AbstractOperator, ::State)
    throw(ArgumentError("Unsupported operator type: $(typeof(ao))"))
end

function apply!(out::Vector{State}, ao::AbstractOperator, states::Vector{State})
    for state in states
        apply!(out, ao, state)
    end
    return out
end

function apply!(out::Vector{State}, ao::AbstractOperator, bits::Int)
    return apply!(out, ao, State(bits))
end

function apply(ao::AbstractOperator, state::State)
    out = Vector{State}()
    apply!(out, ao, state)
    return out
end

function (ao::AbstractOperator)(state::State)
    return apply(ao, state)
end

function reverse_bits(bits::Int, num_bits::Int)
    bits_new = 0
    for i in 0:(num_bits-1)
        bits_new |= ((bits >> i) & 1) << (num_bits - 1 - i)
    end
    return bits_new
end

function sparse_matrix(space::Space{T}, ao::AbstractOperator{T}) where {T<:AbstractSystemTag}
    num_bits = nindices(space)
    dim_space = dim(space)
    basis_space = basis(space)
    state_to_col = Dict{Int,Int}(
        s => i for (i, s) in enumerate(basis_space)
    )

    out = Vector{State}()
    rows = Vector{Int}()
    cols = Vector{Int}()
    vals = Vector{ComplexF64}()

    for ket in basis_space
        empty!(out)

        ket_reversed = reverse_bits(ket, num_bits)
        apply!(out, ao, ket_reversed)

        for state in out
            bra = reverse_bits(state.bits, num_bits)

            if bra in basis_space
                bra_sector = state_to_col[bra]
                ket_sector = state_to_col[ket]

                push!(rows, bra_sector)
                push!(cols, ket_sector)
                push!(vals, state.coeff)
            end
        end
    end

    return sparse(rows, cols, vals, dim_space, dim_space)
end

function dense_matrix(space::Space{T}, ao::AbstractOperator{T}) where {T<:AbstractSystemTag}
    num_bits = nindices(space)
    dim_space = dim(space)
    basis_space = basis(space)
    state_to_col = Dict{Int,Int}(
        s => i for (i, s) in enumerate(basis_space)
    )

    out = Vector{State}()
    mat = zeros(ComplexF64, dim_space, dim_space)

    for ket in basis_space
        empty!(out)

        ket_reversed = reverse_bits(ket, num_bits)
        apply!(out, ao, ket_reversed)

        for state in out
            bra = reverse_bits(state.bits, num_bits)

            if bra in basis_space
                bra_sector = state_to_col[bra]
                ket_sector = state_to_col[ket]

                mat[bra_sector, ket_sector] += state.coeff
            end
        end
    end

    return mat
end
