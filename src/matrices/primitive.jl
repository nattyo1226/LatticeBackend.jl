function apply(pr::AbstractOperatorPrimitive, ::State, ::Int)
    throw(ArgumentError("Unsupported operator type: $(typeof(pr))"))
end

function apply(pr::AbstractOperatorPrimitive, states::Vector{State}, id::Int)
    return reduce(vcat, [apply(pr, state, id) for state in states]; init=Vector{State}())
end

function apply(pr::AbstractOperatorPrimitive, bits::Int, id::Int)
    return apply(pr, State(bits), id)
end

function apply(::Identity, state::State, ::Int)
    return [state]
end

function apply(::PauliX, state::State, id::Int)
    bits_new = state.bits ⊻ (1 << id)
    return [State(bits_new, state.coeff)]
end

function apply(::PauliY, state::State, id::Int)
    bits_new = state.bits ⊻ (1 << id)
    coeff_new = iszero((state.bits >> id) & 1) ? state.coeff * im : -state.coeff * im
    return [State(bits_new, coeff_new)]
end

function apply(::PauliZ, state::State, id::Int)
    c = iszero((state.bits >> id) & 1) ? 1.0 : -1.0
    return [c * state]
end

function apply(::Creation, state::State, id::Int)
    occ_id = (state.bits >> id) & 1
    if isone(occ_id)
        return Vector{State}()
    end

    state_new = state.bits | (1 << id)
    coeff_new = begin
        parity = count_ones(state.bits & ((1 << id) - 1))
        isodd(parity) ? -state.coeff : state.coeff
    end

    return [State(state_new, coeff_new)]
end

function apply(::Annihilation, state::State, id::Int)
    occ_id = (state.bits >> id) & 1
    if iszero(occ_id)
        return Vector{State}()
    end

    state_new = state.bits ⊻ (1 << id)
    coeff_new = begin
        parity = count_ones(state.bits & ((1 << id) - 1))
        isodd(parity) ? -state.coeff : state.coeff
    end

    return [State(state_new, coeff_new)]
end

function apply(::Occupation, state::State, id::Int)
    occ_id = (state.bits >> id) & 1

    if iszero(occ_id)
        return Vector{State}()
    else
        return [state]
    end
end

function apply(::MajoranaX, state::State, id::Int)
    return [apply(Creation(), state, id); apply(Annihilation(), state, id)]
end

function apply(::MajoranaY, state::State, id::Int)
    return [im * apply(Creation(), state, id); (-im) * apply(Annihilation(), state, id)]
end

function apply(::MajoranaZ, state::State, id::Int)
    return [state; (-2) * apply(Occupation(), state, id)]
end

function apply(pr::SummedOperatorPrimitive, state::State, id::Int)
    return reduce(vcat, [apply(pr, state, id) for pr in pr.prs]; init=Vector{State}())
end

function apply(pr::ProductedOperatorPrimitive, state::State, id::Int)
    states = Vector{State}([pr.coeff * state])
    for pr_tmp in reverse(pr.prs)
        states = apply(pr_tmp, states, id)
    end
    return states
end
