function apply(pr::AbstractOperatorPrimitive, ::Int, ::Int)
    throw(ArgumentError("Unsupported operator type: $(typeof(pr))"))
end

function apply(::Identity, state::Int, ::Int)
    return state, 1.0
end

function apply(::PauliX, state::Int, id::Int)
    state_new = state ⊻ (1 << id)
    return state_new, 1.0
end

function apply(::PauliY, state::Int, id::Int)
    state_new = state ⊻ (1 << id)
    coeff_new = ((state >> id) & 1) == 0 ? 1.0 : -1.0
    return state_new, coeff_new * im
end

function apply(::PauliZ, state::Int, id::Int)
    coeff_new = ((state >> id) & 1) == 0 ? 1.0 : -1.0
    return state, coeff_new
end

function apply(::Creation, state::Int, id::Int)
    occ_id = (state >> id) & 1
    if occ_id == 1
        return nothing, 0.0
    end

    state_new = state | (1 << id)
    coeff_new = begin
        parity = count_ones(state & ((1 << id) - 1))
        isodd(parity) ? -1.0 : 1.0
    end

    return state_new, coeff_new
end

function apply(::Annihilation, state::Int, id::Int)
    occ_id = (state >> id) & 1
    if occ_id == 0
        return nothing, 0.0
    end

    state_new = state ⊻ (1 << id)
    coeff_new = begin
        parity = count_ones(state & ((1 << id) - 1))
        isodd(parity) ? -1.0 : 1.0
    end

    return state_new, coeff_new
end

function apply(::Occupation, state::Int, id::Int)
    occ_id = (state >> id) & 1
    coeff_new = occ_id == 0 ? 0.0 : 1.0
    return state, coeff_new
end
