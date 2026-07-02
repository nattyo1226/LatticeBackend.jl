function apply!(out::Vector{State}, ::Int, ::Identity, state::State)
    push!(out, state)
    return out
end

function apply!(out::Vector{State}, id::Int, ::PauliX, state::State)
    mask = 1 << id
    bits_new = state.bits ⊻ mask

    push!(out, State(bits_new, state.coeff))
    return out
end

function apply!(out::Vector{State}, id::Int, ::PauliY, state::State)
    mask = 1 << id
    bits_new = state.bits ⊻ mask

    coeff_new = iszero((state.bits >> id) & 1) ? state.coeff * im : -state.coeff * im

    push!(out, State(bits_new, coeff_new))
    return out
end

function apply!(out::Vector{State}, id::Int, ::PauliZ, state::State)
    c = iszero((state.bits >> id) & 1) ? 1.0 : -1.0
    push!(out, c * state)
    return out
end

function apply!(out::Vector{State}, id::Int, ::MajoranaX, state::State)
    mask = 1 << id
    bits_new = state.bits ⊻ mask

    parity = isodd(count_ones(state.bits & (mask - 1))) ? -1.0 : 1.0
    coeff_new = parity * state.coeff

    push!(out, State(bits_new, coeff_new))
    return out
end

function apply!(out::Vector{State}, id::Int, ::MajoranaY, state::State)
    mask = 1 << id
    bits_new = state.bits ⊻ mask

    parity = isodd(count_ones(state.bits & (mask - 1))) ? -1.0 : 1.0
    occupied = isone((state.bits >> id) & 1)
    coeff_new = occupied ? -im * parity * state.coeff : im * parity * state.coeff

    push!(out, State(bits_new, coeff_new))
    return out
end
