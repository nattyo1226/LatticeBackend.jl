function apply!(
    out::Vector{State},
    lo::LocalOperator,
    state::State,
)
    id = to_bit(lo.id)
    return apply!(out, id, lo.pr, state)
end

function apply!(
    out::Vector{State},
    po::ProductOperator,
    state::State,
)
    states = Vector{State}([state])
    buffer = Vector{State}()

    for lo in Iterators.reverse(po.los)
        empty!(buffer)
        apply!(buffer, lo, states)

        states, buffer = buffer, states
        isempty(states) && break
    end

    for s in states
        push!(out, po.coeff * s)
    end

    return out
end

function combine_like_terms!(out::Vector{State}, states::Vector{State})
    sort!(states, by=s -> s.bits)

    i = firstindex(states)
    while i <= lastindex(states)
        bits = states[i].bits
        coeff = zero(ComplexF64)

        while i <= lastindex(states) && states[i].bits == bits
            coeff += states[i].coeff
            i += 1
        end

        if !iszero(coeff)
            push!(out, State(bits, coeff))
        end
    end

    return out
end

function apply!(
    out::Vector{State},
    so::SumOperator,
    state::State,
)
    buffer = Vector{State}()

    for po in so.pos
        apply!(buffer, po, state)
    end

    combine_like_terms!(out, buffer)

    return out
end
