struct State
    bits::Int
    coeff::ComplexF64

    function State(bits::Int, coeff::Number=1.0)
        return new(bits, ComplexF64(coeff))
    end
end

function Base.:(*)(c::Number, s::State)
    return State(s.bits, c * s.coeff)
end
