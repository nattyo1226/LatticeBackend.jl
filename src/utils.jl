function identity_matrix(::Type{T}) where T<:Number
    return T[1 0; 0 1]
end

function pauli_x_matrix(::Type{T}) where T<:Number
    return T[0 1; 1 0]
end

function pauli_y_matrix(::Type{T}) where T<:Number
    return T[0 -1im; 1im 0]
end

function pauli_z_matrix(::Type{T}) where T<:Number
    return T[1 0; 0 -1]
end

function op_to_matrix(::Type{T}, operator::AbstractOperator) where T<:Number
    if operator isa Identity
        return identity_matrix(T)
    elseif operator isa PauliX
        return pauli_x_matrix(T)
    elseif operator isa PauliY
        return pauli_y_matrix(T)
    elseif operator isa PauliZ
        return pauli_z_matrix(T)
    else
        error("Unsupported operator type: $(typeof(operator))")
    end
end

function op_to_matrix(operator::AbstractOperator)
    return op_to_matrix(Float64, operator)
end
