using LatticeBackends
using Test

using LatticeSpaces
using LatticeOperators
using LinearAlgebra
using SparseArrays

function st(bits::Integer, coeff::Number=1.0)
    return State(Int(bits), coeff)
end

include("matrices/primitive.jl")
include("matrices/operator.jl")
include("matrices/matrix.jl")
