module Matrices

using LatticeSpaces
using LatticeOperators

include("primitive.jl")
include("operator.jl")

export apply, reverse_bits, build_matrix

end
