module Matrices

using LatticeSpaces
using LatticeOperators

include("primitive.jl")
include("operator.jl")

export apply, build_matrix

end
