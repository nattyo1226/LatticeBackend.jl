module Matrices

using LatticeSpaces
using LatticeOperators
using SparseArrays

include("type.jl")
export State, apply, apply!, dense_matrix, sparse_matrix

include("primitive.jl")

include("operator.jl")

end
