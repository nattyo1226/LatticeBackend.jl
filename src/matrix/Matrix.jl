module Matrix

using LatticeGeometry
using LatticeModel
using SparseArrays

include("../utils.jl")

include("common.jl")
export build_local_matrix

include("operator.jl")
export build_term_matrix, build_hamiltonian_matrix

end
