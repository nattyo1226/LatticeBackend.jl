module LatticeOperator

using LatticeModel

include("matrix/Matrix.jl")
using .Matrix
export build_local_matrix, build_term_matrix, build_hamiltonian_matrix
end
