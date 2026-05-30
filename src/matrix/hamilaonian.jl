function build_term_type(term::OnsiteTerm)
    return eltype(op_to_matrix(term.op))
end

function build_term_type(term::PairTerm)
    T1 = eltype(op_to_matrix(term.op1))
    T2 = eltype(op_to_matrix(term.op2))
    return promote_type(T1, T2)
end

function build_local_matrix(::Type{T}, num_sites::Int, id_site::Int, term::OnsiteTerm) where T<:Number
    matrix = term.coeff * op_to_matrix(term.op)
    return build_local_matrix(T, num_sites, id_site, matrix)
end

function build_local_matrix(::Type{T}, lattice::Lattice, id_site::Int, term::OnsiteTerm) where T<:Number
    num_sites = Int(nsites(lattice))
    return build_local_matrix(T, num_sites, id_site, term)
end

function build_local_matrix(num_sites::Int, id_site::Int, term::OnsiteTerm)
    T = build_term_type(term)
    return build_local_matrix(T, num_sites, id_site, term)
end

function build_local_matrix(lattice::Lattice, id_site::Int, term::OnsiteTerm)
    T = build_term_type(term)
    return build_local_matrix(T, lattice, id_site, term)
end

function build_local_matrix(::Type{T}, num_sites::Int, ids_site::NTuple{2,Int}, term::PairTerm) where T<:Number
    matrices = (term.coeff * op_to_matrix(term.op1), op_to_matrix(term.op2))
    return build_local_matrix(T, num_sites, ids_site, matrices)
end

function build_local_matrix(::Type{T}, lattice::Lattice, ids_site::NTuple{2,Int}, term::PairTerm) where T<:Number
    num_sites = Int(nsites(lattice))
    return build_local_matrix(T, num_sites, ids_site, term)
end

function build_local_matrix(num_sites::Int, ids_site::NTuple{2,Int}, term::PairTerm)
    T = build_term_type(term)
    return build_local_matrix(T, num_sites, ids_site, term)
end

function build_local_matrix(lattice::Lattice, ids_site::NTuple{2,Int}, term::PairTerm)
    T = build_term_type(term)
    return build_local_matrix(T, lattice, ids_site, term)
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::OnsiteTerm) where T<:Number
    num_sites = Int(nsites(lattice))
    matrix_result = build_zero_matrix(T, num_sites)

    for id in 1:nsites(lattice)
        matrix_result += build_local_matrix(T, num_sites, id, term)
    end

    return matrix_result
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::PairTerm) where T<:Number
    num_sites = Int(nsites(lattice))

    matrix_result = build_zero_matrix(T, num_sites)
    for (id1, id2) in neighbor_pairs(lattice, term.shell)
        matrix_result += build_local_matrix(T, num_sites, (id1, id2), term)
    end

    return matrix_result
end

function build_term_matrix(lattice::Lattice, term::AbstractTerm)
    T = build_term_type(term)
    return build_term_matrix(T, lattice, term)
end

function build_hamiltonian_matrix(lattice::Lattice, model::GenericModel)
    num_sites = nsites(lattice)
    T = promote_type(build_term_type.(model.terms)...)
    matrix_result = build_zero_matrix(T, num_sites)

    for term in model.terms
        matrix_result += build_term_matrix(T, lattice, term)
    end

    return matrix_result
end
