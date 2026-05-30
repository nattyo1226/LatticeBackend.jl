function build_term_type(term::OnsiteTerm)
    return eltype(op_to_matrix(term.op))
end

function build_term_type(term::PairTerm)
    T1 = eltype(op_to_matrix(term.op1))
    T2 = eltype(op_to_matrix(term.op2))
    return promote_type(T1, T2)
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::OnsiteTerm) where T<:Number
    matrix = op_to_matrix(term.op)

    num_sites = nsites(lattice)
    matrix_result = build_zero_matrix(T, num_sites)

    for site in 1:nsites(lattice)
        matrix_result += term.coeff * build_local_matrix(T, num_sites, site, matrix)
    end

    return matrix_result
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::PairTerm) where T<:Number
    matrix1 = op_to_matrix(term.op1)
    matrix2 = op_to_matrix(term.op2)

    num_sites = nsites(lattice)

    matrix_result = build_zero_matrix(T, num_sites)
    for (site1, site2) in neighbor_pairs(lattice, term.shell)
        matrix_result += term.coeff * build_local_matrix(T, num_sites, (site1, site2), (matrix1, matrix2))
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
