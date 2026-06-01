function build_term_type(term::LocalOnsiteTerm)
    return eltype(op_to_matrix(term.op))
end

function build_term_type(term::OnsiteTerm)
    return eltype(op_to_matrix(term.op))
end

function build_term_type(term::LocalPairTerm)
    T1 = eltype(op_to_matrix(term.op1))
    T2 = eltype(op_to_matrix(term.op2))
    return promote_type(T1, T2)
end

function build_term_type(term::PairTerm)
    T1 = eltype(op_to_matrix(term.op1))
    T2 = eltype(op_to_matrix(term.op2))
    return promote_type(T1, T2)
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::LocalOnsiteTerm) where T<:Number
    num_sites = Int(nsites(lattice))
    matrix_tmp = op_to_matrix(T, term.op)
    return build_local_matrix(T, num_sites, term.id_site, matrix_tmp)
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::OnsiteTerm) where T<:Number
    num_sites = Int(nsites(lattice))
    matrix_result = build_zero_matrix(T, num_sites)

    for id in 1:nsites(lattice)
        matrix_tmp = op_to_matrix(T, term.op)
        matrix_result += build_local_matrix(T, num_sites, id, matrix_tmp)
    end

    return matrix_result
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::LocalPairTerm) where T<:Number
    num_sites = Int(nsites(lattice))

    matrix_tmp1 = op_to_matrix(T, term.op1)
    matrix_tmp2 = op_to_matrix(T, term.op2)
    return build_local_matrix(T, num_sites, (term.id_site1, term.id_site2), (matrix_tmp1, matrix_tmp2))
end

function build_term_matrix(::Type{T}, lattice::Lattice, term::PairTerm) where T<:Number
    num_sites = Int(nsites(lattice))

    matrix_result = build_zero_matrix(T, num_sites)
    for (id1, id2) in neighbor_pairs(lattice, term.shell)
        matrix_tmp1 = op_to_matrix(T, term.op1)
        matrix_tmp2 = op_to_matrix(T, term.op2)
        matrix_result += build_local_matrix(T, num_sites, (id1, id2), (matrix_tmp1, matrix_tmp2))
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
