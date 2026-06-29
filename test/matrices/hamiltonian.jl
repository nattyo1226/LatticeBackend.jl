function test_tfi()
    j = -1.0
    h = -2.0

    space = Space(
        SpinHalfSpace(),
        Hypercubic((2, 2), OpenBoundary),
    )
    model = TFIHamiltonian(space, j, h)

    hamiltonian_matrix = build_matrix(Float64, space, model)

    hamiltonian_matrix_true = zeros(Float64, 16, 16)

    for i in 1:16
        ps = [(i - 1) >> r & 1 for r in 0:3]
        ps = 2 .* ps .- 1
        hamiltonian_matrix_true[i, i] += j * (ps[1] * ps[2] + ps[1] * ps[3] + ps[2] * ps[4] + ps[3] * ps[4])
    end

    for i in 1:16
        for r in 0:3
            j = ((i - 1) ⊻ (1 << r)) + 1
            hamiltonian_matrix_true[i, j] += h
        end
    end

    @test size(hamiltonian_matrix) == (16, 16)
    @test hamiltonian_matrix == hamiltonian_matrix_true
end

function test_hubbard_1()
    t = 1.0
    u = 2.0

    space = Space(
        SpinfulFermionSpace(),
        Hypercubic(2, OpenBoundary),
    )
    model = HubbardHamiltonian(space, t, u)

    hamiltonian_matrix = Float64.(build_matrix(ComplexF64, space, model))

    hamiltonian_matrix_true = zeros(Float64, 16, 16)

    for i in 1:16
        if (((i - 1) >> 0) & 1 == 1) && (((i - 1) >> 1) & 1 == 1)
            hamiltonian_matrix_true[i, i] += u
        end
        if (((i - 1) >> 2) & 1 == 1) && (((i - 1) >> 3) & 1 == 1)
            hamiltonian_matrix_true[i, i] += u
        end
    end

    hamiltonian_matrix_true[2, 5] += t
    hamiltonian_matrix_true[5, 2] += t

    hamiltonian_matrix_true[3, 9] += t
    hamiltonian_matrix_true[9, 3] += t

    hamiltonian_matrix_true[4, 7] += -t
    hamiltonian_matrix_true[7, 4] += -t

    hamiltonian_matrix_true[4, 10] += t
    hamiltonian_matrix_true[10, 4] += t

    hamiltonian_matrix_true[7, 13] += -t
    hamiltonian_matrix_true[13, 7] += -t

    hamiltonian_matrix_true[8, 14] += -t
    hamiltonian_matrix_true[14, 8] += -t

    hamiltonian_matrix_true[10, 13] += t
    hamiltonian_matrix_true[13, 10] += t

    hamiltonian_matrix_true[12, 15] += -t
    hamiltonian_matrix_true[15, 12] += -t

    @test size(hamiltonian_matrix) == (16, 16)
    @test hamiltonian_matrix == hamiltonian_matrix_true
end

function test_hubbard_2()
    t = 1.0
    u = 2.0

    space = Space(
        SpinfulFermionSpace(),
        Hypercubic(2, OpenBoundary),
    )
    model = SymmetricHubbardHamiltonian(space, t, u)

    hamiltonian_matrix = Float64.(build_matrix(ComplexF64, space, model))

    hamiltonian_matrix_true = zeros(Float64, 16, 16)

    for i in 1:16
        if ((i - 1) >> 0) & 1 == ((i - 1) >> 1) & 1
            hamiltonian_matrix_true[i, i] += u / 4
        else
            hamiltonian_matrix_true[i, i] += -u / 4
        end

        if ((i - 1) >> 2) & 1 == ((i - 1) >> 3) & 1
            hamiltonian_matrix_true[i, i] += u / 4
        else
            hamiltonian_matrix_true[i, i] += -u / 4
        end
    end

    hamiltonian_matrix_true[2, 5] += t
    hamiltonian_matrix_true[5, 2] += t

    hamiltonian_matrix_true[3, 9] += t
    hamiltonian_matrix_true[9, 3] += t

    hamiltonian_matrix_true[4, 7] += -t
    hamiltonian_matrix_true[7, 4] += -t

    hamiltonian_matrix_true[4, 10] += t
    hamiltonian_matrix_true[10, 4] += t

    hamiltonian_matrix_true[7, 13] += -t
    hamiltonian_matrix_true[13, 7] += -t

    hamiltonian_matrix_true[8, 14] += -t
    hamiltonian_matrix_true[14, 8] += -t

    hamiltonian_matrix_true[10, 13] += t
    hamiltonian_matrix_true[13, 10] += t

    hamiltonian_matrix_true[12, 15] += -t
    hamiltonian_matrix_true[15, 12] += -t

    basis_sector = basis(space)
    hamiltonian_matrix_true = hamiltonian_matrix_true[basis_sector .+ 1, basis_sector .+ 1]

    @test size(hamiltonian_matrix) == (16, 16)
    @test hamiltonian_matrix == hamiltonian_matrix_true
end

@testset "Matrix builder functions" begin
    test_tfi()
    test_hubbard_1()
    test_hubbard_2()
end
