function test_matrix_1()
    space = Space(SpinHalfSpace(), Hypercubic(1, OpenBoundary))

    x = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    y = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliY())
    z = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliZ())

    X_expected = ComplexF64[
        0 1
        1 0
    ]

    Y_expected = ComplexF64[
        0 -im
        im 0
    ]

    Z_expected = ComplexF64[
        1 0
        0 -1
    ]

    @test dense_matrix(space, x) == X_expected
    @test dense_matrix(space, y) == Y_expected
    @test dense_matrix(space, z) == Z_expected

    @test Matrix(sparse_matrix(space, x)) == X_expected
    @test Matrix(sparse_matrix(space, y)) == Y_expected
    @test Matrix(sparse_matrix(space, z)) == Z_expected
end

function test_matrix_2()
    space = Space(SpinHalfSpace(), Hypercubic(2, OpenBoundary))

    x1 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    z2 = LocalOperator(SiteIndex{SpinHalfTag}(2), PauliZ())

    op = x1 * z2 + 2.0 * z2

    dense = dense_matrix(space, op)
    sparse = sparse_matrix(space, op)

    @test Matrix(sparse) == dense
end

function test_matrix_3()
    space = Space(SpinHalfSpace(), Hypercubic(2, OpenBoundary))

    x1 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    z2 = LocalOperator(SiteIndex{SpinHalfTag}(2), PauliZ())

    op = x1 * z2

    mat = dense_matrix(space, op)

    # X1 Z2 |00⟩ =  |01⟩
    # X1 Z2 |01⟩ =  |00⟩
    # X1 Z2 |10⟩ = -|11⟩
    # X1 Z2 |11⟩ = -|10⟩

    expected = ComplexF64[
        0 0 1 0
        0 0 0 -1
        1 0 0 0
        0 -1 0 0
    ]

    @test mat == expected
    @test Matrix(sparse_matrix(space, op)) == expected
end

function test_matrix_4()
    space = Space(SpinHalfSpace(), Hypercubic(1, OpenBoundary))

    x = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())

    op = x + 2.0 * x

    expected = ComplexF64[
        0 3
        3 0
    ]

    @test dense_matrix(space, op) == expected
    @test Matrix(sparse_matrix(space, op)) == expected
end

function test_matrix_5()
    space = Space(SpinHalfSpace(), Hypercubic(1, OpenBoundary))

    x = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())

    op = x - x

    expected = zeros(ComplexF64, 2, 2)

    @test dense_matrix(space, op) == expected
    @test Matrix(sparse_matrix(space, op)) == expected
end

function test_tfi_1()
    space = Space(SpinHalfSpace(), Hypercubic(2, OpenBoundary))

    j = 1.0
    h = 0.5
    H = tfi(space, j, h)

    Hd = dense_matrix(space, H)
    Hs = sparse_matrix(space, H)

    @test size(Hd) == (4, 4)
    @test Matrix(Hs) == Hd
    @test Hd ≈ Hd'

    expected = ComplexF64[
        -j -h -h 0
        -h j 0 -h
        -h 0 j -h
        0 -h -h -j
    ]

    @test Hd ≈ expected
end

function test_tfi_2()
    space = Space(SpinHalfSpace(), Hypercubic(3, OpenBoundary))

    H = tfi(space, 1.2, 0.7)

    Hd = dense_matrix(space, H)
    Hs = sparse_matrix(space, H)

    @test size(Hd) == (8, 8)
    @test Matrix(Hs) ≈ Hd
    @test Hd ≈ Hd'
end

function test_tfi_3()
    space = Space(SpinHalfSpace(), Hypercubic(3, PeriodicBoundary))

    H = tfi(space, 1.0, 0.3)

    Hd = dense_matrix(space, H)
    Hs = sparse_matrix(space, H)

    @test size(Hd) == (8, 8)
    @test Matrix(Hs) ≈ Hd
    @test Hd ≈ Hd'
end

function test_hubbard_1()
    space = Space(SpinfulFermionSpace(), Hypercubic(2, OpenBoundary))

    t = 1.0
    u = 2.0
    H = hubbard(space, t, u)

    Hd = dense_matrix(space, H)
    Hs = sparse_matrix(space, H)

    @test size(Hd) == (16, 16)
    @test Matrix(Hs) ≈ Hd
    @test Hd ≈ Hd'

    expected = ComplexF64[
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 -t 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 -t 0 0 0 0 0 0 0
        0 0 0 u 0 0 t 0 0 -t 0 0 0 0 0 0
        0 -t 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 t 0 0 0 0 0 0 0 0 t 0 0 0
        0 0 0 0 0 0 0 u 0 0 0 0 0 t 0 0
        0 0 -t 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 -t 0 0 0 0 0 0 0 0 -t 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 u 0 0 1 0
        0 0 0 0 0 0 t 0 0 -t 0 0 u 0 0 0
        0 0 0 0 0 0 0 t 0 0 0 0 0 u 0 0
        0 0 0 0 0 0 0 0 0 0 0 t 0 0 u 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2u
    ]

    @test Hd ≈ expected
end

function test_hubbard_2()
    space = Space(
        SpinfulFermionSpace(),
        Hypercubic(2, OpenBoundary),
        ParticleNumberSector(2),
    )

    t = 1.0
    u = 2.0
    H = hubbard(space, t, u)

    Hd = dense_matrix(space, H)
    Hs = sparse_matrix(space, H)

    @test size(Hd) == (6, 6)
    @test Matrix(Hs) ≈ Hd
    @test Hd ≈ Hd'
end

function test_symmetric_hubbard_1()
    space = Space(
        SpinfulFermionSpace(),
        Hypercubic(2, OpenBoundary),
        ParticleNumberSector(2),
    )

    t = 1.0
    u = 2.0
    H = symmetric_hubbard(space, t, u)

    Hd = dense_matrix(space, H)
    Hs = sparse_matrix(space, H)

    @test size(Hd) == (6, 6)
    @test Matrix(Hs) ≈ Hd
    @test Hd ≈ Hd'
end

function test_symmetric_hubbard_2()
    space = Space(
        SpinfulFermionSpace(),
        Hypercubic(2, OpenBoundary),
        ParticleNumberSector(2),
    )

    t = 1.0
    u = 2.0

    H1 = dense_matrix(space, hubbard(space, t, u))
    H2 = dense_matrix(space, symmetric_hubbard(space, t, u))

    D = H1 - H2

    @test D ≈ Diagonal(diag(D))
    @test all(isapprox.(diag(D), diag(D)[1]))
end

@testset "Matrix builder" begin
    test_matrix_1()
    test_matrix_2()
    test_matrix_3()
    test_matrix_4()
    test_matrix_5()
end

@testset "Build Hamiltonians" begin
    test_tfi_1()
    test_tfi_2()
    test_tfi_3()
    test_hubbard_1()
    test_hubbard_2()
    test_symmetric_hubbard_1()
    test_symmetric_hubbard_2()
end
