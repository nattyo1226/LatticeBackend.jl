function test_matrix()
    geometry = Hypercubic(2, 2)
    boundary = OpenBoundary()
    lattice = Lattice(geometry, boundary)
    model = TFIHamiltonian(lattice, -1.0, -2.0)

    hamiltonian_matrix = build_matrix(lattice, model)
    hamiltonian_matrix_true = Float64[
        -4 -2 -2 0 -2 0 0 0 -2 0 0 0 0 0 0 0
        -2 0 0 -2 0 -2 0 0 0 -2 0 0 0 0 0 0
        -2 0 0 -2 0 0 -2 0 0 0 -2 0 0 0 0 0
        0 -2 -2 0 0 0 0 -2 0 0 0 -2 0 0 0 0
        -2 0 0 0 0 -2 -2 0 0 0 0 0 -2 0 0 0
        0 -2 0 0 -2 0 0 -2 0 0 0 0 0 -2 0 0
        0 0 -2 0 -2 0 4 -2 0 0 0 0 0 0 -2 0
        0 0 0 -2 0 -2 -2 0 0 0 0 0 0 0 0 -2
        -2 0 0 0 0 0 0 0 0 -2 -2 0 -2 0 0 0
        0 -2 0 0 0 0 0 0 -2 4 0 -2 0 -2 0 0
        0 0 -2 0 0 0 0 0 -2 0 0 -2 0 0 -2 0
        0 0 0 -2 0 0 0 0 0 -2 -2 0 0 0 0 -2
        0 0 0 0 -2 0 0 0 -2 0 0 0 0 -2 -2 0
        0 0 0 0 0 -2 0 0 0 -2 0 0 -2 0 0 -2
        0 0 0 0 0 0 -2 0 0 0 -2 0 -2 0 0 -2
        0 0 0 0 0 0 0 -2 0 0 0 -2 0 -2 -2 -4
    ]
    @test size(hamiltonian_matrix) == (16, 16)
    @test hamiltonian_matrix == hamiltonian_matrix_true
end

@testset "Matrix builder functions" begin
    test_matrix()
end
