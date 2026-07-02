function test_apply_primitive_1()
    s0 = st(0b00)
    s1 = st(0b01)
    s2 = st(0b10)

    @test apply(0, Identity{SpinHalfTag}(), s2) == [s2]

    @test apply(0, PauliX(), s0) == [st(0b01)]
    @test apply(0, PauliX(), s1) == [st(0b00)]
    @test apply(1, PauliX(), s0) == [st(0b10)]

    @test apply(0, PauliY(), s0) == [st(0b01, im)]
    @test apply(0, PauliY(), s1) == [st(0b00, -im)]
    @test apply(1, PauliY(), s0) == [st(0b10, im)]
    @test apply(1, PauliY(), s2) == [st(0b00, -im)]

    @test apply(0, PauliZ(), s0) == [st(0b00, 1)]
    @test apply(0, PauliZ(), s1) == [st(0b01, -1)]
    @test apply(0, PauliZ(), s2) == [st(0b10, 1)]
    @test apply(1, PauliZ(), s2) == [st(0b10, -1)]
end

function test_apply_primitive_2()
    @test apply(0, MajoranaX(), st(0b00)) == [st(0b01, 1)]
    @test apply(0, MajoranaX(), st(0b01)) == [st(0b00, 1)]

    @test apply(0, MajoranaY(), st(0b00)) == [st(0b01, im)]
    @test apply(0, MajoranaY(), st(0b01)) == [st(0b00, -im)]

    @test apply(1, MajoranaX(), st(0b00)) == [st(0b10, 1)]
    @test apply(1, MajoranaX(), st(0b01)) == [st(0b11, -1)]

    @test apply(1, MajoranaY(), st(0b00)) == [st(0b10, im)]
    @test apply(1, MajoranaY(), st(0b01)) == [st(0b11, -im)]
    @test apply(1, MajoranaY(), st(0b10)) == [st(0b00, -im)]
    @test apply(1, MajoranaY(), st(0b11)) == [st(0b01, im)]
end

function test_apply_primitive_3()
    out = State[st(0b11, 2.0)]

    apply!(out, 0, PauliX(), st(0b00))

    @test out == [
        st(0b11, 2.0),
        st(0b01, 1.0),
    ]
end

function test_apply_primitive_4()
    states = [
        st(0b00, 2.0),
        st(0b01, 3.0),
    ]

    out = State[]
    apply!(out, 0, PauliX(), states)

    @test out == [
        st(0b01, 2.0),
        st(0b00, 3.0),
    ]
end

@testset "Primitives application" begin
    test_apply_primitive_1()
    test_apply_primitive_2()
    test_apply_primitive_3()
    test_apply_primitive_4()
end
