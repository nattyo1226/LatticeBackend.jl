function test_apply_operator_1()
    lo_x = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    lo_z = LocalOperator(SiteIndex{SpinHalfTag}(2), PauliZ())

    @test lo_x(st(0b00)) == [st(0b01)]
    @test lo_x(st(0b01)) == [st(0b00)]

    @test lo_z(st(0b00)) == [st(0b00)]
    @test lo_z(st(0b10)) == [st(0b10, -1)]
end

function test_apply_operator_2()
    # X_1 Z_2 |00⟩ = |01⟩
    lo1 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    lo2 = LocalOperator(SiteIndex{SpinHalfTag}(2), PauliZ())

    po = lo1 * lo2
    @test po(st(0b00)) == [st(0b01)]
    @test po(st(0b10)) == [st(0b11, -1)]

    # coefficient
    po2 = 2.0im * po
    @test po2(st(0b10)) == [st(0b11, -2.0im)]
end

function test_apply_operator_3()
    lo1 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliZ())
    lo2 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())

    po = lo1 * lo2

    @test po(st(0b0)) == [st(0b1, -1)]
end

function test_apply_operator_4()
    lo1 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    lo2 = LocalOperator(SiteIndex{SpinHalfTag}(1), PauliZ())

    so = lo1 + lo2

    # (X + Z)|0⟩ = |1⟩ + |0⟩
    @test Set(so(st(0b0))) == Set([st(0b1), st(0b0)])

    # (X + Z)|1⟩ = |0⟩ - |1⟩
    @test Set(so(st(0b1))) == Set([st(0b0), st(0b1, -1)])
end

function test_apply_operator_5()
    po1 = 1.0 * LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    po2 = 2.0 * LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())

    so = po1 + po2

    # (X + 2X)|0⟩ = 3|1⟩
    @test so(st(0b0)) == [st(0b1, 3.0)]
end

function test_apply_operator_6()
    po1 = 1.0 * LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())
    po2 = -1.0 * LocalOperator(SiteIndex{SpinHalfTag}(1), PauliX())

    so = po1 + po2

    # (X - X)|0⟩ = 0
    @test so(st(0b0)) == State[]
end

@testset "Operators application" begin
    test_apply_operator_1()
    test_apply_operator_2()
    test_apply_operator_3()
    test_apply_operator_4()
    test_apply_operator_5()
    test_apply_operator_6()
end
