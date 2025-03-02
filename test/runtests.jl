using TrivialNumbers
using Test

@testset "TrivialNumbers.jl" begin
    ╱ = TrivialNumbers.╱
    ╲ = TrivialNumbers.╲
    x = (1+2╱)
    y = (2+3╱)
    z = Trivial(2,-2)
    @test TrivialNumbers.conj(1+2TrivialNumbers.╱) == 1+2TrivialNumbers.╲
    @test TrivialNumbers.verso(Trivial(0, 2, 3)) == Trivial(3, 0, 2)
    @test -(x,y,z) == -(TrivialNumbers.dual(x,y,z)...)
    @test exp(0.5+╲) |> abs ≈ 1
end
