using TrivialNumbers
using Test

@testset "TrivialNumbers.jl" begin
    @test TrivialNumbers.conj(1+2TrivialNumbers.╱) == 1+2TrivialNumbers.╲
    @test TrivialNumbers.verso(Trivial(0, 2, 3)) == Trivial(3, 0, 2)
end
