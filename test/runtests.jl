using TrivialNumbers
using Test

const x = (1+2╱)
const y = (2+3╱)
const z = Trivial(2,-2)

@testset "TrivialNumbers.jl" begin
    @test TrivialNumbers.conj(1+2TrivialNumbers.╱) == 1+2TrivialNumbers.╲
    @test TrivialNumbers.verso(Trivial(0, 2, 3)) == Trivial(3, 0, 2)
    @test -(x,y,z) == -(dual(x,y,z)...)
end
