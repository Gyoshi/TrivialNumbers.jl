using Symbolics

include("trivial.jl")

@variables xa, xb, xc, ya, yb
product = Trivial(xa, xb) * Trivial(ya, yb);
simplify(product.a; expand=true)
simplify(product.b; expand=true)

Symbolics.Num(x::Trivial) = isreal(x) ? x.a : throw(InexactError(:Num, Trivial, x))
division = Trivial(xa, xb) / Trivial(ya, yb);
simplify(division.a; expand=true)
simplify(division.b; expand=true)

conjugate = conj(Trivial(xa, xb));
norm2 = Trivial(xa, xb)*conjugate;
simplify(norm2.a; expand=true)

conjugate = conj(Trivial(xa, xb, xc));
norm2 = Trivial(xa, xb, xc)*conjugate;
simplify(norm2.a; expand=true)

include("quadrivial.jl")
Symbolics.Num(x::Quadrivial) = isreal(x) ? x.a : throw(InexactError(:Num, Quadrivial, x))

@variables a, b, c
x = Quadrivial(a, b, c);
product = abs3(x);
simplify(product; expand=true)

real(x).a

using LinearAlgebra
@variables a, b, c
quad = reshape([a a-c a-b b-c b b-a c-b c-a c], (3,3))
det(quad) |> expand

# V4 surfaces
@variables a, b, c, n, ϵ
-6a ~ -2(b+c) + sqrt

ex = -3a^2 + 2(b+c)*a + (b-c)^2
ex = -3b^2 + 2(c+a)*b + (c-a)^2
ex = -3c^2 + 2(a+b)*c + (a-b)^2
substitute(ex |> expand, [a^2 => 4n^2 - 4n*ϵ, a => 2n-ϵ, b => n, c => n]) |> expand

# V_C5
# @variables x[0:4], r

# id1 = r^5 => 1
# id2 = sum(r^i for i = 0:4) => 0
# rules = Dict(id1, id2)

# # ex = 2*sum(r^i for i = 0:4)
# # ex = (1+r^2)*(r^4)
# ex = prod(sum(x[i]*r^e for (i, e) in zip(0:4, g)) for g in [[0, 1, 2, 3, 4], [0, 3, 1, 4, 2], [0, 4, 3, 2, 1], [0, 2, 4, 1, 3]])
# ex = simplify(ex, expand=true)
# # substitute(ex, rules) # doesn't work =<

# cycle = Dict(r^i => r^mod(i, 5) for i in 5:16)
# ex2 = substitute(ex, cycle)
# ex2 = substitute(ex2, r + r^2 + r^3 + r^4 => -1)
# ex2 = simplify(ex2, expand=true)
# ex2 = substitute(ex2, r + r^2 + r^3 + r^4 => -1)

# sumrules = Dict{Num, Num}()
# for a in 0:1, b in 0:1, c in 0:1, d in 0:1, e in 0:1
#     if a + b + c + d + e == 3
#         push!(sumrules, a + b*r + c*r^2 + d*r^3 + e*r^4 => (a-1) + (b-1)*r + (c-1)*r^2 + (d-1)*r^3 + (e-1)*r^4)
#     end
# end
# ex2 = substitute(ex2, sumrules)
# ex2 = simplify(ex2, expand=true)
# ex2 = substitute(ex2, 2 + r + r^2 + r^3 + r^4 => 1)
# ex2 = substitute(ex2, 2 + r^3 + r^4 => 1 - r - r^2)
# ex2 = substitute(ex2, (2 + r^2 + r^3 + r^4) => 1 - r)
# ex2 = simplify(ex2, expand=true)

# ex3 = include("prod_C5-canceled.txt")
# ex3 = simplify(ex3, expand=true)

