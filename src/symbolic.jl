include("trivial.jl")

using Symbolics
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

