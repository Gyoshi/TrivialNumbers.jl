using StaticArrays

struct Quadrivial{T<:Real} <: Number
    triplet :: MVector{3,T}
    diffs :: MVector{3,T}
end

# Constructors
function Quadrivial(triplet::T1, diffs::T2) where {T1<:Vector, T2<:Vector}
    if !all(triplet-circshift(triplet, -1) .≈ diffs)
        throw(DomainError((triplet=triplet, diffs=diffs), "Difference array must be equal to the circular-sequential differences of the triplet."))
    end
    return Quadrivial{promote_type(T1.parameters[1], T2.parameters[1])}(triplet, diffs)
end
Quadrivial(x::Array) = Quadrivial(x, x-circshift(x, -1))
Quadrivial(a::Real, b::Real, c::Real) = Quadrivial([a,b,c])
Quadrivial(a::Real, b::Real, c::Real, d::Real) = Quadrivial(a-b, a-c, a-d)
Quadrivial(x::Real) = Quadrivial(x, x, x)

# Promotion
Base.convert(::Type{Quadrivial{T}}, z::Quadrivial{T}) where {T<:Real} = z
Base.convert(::Type{Quadrivial{T}}, z::Quadrivial) where {T<:Real} = Quadrivial{T}(convert(T, triplet(z)[1]), convert(T, z.b), convert(T, z.c))
Base.convert(::Type{Quadrivial{T}}, x::Real) where {T<:Real} = Quadrivial{T}(convert(T, x), convert(T, x), convert(T, x))
(::Type{T})(z::Quadrivial) where {T<:Real} = (isreal(z) || isnan(z) ? triplet(z)[1] : throw(InexactError(:convert, T, z)))
Base.promote_rule(::Type{Quadrivial{T}}, ::Type{S}) where {T<:Real, S<:Real} = Quadrivial{promote_type(T,S)}
Base.promote_rule(::Type{Quadrivial{T}}, ::Type{Quadrivial{S}}) where {T<:Real, S<:Real} = Quadrivial{promote_type(T,S)}

# Util
Base.getproperty(x::Quadrivial, _::Symbol) = error("Quadrivial struct fields are private to avoid passing of mutables. Use `triplet()` instead.")
triplet(x::Quadrivial) = Tuple(getfield(x, :triplet))
function quadruplet(x::Quadrivial)
    result = (0, -triplet(x)[1], -triplet(x)[2], -triplet(x)[3])
    return result .- min(result...)
end
quadruplet(x::Number) = quadruplet(Quadrivial(x))

verso(x::Quadrivial) = Quadrivial(triplet(x)[SVector(2,3,1)])
recto(x::Quadrivial) = Quadrivial(triplet(x)[SVector(3,1,2)])
verso(x::Number) = verso(Quadrivial(x))
recto(x::Number) = recto(Quadrivial(x))
verso(1)

x = Quadrivial(1,2,3)

Base.isreal(x::Quadrivial) = x.a ≈ x.b ≈ x.c 
Base.real(x::Quadrivial{T}) where T = convert(T, (x + verso(x) + recto(x)).triplet[1])/3
Base.imag(x::Quadrivial) = x - real(x)

Base.isnan(z::Quadrivial) = isnan.(triplet(z)) |> any
Base.isinf(z::Quadrivial) = isinf.(triplet(z)) |> any
Base.isfinite(z::Quadrivial) = isfinite.(triplet(z)) |> all

Base.:≈(a::Quadrivial, b::Quadrivial) = all(triplet(a) .≈ triplet(b))

# Operators
Base.:+(x::Quadrivial, y::Quadrivial) = Quadrivial(x.a + y.a, x.b + y.b, x.c + y.c)
Base.:-(x::Quadrivial) = Quadrivial(-x.a, -x.b, -x.c)
Base.:-(x::Quadrivial, y::Quadrivial) = x + -y

const i = Quadrivial(-1,0,0)
const j = Quadrivial(0,-1,0)
const k = Quadrivial(0,0,-1)

Base.:*(x::Quadrivial, y::Quadrivial) = Quadrivial(
    x.a*y.a + (x.c-x.b)*(y.c-y.b),
    x.b*y.b + (x.c-x.a)*(y.c-y.a),
    x.c*y.c + (x.a-x.b)*(y.a-y.b),
)

inner_product(x, y, z) = x*recto(y)*verso(z)
abs3(x::Quadrivial) = (x*verso(x)*recto(x)).a
abs3(x::Number) = abs3(Quadrivial(x))
Base.abs(x::Quadrivial) = cbrt(abs3(x))
Base.inv(x::Quadrivial) = Quadrivial(triplet(verso(x)*recto(x))./abs3(x)...)
Base.:/(x::Quadrivial, y::Quadrivial) = x*inv(y)

# Math
function Base.exp(x::Quadrivial)
    term = one(Quadrivial)
    n = 0
    rsum = 0.

    for i in 1:1000
        n += 1
        rsum += term
        term *= x/n        
    end
    return rsum
end

function Base.cos(x::Quadrivial)
    x2 = x*x
    term = one(Quadrivial)
    n = 0
    rsum = 0.

    for i in 1:1000
        n += 2
        rsum += term
        term *= -x2/(n*(n-1))      
    end
    return rsum
end
function Base.sin(x::Quadrivial)
    x2 = x*x
    term = x
    n = 1
    rsum = 0.

    for i in 1:1000
        n += 2
        rsum += term
        term *= -x2/(n*(n-1))      
    end
    return rsum
end

# Output
function quadrivial_show(io::IO, x::Quadrivial{T}) where {T}
    first_plus = true
    if all(triplet(x) .== zero(T))
        print(io, 0)
        return
    end
    coefficients = quadruplet(x)
    reps, unique_index = [count(==(element),coefficients) for element in unique(coefficients)]|> findmax
    if sum(coefficients.>0) > 2 && reps > 1
        coefficients = coefficients .- unique(coefficients)[unique_index]
    end
    for (nat, symb) in zip(coefficients, ["", "i", "j", "k"])
        if nat == zero(T)
            continue
        end
        if !first_plus
            if signbit(nat)==1 && !isnan(nat)
                nat = -nat
                print(io, " - ")
            else
                print(io, " + ")
            end
        else
            if signbit(nat)==1 && !isnan(nat)
                nat = -nat
                print(io, "-")
            end
            first_plus = false
        end
        if symb == "" || nat != one(T)
            show(io, nat)
        end
        print(io, symb)
    end
end
Base.show(io::IO, x::Quadrivial) = quadrivial_show(io, x)

