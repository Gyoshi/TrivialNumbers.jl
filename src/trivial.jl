struct Trivial{T<:Real} <: Number
    # triplet : Tuple{T, T, T}
    # duplet :: Tuple{T, T}
    a :: T
    b :: T
end

# Constructors
# Trivial(a::S, b::T) where {S<:Number,T<:Number} = Trivial(promote(a,b)...)
Trivial(a, b) = Trivial{promote_type(typeof(a), typeof(b))}(a,b)
Trivial(a, b, c) = Trivial(a-b, b-c)
Trivial(x::Real) = Trivial(x, zero(x))
Trivial(x::T) where T<:Real = Trivial{T}(x, zero(T))

# Promotion
Base.convert(::Type{Trivial{T}}, z::Trivial{T}) where {T<:Real} = z
Base.convert(::Type{Trivial{T}}, z::Trivial) where {T<:Real} = Trivial{T}(convert(T, z.a), convert(T, z.b))
Base.convert(::Type{Trivial{T}}, x::Real) where {T<:Real} = Trivial{T}(convert(T, x), convert(T, 0))
(::Type{T})(z::Trivial) where {T<:Real} = (isreal(z) ? z.a : throw(InexactError(:convert, T, z)))
Base.promote_rule(::Type{Trivial{T}}, ::Type{S}) where {T<:Real, S<:Real} = Trivial{promote_type(T,S)}
Base.promote_rule(::Type{Trivial{T}}, ::Type{Trivial{S}}) where {T<:Real, S<:Real} = Trivial{promote_type(T,S)}

# Util
duplet(x::Trivial) = (x.a, x.b)
function triplet(x::Trivial)
    triplet = (0, -x.a, -x.a-x.b)
    return triplet .- min(triplet...)
end
isreal(x::Trivial) = x.b == 0
real(x::Trivial) = x.a

# Operators
Base.:+(x::Trivial, y::Trivial) = Trivial(x.a + y.a, x.b + y.b)
Base.:-(x::Trivial) = Trivial(-x.a, -x.b)
Base.:-(x::Trivial, y::Trivial) = x + -y

verso(x::Trivial) = Trivial(-x.a-x.b, x.a) 
recto(x::Trivial) = Trivial(x.b, -x.a-x.b)
verso(x::Real) = verso(Trivial(x))
recto(x::Real) = recto(Trivial(x))
const ╱ = verso(1)
const ╲ = recto(1)

Base.:*(x::Trivial, y::Trivial) = Trivial(x.a*y.a-x.b*y.b, (x.a + x.b)*y.b + x.b*y.a)
# Base.:*(x::Trivial, y::Trivial) = Trivial(reduce(.+, n.*triplet(op(y)) for (n, op) in zip(triplet(x), [x -> x, verso, recto]))...)
conj(x::Trivial) = Trivial(x.a+x.b, -x.b)
Base.abs2(x::Trivial{T}) where {T} = convert(T, x*conj(x))
Base.abs(x::Trivial) = sqrt(abs2(x))
Base.inv(x::Trivial) = Trivial(duplet(conj(x))./abs2(x)...)
Base.:/(x::Trivial, y::Trivial) = x*inv(y)

# Math
function Base.exp(x::Trivial)
    term = 1
    n = 0
    sum = 1.

    while abs(term) >= eps(sum|>abs)
        n += 1
        term *= x/n
        sum += term
    end
    return sum
end

# Output
function trivial_show(io::IO, x::Trivial{T}) where {T}
    if all(triplet(x) .== zero(T))
        print(io, 0)
        return
    end
    for (nat, symb) in zip(triplet(x), ["", " ╱", " ╲"])
        if nat == zero(T)
            continue
        end
        print(io, symb, nat)
    end
end
Base.show(io::IO, x::Trivial) = trivial_show(io, x)

# ⦦ ⦧

# ⧄ ⧅

# ╱ ╲ 

# ⨫ ⨬
