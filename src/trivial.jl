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
# Trivial(x::Real) = Trivial(x, zero(x))
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
triplet(x::Number) = triplet(Trivial(x))
Base.isreal(x::Trivial) = x.b == 0
Base.real(x::Trivial) = x.a + 0.5*x.b

Base.isnan(z::Trivial) = isnan.(duplet(z)) |> any
Base.isinf(z::Trivial) = isinf.(duplet(z)) |> any
Base.isfinite(z::Trivial) = isfinite.(duplet(z)) |> all

Base.:≈(a::Trivial, b::Trivial) = all(duplet(a) .≈ duplet(b))

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

∥(x::Trivial, y::Trivial) = -(verso(x)+recto(y))

Base.:*(x::Trivial, y::Trivial) = Trivial(x.a*y.a-x.b*y.b, (x.a + x.b)*y.b + x.b*y.a)
# Base.:*(x::Trivial, y::Trivial) = Trivial(reduce(.+, n.*triplet(op(y)) for (n, op) in zip(triplet(x), [x -> x, verso, recto]))...)
conj(x::Trivial) = Trivial(x.a+x.b, -x.b)
Base.abs2(x::Trivial{T}) where {T} = convert(T, x*conj(x))
Base.abs(x::Trivial) = sqrt(abs2(x))
Base.inv(x::Trivial) = Trivial(duplet(conj(x))./abs2(x)...)
Base.:/(x::Trivial, y::Trivial) = x*inv(y)
# Base.:/(x::Trivial, y::Trivial) = Trivial(((x.a*y.a + x.a*y.b + x.b*y.b, -x.a*y.b + x.b*y.a)./((y.a+y.b)^2-y.a*y.b))...)

# Math
function Base.exp(x::Trivial)
    term = one(Trivial)
    n = 0
    sum = 0.

    while abs(term) >= eps(sum|>abs)
        n += 1
        sum += term
        term *= x/n
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

# Relations
import Base.sign
function sign(x::Trivial)
    diffs = [x.a, x.b, -x.a-x.b]
    nzeros = sum(diffs.==0...)

    kind = nzeros == 0 ? :spiral : nzeros == 1 ? :isosceles : :zero
    if kind == :zero
        return (;kind, polarity=0, rotation=1)
    end

    polarity = 0
    rotation = 1
    if kind == :spiral
        npos = sum(diffs .> 0) # number of counterclockwise increases
        polarity = npos == 2 ? +1 : -1
        rotation = findfirst(triplet(x).==maximum(triplet(x)))
    elseif kind == :isosceles
        npos = sum(triplet(x) .> 0) # number of positive values in triplet(x)t
        polarity = npos == 2 ? -1 : +1
        # rotation = findfirst(triplet(x) .> 0)
        
        if polarity > 0
        rotation = findfirst(triplet(x).==maximum(triplet(x)))
        else
            rotation = mod(findfirst(triplet(x) .== 0), 3) + 1
        end
        # rotation = mod(findfirst(diffs .== 0) - 2, 3) + 1
    end 

    return (;kind, polarity, rotation)
end
sign(x::Trivial, y::Trivial, z::Trivial) = sign(x+recto(y)+verso(z))
sign(Trivial(1+0.1╱))

Base.:-(a::Number, b::Number, c::Number) = a+verso(b)+recto(c)

function dual(a::Number, b::Number, c::Number)
    i = triplet(a)
    j = triplet(b)
    k = triplet(c)
    return Trivial(i[1], j[1], k[1]), Trivial(i[2], j[2], k[2]), Trivial(i[3], j[3], k[3])
end

# ⦦ ⦧

# ⧄ ⧅

# ╱ ╲ 

# ⨫ ⨬
