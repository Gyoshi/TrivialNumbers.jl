using StaticArrays, LinearAlgebra

struct Quinquevial{T<:Real} <: Number
    quadruplet :: MVector{4,T}
end

# Constructors
function Quinquevial(quadruplet::T1) where {T1<:AbstractArray}
    i = T1 <: StaticArray ? 2 : 1
    return Quinquevial{T1.parameters[i]}(quadruplet)
end
Quinquevial(a1::Real, a2::Real, a3::Real, a4::Real) = Quinquevial([a1,a2,a3,a4])
Quinquevial(a0::Real, a1::Real, a2::Real, a3::Real, a4::Real) = Quinquevial([a1-a0, a2-a0, a3-a0, a4-a0])

Quinquevial(x::Real) = Quinquevial(repeat([-x], 4))

# Promotion
Base.convert(::Type{Quinquevial{T}}, z::Quinquevial{T}) where {T<:Real} = z
Base.convert(::Type{Quinquevial{T}}, z::Quinquevial) where {T<:Real} = Quinquevial(convert.(T, quadruplet(z))...)
Base.convert(::Type{Quinquevial{T}}, x::Real) where {T<:Real} = Quinquevial(x)
(::Type{T})(z::Quinquevial) where {T<:Real} = (isreal(z) || isnan(z) ? quadruplet(z)[1] : throw(InexactError(:convert, T, z)))
Base.promote_rule(::Type{Quinquevial{T}}, ::Type{S}) where {T<:Real, S<:Real} = Quinquevial{promote_type(T,S)}
Base.promote_rule(::Type{Quinquevial{T}}, ::Type{Quinquevial{S}}) where {T<:Real, S<:Real} = Quinquevial{promote_type(T,S)}

# Util
Base.getproperty(x::Quinquevial, _::Symbol) = error("Quinquevial struct fields are private to avoid passing of mutables. Use `quadruplet()` instead.")
quadruplet(x::Quinquevial) = copy(getfield(x, :quadruplet))
diffsquadruplet(x::Quinquevial) = copy(getfield(x, :diffs))
function quintuplet(x::Quinquevial)
    q4 = quadruplet(x)
    result = (0, q4[1], q4[2], q4[3], q4[4])
    return result .- min(result...)
end
quintuplet(x::Number) = quintuplet(Quinquevial(x))

lconj(x::Quinquevial) = Quinquevial(quadruplet(x)[SVector(3,1,4,2)])
bconj(x::Quinquevial) = Quinquevial(quadruplet(x)[SVector(4,3,2,1)])
rconj(x::Quinquevial) = Quinquevial(quadruplet(x)[SVector(2,4,1,3)])
lconj(x::Number) = Quinquevial(lconj(Quinquevial(x)))
bconj(x::Number) = Quinquevial(bconj(Quinquevial(x)))
rconj(x::Number) = Quinquevial(rconj(Quinquevial(x)))
quinq_conjugates = [x->x, lconj, bconj, rconj]

Base.isreal(x::Quinquevial) = reduce(≈, quadruplet(x))
Base.real(x::Quinquevial{T}) where T = convert(T, -quadruplet(sum(c(x) for c in quinq_conjugates)) |> first)/4
Base.imag(x::Quinquevial) = x - real(x)

Base.isnan(z::Quinquevial) = isnan.(quadruplet(z)) |> any
Base.isinf(z::Quinquevial) = isinf.(quadruplet(z)) |> any
Base.isfinite(z::Quinquevial) = isfinite.(quadruplet(z)) |> all

Base.:≈(a::Quinquevial, b::Quinquevial) = all(quadruplet(a) .≈ quadruplet(b))

# Operators
Base.:+(x::Quinquevial, y::Quinquevial) = Quinquevial(quadruplet(x)+quadruplet(y))
Base.:-(x::Quinquevial) = Quinquevial(-quadruplet(x))
Base.:-(x::Quinquevial, y::Quinquevial) = x + -y

const r = Quinquevial(1,0,0,0)

# function mulmatrix(x::Quinquevial) # 7 times slower when used with mul
#     (a, b, c) = quadruplet(x)
#     (α, β, γ) = diffsquadruplet(x)
#     return [a β -β; -γ b γ; α -α c]
# end

function mulmatrix(x::Quinquevial)
    a, b, c, d = quadruplet(x)
    return [    -d  d-c c-b b-a
            ;   a-d -c  d-b c-a
            ;   b-d a-c -b  d-a
            ;   c-d b-c a-b -a
            ]
end

Base.:*(x::Real, y::Quinquevial) = Quinquevial(x*quadruplet(y))
Base.:*(x::Quinquevial, y::Real) = y*x
function Base.:*(x::Quinquevial, y::Quinquevial)
    w = mulmatrix(x)*mulmatrix(y)
    return bconj(Quinquevial(-diag(w)))
end
(1+r)^2
(1 + 2r + r^2)

inner_product(w, x, y, z) = w*lconj(x)*bconj(y)*rconj(z)
abs4(x::Quinquevial) = real(inner_product(x,x,x,x))
abs4(x::Number) = abs4(Quinquevial(x))
Base.abs(x::Quinquevial) = sqrt(sqrt(abs4(x)))

Base.:/(x::Quinquevial, y::Real) = Quinquevial(quadruplet(x)/y)
Base.inv(x::Quinquevial) = inner_product(1/abs4(x), x, x, x)
Base.:/(x::Quinquevial, y::Quinquevial) = x*inv(y)

# Math
function Base.exp(x::Quinquevial)
    term = one(Quinquevial)
    n = 0
    rsum = 0.

    for i in 1:1000
        n += 1
        rsum += term
        term *= x/n        
    end
    return rsum
end

x = r
exp(x)

function Base.cos(x::Quinquevial)
    x2 = x*x
    term = one(Quinquevial)
    n = 0
    rsum = 0.

    for i in 1:1000
        n += 2
        rsum += term
        term *= -x2/(n*(n-1))      
    end
    return rsum
end
function Base.sin(x::Quinquevial)
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
function Quinquevial_show(io::IO, x::Quinquevial{T}) where {T}
    first_plus = true
    if all(quintuplet(x) .== zero(T))
        print(io, 0)
        return
    end
    coefficients = quintuplet(x)
    reps, unique_index = [count(==(element),coefficients) for element in unique(coefficients)]|> findmax
    if reps > 2
        coefficients = coefficients .- unique(coefficients)[unique_index]
    end
    for (nat, symb) in zip(coefficients, ["", "r", "r²", "r⁻²", "r⁻¹"]) 
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
Base.show(io::IO, x::Quinquevial) = Quinquevial_show(io, x)
