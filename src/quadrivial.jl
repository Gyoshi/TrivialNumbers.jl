using StaticArrays, LinearAlgebra

struct Quadrivial{T<:Real} <: Number
    triplet :: MVector{3,T}
    diffs :: MVector{3,T}
end

# Constructors
function Quadrivial(triplet::T1, diffs::T2) where {T1<:AbstractArray, T2<:AbstractArray}
    i = T1 <: StaticArray ? 2 : 1
    j = T2 <: StaticArray ? 2 : 1
    if !all(triplet .≈ circshift(triplet, -1) + diffs) && !all(triplet - circshift(triplet, -1) .≈ diffs) # double check because 0 ≈ only zero
        throw(DomainError((triplet=triplet, diffs=diffs), "Difference array must be equal to the circular-sequential differences of the triplet."))
    end
    return Quadrivial{promote_type(T1.parameters[i], T2.parameters[j])}(triplet, diffs)
end
Quadrivial(x::AbstractArray) = Quadrivial(x, x-circshift(x, -1))
Quadrivial(a::Real, b::Real, c::Real) = Quadrivial([a,b,c])
Quadrivial(a::Real, b::Real, c::Real, d::Real) = Quadrivial(a-b, a-c, a-d)
Quadrivial(x::Real) = Quadrivial(x, x, x)

# Promotion
Base.convert(::Type{Quadrivial{T}}, z::Quadrivial{T}) where {T<:Real} = z
Base.convert(::Type{Quadrivial{T}}, z::Quadrivial) where {T<:Real} = Quadrivial(convert.(T, triplet(z))...)
Base.convert(::Type{Quadrivial{T}}, x::Real) where {T<:Real} = Quadrivial(convert(T, x), convert(T, x), convert(T, x))
(::Type{T})(z::Quadrivial) where {T<:Real} = (isreal(z) || isnan(z) ? triplet(z)[1] : throw(InexactError(:convert, T, z)))
Base.promote_rule(::Type{Quadrivial{T}}, ::Type{S}) where {T<:Real, S<:Real} = Quadrivial{promote_type(T,S)}
Base.promote_rule(::Type{Quadrivial{T}}, ::Type{Quadrivial{S}}) where {T<:Real, S<:Real} = Quadrivial{promote_type(T,S)}

# Util
Base.getproperty(x::Quadrivial, _::Symbol) = error("Quadrivial struct fields are private to avoid passing of mutables. Use `triplet()` instead.")
triplet(x::Quadrivial) = copy(getfield(x, :triplet))
diffstriplet(x::Quadrivial) = copy(getfield(x, :diffs))
function quadruplet(x::Quadrivial)
    result = (0, -triplet(x)[1], -triplet(x)[2], -triplet(x)[3])
    return result .- min(result...)
end
quadruplet(x::Number) = quadruplet(Quadrivial(x))

verso(x::Quadrivial) = Quadrivial(triplet(x)[SVector(2,3,1)], diffstriplet(x)[SVector(2,3,1)])
recto(x::Quadrivial) = Quadrivial(triplet(x)[SVector(3,1,2)], diffstriplet(x)[SVector(3,1,2)])
verso(x::Number) = verso(Quadrivial(x))
recto(x::Number) = recto(Quadrivial(x))

Base.isreal(x::Quadrivial) = reduce(≈, triplet(x))
Base.real(x::Quadrivial{T}) where T = convert(T, triplet(x + verso(x) + recto(x))[1])/3
Base.imag(x::Quadrivial) = x - real(x)

Base.isnan(z::Quadrivial) = isnan.(triplet(z)) |> any
Base.isinf(z::Quadrivial) = isinf.(triplet(z)) |> any
Base.isfinite(z::Quadrivial) = isfinite.(triplet(z)) |> all

Base.:≈(a::Quadrivial, b::Quadrivial) = all(triplet(a) .≈ triplet(b))

# Operators
Base.:+(x::Quadrivial, y::Quadrivial) = Quadrivial(triplet(x)+triplet(y), diffstriplet(x)+diffstriplet(y))
Base.:-(x::Quadrivial) = Quadrivial(-triplet(x), -diffstriplet(x))
Base.:-(x::Quadrivial, y::Quadrivial) = x + -y

const i = Quadrivial(-1,0,0)
const j = Quadrivial(0,-1,0)
const k = Quadrivial(0,0,-1)

# function mulmatrix(x::Quadrivial) # 7 times slower when used with mul
#     (a, b, c) = triplet(x)
#     (α, β, γ) = diffstriplet(x)
#     return [a β -β; -γ b γ; α -α c]
# end

function mulmatrix(x::Quadrivial)
    t = triplet(x)
    d = diffstriplet(x)
    return vcat(t,d,-d)[SA[1 5 8; 9 2 6; 4 7 3]] 
end

Base.:*(x::Real, y::Quadrivial) = Quadrivial(x*triplet(y), x*diffstriplet(y))
Base.:*(x::Quadrivial, y::Real) = y*x
function Base.:*(x::Quadrivial, y::Quadrivial)
    w = mulmatrix(x)*mulmatrix(y)
    diffs = [w[3,1]-w[3,2], w[1,2]-w[1,3], w[2,3]-w[2,1]]/2
    # return Quadrivial((diag(w) + diffs + diag(w)[SA[2,3,1]])/2, diffs)
    return Quadrivial(diag(w), diffs)
end

# @btime *(x,y) setup=(x = 1+0.8i, y = 0.5-j)
# @btime *(x,y) setup=(x = 1+0.8i, y = 0.5-j)

inner_product(x, y, z) = x*recto(y)*verso(z)
abs3(x::Quadrivial) = real(x*verso(x)*recto(x))
abs3(x::Number) = abs3(Quadrivial(x))
Base.abs(x::Quadrivial) = cbrt(abs3(x))

Base.:/(x::Quadrivial, y::Real) = Quadrivial(triplet(x)/y, diffstriplet(x)/y)
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

x = 1 + 0.8i
exp(x)

trip = [1.4109347442680748e-12, 0.0002733125149911817, 0.00027331251499118165]
diffs = [-0.0002733125135802469, 0.0, 0.00027331251358024696]
circshift(trip, -1) + diffs - trip

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

