using Plots
include("trivial.jl")

(1+2╱)/sqrt(3) |> abs
Trivial(0, 1, -1) |> abs2

x = (1+2╱)
y = (2+3╱)
z = Trivial(2,-2)

x,y,z
dual(x, y, z)
-x, -y, -(x∥y)
dual(-x, -y, -(x∥y))

sign(x, y, x∥y)
sign(dual(x, y, x∥y)...)
sign(-x, -y, z)
sign(dual(-x, -y, z)...)
# shouldn't dual preserve minus (or at least sign), from a conceptual viewpoint?

dual(-2, Trivial(0,5,3), -4)

-(x,y,z)
-(dual(x,y,z)...)
dual(x,y,z)

inv(1 + 2╱)
x/y
y/y*x
x∥y

# exp
exp(Trivial(-1,1)) 
exp(╱) |> abs
exp(0.5+╲) |> abs ≈ 1

d = 0:0.01:10
scatter(triplet.(exp.(d.*(╱))), aspect_ratio=:equal)
    plot!(xlims=(0,1), ylims=(0,1), zlims=(0,1), aspect_ratio=:equal)

scatter(triplet.(exp.(d.*(0.5 + ╱))), camera=(135,35))
scatter(duplet.(exp.(d.*(0.5 + ╱))), aspect_ratio=:equal)
    
plot(exp.(d.*(0.5+╲)).|>abs)
plot(exp.(d.*(0.5+╲)).|>real)
plot(exp.(d.*(cis(deg2rad(-120))).|>Base.real))

x = 0.5+╱
exp(x/abs(x)) 
-Trivial(1)
exp(x)
exp(x) |> abs
abs(x)
(x/abs(x))^2
(0.5 + ╲)^2
cos(pi/6)
abs(sqrt(3)/2 + 0.5im)

# unit circle

circle = Set()
for i = 1:20, j = 0:i
    for m = 1:20, n = 0:m
        x = Trivial(j//i, n//m)
        if abs2(x) == 1
            push!(circle, x)
        end
    end
end
circle

as = bs = 1:100
function norm(a, b)
    x = Trivial(a, b, 0)
    return triplet(x/abs(x))[1:2]
end
unitpoints = [norm(a, b) for a in as for b in bs]
scatter(unitpoints, aspect_ratio=:equal)
scatter!([triplet(Trivial(a, b, 0))[1:2] for a in as for b in bs])

function scaleto120(a,b)
    angle = atan(a,b)
    angle *= 120/90
    return hypot(a,b) .* (sin(deg2rad(0)+angle), cos(deg2rad(0)+angle))
end

scatter!([scaleto120(norm(a, b)...) for a in as for b in bs], aspect_ratio=:equal)


# Sign 
Trivial(1, 0) |> sign
Trivial(1,-1) |> sign
Trivial(2,-1) |> sign
Trivial(0,0) |> sign

as = bs = 0:12
turns = 0:2
xs = [Trivial(circshift([a, b, 0], r)...) for a in as for b in bs for r in turns]
triangles = [(x, Trivial(3), Trivial(0,1,0)) for x in xs]
kind = [sign(t...).kind |> x -> findfirst([:equal, :isosceles, :spiral] .== x) for t in triangles]
polarity = [sign(t...).polarity for t in triangles]
rotation = [sign(t...).rotation for t in triangles]
scatter([Tuple(circshift([triplet(x)...], 2)) for x in xs], zcolor=rotation, camera=(135, 35))
    plot!(xlims=(0,12), ylims=(0,12), zlims=(0,12), aspect_ratio=:equal)

sign(Trivial(1,1,0))
sign(Trivial(1,1.1,0))

Trivial(3)∥Trivial(0,1,0)

x |> sign
conj(x) |> sign

### Quadrivial numbers
include("quadrivial.jl")

x = 3+i + 1+2k
x = 3-j+k

abs(x)
abs(x * x) |> sqrt
abs(1.5+i+j)

exp(i) # cosh(1) + sinh(1)i
exp(-i)
exp(-i) |> abs
ℯ

exp(i+j) # cosh(1)/ℯ - k⋅sinh(1)/ℯ
i+j |> abs
i+j |> exp |> abs
(1+i)^2 |> abs

exp(1+i)
abs(1+i)
exp(1+i) |>abs
exp(0.5 + i)

x
3/4*Quadrivial(real(x), real(x*i), real(x*j), real(x*k))

# Plots
using Plots
gr(size=(700,600),markerstrokewidth=0,markersize=4, linewidth=2.5)

xs = ys = zs = -0.:0.1:1
xarray = [x for x in xs for y in ys for z in zs if x+y+z <= 3]
yarray = [y for x in xs for y in ys for z in zs if x+y+z <= 3]
zarray = [z for x in xs for y in ys for z in zs if x+y+z <= 3]

points = [Quadrivial(x, y, z) for x in xs for y in ys for z in zs if x+y+z <= 3]

scatter(xarray, yarray, zarray, zcolor=abs.(points), camera=(135, 35))
scatter(xarray, yarray, zarray, zcolor=abs.(exp.(points)), camera=(10,-20))
scatter(xarray, yarray, zarray, zcolor=real.(exp.(points)), camera=(10,-20))

# taking first term is not sufficient for identity: need real(⋅) value
scatter(abs.(exp.(points)), exp.((real.(points))))
scatter(abs.(exp.(points)), exp.((quadruplet(p)[1] for p in points)))
scatter(xarray, yarray, zarray, zcolor=real.(i.*exp.(points .- real.(points))), camera=(134,35))

magnitude = abs.(points)
scatter(real.(points), abs.(exp.(points)), zcolor=magnitude)

scaled_points = [p/abs(abs(p)) for p in points]
plot(abs.(scaled_points))
# scaled_points[(abs.(scaled_points) .≈ 1)] .= one(Quadrivial)
scatter([s.a for s in scaled_points], [s.b for s in scaled_points], [s.c for s in scaled_points], zcolor=real.(points), camera = (-45, 350), xlims = (-5,5), ylims = (-5,5), zlims = (-5,5))

ω = one(Quadrivial) + 0.1i
t = 0:0.01:10pi
y = exp.(ω*t)
plot(t, abs.(y))
plot(t, [real.(y), real.(y.*i), real.(y.*j), real.(y.*k)], labels=["Re" "Im_i" "Im_j" "Im_k"])

y = exp.(ω*t .- real.(ω*t))
plot(t, [real.(y), real.(y.*i), real.(y.*j), real.(y.*k)], labels=["Re" "i" "j" "k"])

y = cos.(ω*t)
# plot(t, [real.(y), real.(y.*i), real.(y.*j), real.(y.*k)], labels=["Re" "i" "j" "k"])
plot(t, [[p.a for p in y], [p.b for p in y], [p.c for p in y]], labels=["-i" "-j" "-k"])

y = sin.(ω*t)
plot(t, [real.(y), real.(y.*i), real.(y.*j), real.(y.*k)], labels=["Re" "i" "j" "k"])
plot(t, [[p.a for p in y], [p.b for p in y], [p.c for p in y]], labels=["-i" "-j" "-k"])

xs = 0:.001:1
tween = [x .- i .* (1-x) for x in xs]
# tween = [2x - 1 for x in xs]
plot(xs, abs3.(tween))
plot!(xs, cbrt.(2(xs.-0.5)))

# para-conjugates
@assert abs3(x) == (Quadrivial(x.a, x.c, x.b)*Quadrivial(x.b, x.a, x.c)*Quadrivial(x.c, x.b, x.a))
@assert real(x) == (Quadrivial(x.a, x.c, x.b)+Quadrivial(x.b, x.a, x.c)+Quadrivial(x.c, x.b, x.a))/3


log(abs(exp(i)))
i |> real

abs(1+0.9999i)
abs(1-0.99999999i)

triplet(1+0.9999i)
triplet(1-0.99999999i)

a = imag(i)/imag(i)*i
b = imag(j)/imag(i)*i
c = imag(k)/imag(i)*i

imag(a*a)/imag(i)

real(a)
a^2
a^2 |> imag

## Orthonormal basis

e0 = 3/4*Quadrivial(1, -1/3, -1/3, -1/3)
3*imag(2i-j) |> quadruplet |> x -> x.-x[1]
e1 = Quadrivial(0, 5, -4, -1)
e2 = Quadrivial(0, -1, -2, 3)
abs3(e1), abs3(e2)
e1 = e1/(2*cbrt(20))
e2 = e2/(2*cbrt(6))

x = (1 + e1 - e2) 

y = 4 + j - 2k
y = [quadruplet(y)...]

C = [ # collects 4 values into 3 (setting j to zero)
    1 0 -1 0;
    0 1 -1 0;
    0 0 -1 1
]
Cinv = [ # collects 3 values into 4 (assuming j is zero)
    1 0 0;
    0 1 0;
    0 0 0;
    0 0 1
]

A = hcat( # projects onto Orthonormal basis
    [1, 0, 0, 0],
    [quadruplet(e1)...],
    [quadruplet(e2)...],
) |> transpose
A = A*Cinv

ybase = A*C*y
recovered_y = Cinv*inv(A)*ybase
Quadrivial(recovered_y...) .≈ Quadrivial(y...)

Cinv*inv(A)[:,1]

# conjugate basis
f(x) = abs3(3x + i + j + k) - (abs3(x+i) + abs3(x+j) + abs3(x+k))
plot(-1:0.01:1, f)

g = 0.6265382932707997 # ≈ 1/(1 - 1/(1 + sqrt(2))^(1/3) + (1 + sqrt(2))^(1/3))
f(g)

e1 = g+i
e2 = g+j
e3 = g+k

abs3(e1 + e2 + e3) - (abs3(e1) + abs3(e2) + abs3(e3))# duh
e1/abs(e1)
e1*e2*e3
e1/abs(e1)+e2/abs(e2)+e3/abs(e3) ≈ cbrt(3)
abs(e1-e2/3-e3/3)
real(e1/abs(e1))
1/cbrt(9)
imag(e1/abs(e1)) |> abs
e1/abs(e1)
e1-e2 |> abs # cool property!

(1+1e-15)*e1 - e2 |> abs
(1+1e-15) + i |> abs

f(x) = abs3(3x - i - j - k) - (abs3(x-i) + abs3(x-j) + abs3(x-k))
plot(-1:0.01:1, f)
f(-g)

q1 = e1-e2
q2 = e2-e3
q3 = e3-e1

q1+q2+q3

e1 = (g+i)/abs(g+i)
e2 = (g+j)/abs(g+j)
e3 = (g+k)/abs(g+k)

# inner_product(x, y, z) = (abs3(x + y + z) - (abs3(x) + abs3(y) + abs3(z)))/24
# inner_product(e1,e2,e3)
# inner_product(-e1,-e2,-e3)
# inner_product(e3,e2,e1)



1/e1-1/e2
e1-e2

# orthogonal, real + 2x singular basis
# when is a quadrivial number singular? when a=b+c, b=c+a, or c=a+b
# not a "strict orthogonality" in the sense that the orthogonality is only valid for a specific relative scaling 

Quadrivial(-1,1,0)
Quadrivial(0,1,-1)
f(k) = abs3(1 + Quadrivial(-k,k,0) + Quadrivial(0,k,-k)) - (1 + abs3(Quadrivial(-k,k,0)) + abs3(Quadrivial(0,k,-k)))
plot(f)

f(-0.75)

-j+k |> triplet
e0 = one(Quadrivial)
e1 = 0.75(i-j)
e2 = 0.75(i-k)

e0+e1+e2 |> abs3
e1+e2 |> abs3

e1+e2 |> abs

# Is there an orthogonal basis that remains orthogonal under independent scaling of basis vectors?

f(q, x, y) = inner_product(1, q*Quadrivial(x,y,-x-y),0)
plot(-1:0.01:1, x-> f(2,x,2))

2Quadrivial(2,1,3) |> triplet

Quadrivial(1,50, -1-50)