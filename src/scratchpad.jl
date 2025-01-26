using Plots
include("trivial.jl")

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



