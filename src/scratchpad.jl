using Plots
include("trivial.jl")

x = (1+2╱)
y = (2+3╱)

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

scatter(triplet.(exp.(d.*(0.5 + ╱))))
plot(exp.(d.*(0.5+╲)).|>abs)
plot(exp.(d.*(0.5+╲)).|>real)
plot(exp.(d.*(cis(deg2rad(-120))).|>Base.real))

x = 0.5+╲
exp(x)
exp(x) |> abs

(0.5 + ╱)^2

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
a = b = 2
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
