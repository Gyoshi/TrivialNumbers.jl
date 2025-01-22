include("trivial.jl")

x = (1+2╱)
y = (2+3╱)

inv(1 + 2╱)
x/y
y/y*x

# exp
exp(Trivial(-1,1)) 
exp(╱) |> abs

d = 0:0.01:10

scatter(exp.(d.*(╱)))

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

using Plots
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
