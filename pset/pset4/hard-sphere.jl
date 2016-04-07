#!/usr/bin/julia -f

using GSL
using PyPlot

function gen_scatter(θs, δs, kr)
    lmax = length(δs) - 1
    calc = θ->begin
        a = 0.0im
        @inbounds for l in 0:lmax
            δ = δs[l + 1]
            a += (2l + 1) * exp(im * δ) * sin(δ) * sf_legendre_Pl(l, cos(θ))
        end
        abs2(a)
    end
    println(4π / kr^2 * sum((2l + 1) * sin(δs[l + 1])^2 for l in 0:lmax))
    Float64[calc(θ) / kr^2 for θ in θs]
end

function hard_sphere_δs(n, kr)
    calc = l->begin
        tanδ = besselj(l + 0.5, kr) / bessely(l + 0.5, kr)
        atan(tanδ)
    end
    Float64[calc(l - 1) for l in 1:n]
end

θs = linspace(0, 2π, 1000)
rs = gen_scatter(θs, hard_sphere_δs(3, 1), 1)

ax = subplot(111, projection="polar")
ax[:plot](θs, rs)
ax[:grid](true)
savefig("p5_3wave_kr1.png")
show()
