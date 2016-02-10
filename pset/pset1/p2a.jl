#!/usr/bin/julia -f

function calc_eigvals(t)
    H = Float64[0 -1 -t 0 -t -1 -t -1
                -1 0 -1 -t -1 -t 0 -t
                -t -1 0 -1 -t 0 -t -1
                0 -t -1 0 -1 -t -1 -t
                -t -1 -t -1 0 -1 -t 0
                -1 -t 0 -t -1 0 -1 -t
                -t 0 -t -1 -t -1 0 -1
                -1 -t -1 -t 0 -t -1 0]
    eigvals(H)
end

function all_eigvals(ts)
    res = Matrix{Float64}(8, length(ts))
    @inbounds for i in 1:length(ts)
        D = calc_eigvals(ts[i])
        for j in 1:8
            res[j, i] = D[j]
        end
    end
    res
end

ts = linspace(0, 10, 10000)
res = all_eigvals(ts)

using PyPlot
for i in 1:8
    plot(ts, res[i, :])
end
grid()
xlabel("\$t'/t\$")
ylabel("\$f_j\$")
savefig("p2c.png")
