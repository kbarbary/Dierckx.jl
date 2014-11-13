#!/usr/bin/env julia
using Dierckx
using Grid

println("Grid.jl     ", Pkg.installed("Grid"))
println("Dierckx.jl  ", Pkg.installed("Dierckx"))
println()

# --------------------------------------------------------------
println("                         | Grid.jl      | Dierckx.jl")
println("-------------------------|--------------|--------------")

for n=[10, 100, 1000]
    @printf " eval 2-d grid %-9s | " "$(n)x$(n)"

    x = Float64[1:n]
    z = Float64[sin(i)+sin(j) for i in x, j in x]
    xi = Float64[1.5:0.5:n/2. + 1.]

    # Grid.jl
    g = CoordInterpGrid((1:n, 1:n), z, BCnearest, InterpQuadratic)
    g[xi, xi]
    t = @elapsed g[xi, xi]
    @printf "%9.5f ms | " t*1000.

    # Dierckx.jl
    spl = Spline2D(x, x, z; kx=2, ky=2)
    evalgrid(spl, xi, xi)
    t = @elapsed evalgrid(spl, xi, xi)
    @printf "%9.5f ms \n" t*1000.
end
