#!/usr/bin/env julia
using Dierckx
using Grid

macro timen(ex, n)
    quote
        local t0 = time_ns()
        for i = 1:$(esc(n))
            local val = $(esc(ex))
        end
        local t1 = time_ns()
        (t1 - t0) / 1.e9 / $(esc(n))
    end
end

function fmt_time(t)
    if t < 1.e-6
        return @sprintf "%7.2f ns" t*1.e9
    elseif t < 1.e-3
        return @sprintf "%7.2f us" t*1.e6
    elseif t < 1.
        return @sprintf "%7.2f ms" t*1.e3
    end
    @sprintf "%7.2f s" t
end 

println("Grid.jl     ", Pkg.installed("Grid"))
println("Dierckx.jl  ", Pkg.installed("Dierckx"))
println()
println("| benchmark               | Grid.jl    | Dierckx.jl | ratio |")
println("|-------------------------|------------|------------|-------|")

for n=[10, 1000, 100_000]
    @printf "| eval 1-d grid n=%-7s | " "$(n)"
    nloops = div(10_000_000, n)
    x = Float64[1:n;]
    y = sin(x)
    xi = Float64[1.5:0.5:(n/2. + 1.);]

    g = CoordInterpGrid(1:n, y, BCnearest, InterpQuadratic)
    g[xi]
    t_grid = @timen g[xi] nloops
    print(fmt_time(t_grid), " | ")
    
    # Dierckx.jl
    spl = Spline1D(x, y; k=2)
    evaluate(spl, xi)
    t_dierckx = @timen evaluate(spl, xi) nloops
    print(fmt_time(t_dierckx), " | ")

    # ratio
    @printf "%5.2f |\n" t_grid/t_dierckx
end

for n=[10, 100, 1000]
    @printf "| eval 2-d grid %-9s | " "$(n)x$(n)"

    nloops = div(10_000_000, n*n)

    x = Float64[1:n;]
    z = Float64[sin(i)+sin(j) for i in x, j in x]
    xi = Float64[1.5:0.5:(n/2. + 1.);]

    # Grid.jl
    g = CoordInterpGrid((1:n, 1:n), z, BCnearest, InterpQuadratic)
    g[xi, xi]
    t_grid = @timen g[xi,xi] nloops
    print(fmt_time(t_grid), " | ")
    
    # Dierckx.jl
    spl = Spline2D(x, x, z; kx=2, ky=2)
    evalgrid(spl, xi, xi)
    t_dierckx = @timen evalgrid(spl, xi, xi) nloops
    print(fmt_time(t_dierckx), " | ")

    # ratio
    @printf "%5.2f |\n" t_grid/t_dierckx
end
