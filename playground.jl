using SumOfSquares
using DynamicPolynomials
using MosekTools
using CSDP

optimizer = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(optimizer)
@polyvar x[1:2] 

p = 1;

f = [ x[2],
     -x[1] + (p/3)*x[1]^3 - x[2]]

g₁ = -(x[1]+1)^2 - (x[2]+1)^2 + 0.16  # 𝒳ᵤ = {x ∈ R²: g₁(x) ≥ 0}
h₁ = -(x[1]-1.5)^2 - x[2]^2 + 0.25    # 𝒳₀ = {x ∈ R²: h₁(x) ≥ 0}

X = monomials(x, 0:4)
@variable(model, B, Poly(X))

ε = 0.001
@constraint(model, B >= ε, domain = @set(g₁ >= 0))
@constraint(model, B <= 0, domain = @set(h₁ >= 0))

using LinearAlgebra # Needed for `dot`
dBdt = dot(differentiate(B, x), f)
@constraint(model, -dBdt >= 0)

optimize!(model)

import DifferentialEquations, Plots, ImplicitPlots
function phase_plot(f, B, g₁, h₁, quiver_scaling, Δt, X0, solver = DifferentialEquations.Tsit5())
    X₀plot = ImplicitPlots.implicit_plot(h₁; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="X₀", linecolor=:blue)
    Xᵤplot = ImplicitPlots.implicit_plot!(g₁; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="Xᵤ", linecolor=:teal)
    Bplot  = ImplicitPlots.implicit_plot!(B; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="B = 0", linecolor=:red)
    Plots.plot(X₀plot)
    Plots.plot!(Xᵤplot)
    Plots.plot!(Bplot)
    ∇(vx, vy) = [fi(x[1] => vx, x[2] => vy) for fi in f]
    ∇pt(v, p, t) = ∇(v[1], v[2])
    function traj(v0)
        tspan = (0.0, Δt)
        prob = DifferentialEquations.ODEProblem(∇pt, v0, tspan)
        return DifferentialEquations.solve(prob, solver, reltol=1e-8, abstol=1e-8)
    end
    ticks = -5:0.5:5
    X = repeat(ticks, 1, length(ticks))
    Y = X'
    Plots.quiver!(X, Y, quiver = (x, y) -> ∇(x, y) / quiver_scaling, linewidth=0.5)
    for x0 in X0
        Plots.plot!(traj(x0), idxs=(1, 2), label = nothing)
    end
    Plots.plot!(xlims = (-2, 3), ylims = (-2.5, 2.5), xlabel = "x₁", ylabel = "x₂")
end

phase_plot(f, value(B), g₁, h₁, 10, 30.0, [[x1, x2] for x1 in 1.2:0.2:1.7, x2 in -0.35:0.1:0.35])


x = [1, 3, 4, 5]
y = [2, 18, 32, 50]

k = -2
p(x) = 2x^2+k*(x-1)*(x-3)*(x-4)*(x-5)


yp = p.(x)
y6 = p(6)
