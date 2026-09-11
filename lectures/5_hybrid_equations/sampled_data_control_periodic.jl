using DiffEqCallbacks
using OrdinaryDiffEq
using ControlSystems
using LinearAlgebra
using Plots

#---------------------------------------------------------------------------------------------------------
# Discrete-time design and simulation
#---------------------------------------------------------------------------------------------------------

# System parameters
m₁ = 1      # Mass 1
m₂ = 0.1    # Mass 2
b₁₂ = 0.03  # Damping coefficient between the two masses
k₁₂ = 0.09  # Spring stiffness between the two masses
b₀₁ = 0.3   # Damping coefficient between the first mass and the wall
k₀₁ = 0.5   # Spring stiffness between the first mass and the wall

# State vector: x = [y'; y; d'; d], 
# where y is the position of the first mass, and d is the position of the second mass.

# State matrices
A = [-(b₁₂+b₀₁)/m₁ -(k₁₂+k₀₁)/m₁ b₁₂/m₁ k₁₂/m₁;
    1 0 0 0;
    b₁₂/m₂ k₁₂/m₂ -b₁₂/m₂ -k₁₂/m₂;
    0 0 1 0]
B = [1/m₁; 0; 0; 0]         # The force (the control input) acts on the first mass.
C = [0 0 0 1]               # We want to control the position of the second mass, d.
D = [0]

# State space model, first continuous-time, then ZOH-discretized
G = ss(A, B, C, D)
Ac, Bc, Cc, Dc = ssdata(G)  # Keeping the continuous-time matrices around for the sampled-data simulation below
h = 0.5                     # Sampling time
Gd = c2d(G, h)
A, B, C, D = ssdata(Gd)     # Rewriting the original matrices by those correspoinding to the discretized system

# Feedforward for reference tracking
N = [A-I B; C 0]\[0; 0; 0; 0; 1]
Nₓ = N[1:4]
Nᵤ = N[5]

# Setting up the cost function for the LQR
Q_y = 10.0 
Q = C'*Q_y*C
rank(obsv(A,Q))             # Check the observability of the system with the new Q matrix. Should be 4.
R = 1.0

# Solving the LQR problem
K = lqr(Discrete, A, B, Q, R)

# Introduce the reference into the state feedback scheme
y_ref(t) = 1.0(t>=1.5)          # Reference output (step function starting at t=1.5)
N = K*Nₓ .+ Nᵤ                  # Dot needed here because K*Nₓ is a 1x1 matrix, while Nᵤ is a scaler.
κ(x,t)  = -K*x .+ N*y_ref(t)    # Control law (u is a function of t and x)

# Simulate the closed-loop system in discrete time
t = 0:h:10                      # Time vector
x_init = [0, 0.0, 0, 0.0]       # Initial condition
res = lsim(G,κ,t,x0=x_init)     # Simulation results

# Plotting the results
p1 = plot(res, ylabel="Tracked output", lab="d", ploty=true, plotx=false, plotu=false, layout=1, sp=1, lw=2, seriestype=:steppost, color=palette(:default)[4])
plot!(t, y_ref.(t), lab="reference", lw=2, ls=:dash, layout=(3,1), sp=1, seriestype=:steppost)
p2 = plot(res, ylabel="State", lab=["y dot" "y" "d dot" "d"], ploty=false, plotx=true, plotu=false, layout=1, sp=1, lw=2, seriestype=:steppost, legend_position = :topright)
p3 = plot(res, ylabel="Control", lab="u", ploty=false, plotx=false, plotu=true, layout=1, sp=1, lw=2, seriestype=:steppost)
plot(p1, p2, p3, layout=(3,1))

#---------------------------------------------------------------------------------------------------------
# Continuous-time simulation for the sampled-data system (continuous-time plant and discrete-time control)
#---------------------------------------------------------------------------------------------------------

# Defining the set-based conditions and functions for the hybrid equations

fₚ(x,u) = Ac*x + Bc*u       # State equation for the continuous dynamics of the system.
hₚ(x,u) = Cc*x + Dc*u       # Output equation.

struct CtrlParams{T}        # Structure for storing the ZOH control and the last sampling time.
    u::Vector{T}
    t_last::T
end

#is_in_C(x,u,t) = (t < h)    # Actually not really needed, just a complement of D.
#is_in_D(x,u,t) = (t >= h)   # The system can jump only if the timer τ is greater than or equal to the sampling time h.

function f!(dx,x,p,t)       # Already in the format for the ODE solver.
    uₖ = p.u                # Control input at the last sampling time.
    dx .= fₚ(x,uₖ)          # Right hand side of the state equation with the ZOH control input.
end

function g(x,u,t)           # It only updates the discrete state variables (the ZOH control input and the last sampling time).
    p = CtrlParams(κ(x,t), t)
end

saved_values = SavedValues(Float64, Vector{Float64})   # (time type, saved-value type)
save_func(x, t, integrator) = copy(integrator.p.u)  

function sampling_condition(x,t,integrator)
    tₖ = integrator.p.t_last
    uₖ = integrator.p.u
    τ = t - tₖ              # Time elapsed since the last sampling time.  
    #return is_in_D(x,uₖ,τ)  # The system can jump only if the timer τ is greater than or equal to the sampling time h.
    return τ-h              # Returns a real and not just a binary value so that zero-crossing detection can be used.
end

function affect!(integrator)
    x = integrator.u        # Indeed, it is as notationally confusing as this: u is the state, 
    u = integrator.p.u      #                                                  the control is accessed through p.u.
    t = integrator.t
    #y = hₚ(x,u)
    integrator.p = g(x,u,t)
end

#sampling_cb = DiscreteCallback(condition,affect!)
sampling_cb = ContinuousCallback(sampling_condition,affect!)
saving_cb = SavingCallback(save_func, saved_values)
cb = CallbackSet(sampling_cb, saving_cb)    # Order matters: the sampling callback's affect! must update integrator.p.u
                                             # before the saving callback records it at that same instant.

# Setting up the simulation problem
p_init = CtrlParams([0.0], 0.0)
tspan = (first(t),last(t))
prob = ODEProblem(f!,x_init,tspan,p_init)

# Solving the simulation problem
sol = solve(prob,Tsit5(),callback=cb,dtmax=0.1) # ContinuousCallback more suitable here

# Plotting the solution
p4 = plot(sol.t, sol[4,:], ylabel="Tracked output", lab="d", ploty=true, plotx=false, plotu=false, layout=1, sp=1, lw=2, color=palette(:default)[4])
plot!(t, y_ref.(t), lab="reference", lw=2, ls=:dash, layout=(3,1), sp=1, seriestype=:steppost)
p5 = plot(sol, ylabel="State", lab=["y dot" "y" "d dot" "d"], ploty=false, plotx=true, plotu=false, layout=1, sp=1, lw=2, legend_position = :topright)
p6 = plot(saved_values.t, first.(saved_values.saveval), ylabel="Control", lab="u", lw=2, seriestype=:steppost)
plot(p4, p5, p6, layout=(3,1))
