#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

copy-paste from ex8
1d copy-paste from TVCR#bh-1dA/1dA_separate.jl 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
#using GeoPythonPlot # will load optional extension
using PyPlot
using ModelingToolkit
using MethodOfLines
using OrdinaryDiffEq
using DomainSets

#close("all")
@parameters t x κ u
@variables c(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

function get_sol(bc1, bc2) 
    T = 100.0 #yr 
    L = 1000.0 #km 
    Δx = 10.0 #km 
    Δt = 10.0 #yr 
    û = L / 15 # km/yr 
    κ̂ = L^2 / 15 #km^2
    eq = Dt(c(t, x)) ~ κ * Dxx(c(t,x)) - u * Dx(c(t, x))
    bcs = [c(t,0) ~ bc1,
           c(t,L)~bc2,
           c(0,x)~0.0]
    domains = [t ∈ Interval(0.0, T),
               x ∈ Interval(0.0, L)]
    @named pdesys = PDESystem(eq, bcs, domains, [t,x],[c(t,x)], [u=>û, κ=>κ̂])

    discret = MOLFiniteDifference([x=>Δx, t=>Δt])
    prob = discretize(pdesys, discret)
    sol = solve(prob, Tsit5())
    return sol
end

sol = get_sol(1.0, 0.0)

figure()
for i in 1:size(sol[c(t,x)])[1]
    plot(sol[x], sol[c(t,x)][i, :], label = "t = " * string(i), linewidth = i, alpha = (12 - i) / 11)
end
legend()
#=
sol2 = get_sol(2.0, 0.0)
subplot(1,2,2) 
for i in 1:size(sol2[c(t,x)])[1]
    plot(sol2[x], sol2[c(t,x)][i, :])
end
=#

function make_field(sol, sym, ln, un)
    v = sol[c(t,x)][end, :]
    wet = fill(true, length(v))
    interior = copy(wet)
    interior[end] = 0
    interior[1] = 0 
    γ = TMI.Grid((collect(sol[x]),), convert(BitVector,wet), convert(BitVector, interior), (false,),  [CartesianIndex(1,),CartesianIndex(-1,)])
    return Field(v, γ, sym, ln, un)
end

θ̄ = make_field(sol, :θ, "temperature", "K")
#S̄ = make_field(get_sol(2.0, 0.0), :S, "salinity", "‰")

y = (θ = θ̄,)
#     S = S̄)
 
#how closely are we fitting data? 
w = (θ =  0.01,
#    S = 0.001,
 #   δ¹⁸O = 0.05,
  #  P★ = 0.05,
   # δ¹³C★ = 0.05
)


@time m̃ = massfractions(y, w);

Ã = watermassmatrix(m̃, y.θ.γ)
γ = θ̄.γ

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
dim =1
b = (west = TMI.getboundarycondition(θ̄, dim, 1, γ),
     east = TMI.getboundarycondition(θ̄, 1, length(θ̄.tracer), γ))

## reconstruct temperature
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ãlu,b,γ)

# compare to c.θ
Base.maximum(y.θ.tracer - θ̃.tracer)
Base.minimum(y.θ.tracer - θ̃.tracer)

figure();
for θ in [θ̄, θ̃]
    plot(θ.tracer)
end

