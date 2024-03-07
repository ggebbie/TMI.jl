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

close("all")
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
subplot(1,2,1) 
for i in 1:size(sol[c(t,x)])[1]
    plot(sol[x], sol[c(t,x)][i, :])
end

#=
sol2 = get_sol(2.0, 0.0)
subplot(1,2,2) 
for i in 1:size(sol2[c(t,x)])[1]
    plot(sol2[x], sol2[c(t,x)][i, :])
end
=#

function make_field(sol, sym, ln, un)
    v = sol[c(t,x)][end, :]
    ss = Array{Float64}(undef, 1,1,length(v))
    ss[1,1,:] = v
    wet = BitArray(undef, 1,1,length(v))
    wet[1,1, :] .= fill(true, length(v))
    interior = wet
    interior[1,1,1] = 0 
    γ = TMI.Grid([0.0], [0.0], collect(sol[x]), wet, interior)
    return Field(ss, γ,sym, ln, un), γ
end

θ̄, γ = make_field(sol, :θ, "temperature", "K")

y = (θ = θ̄,)
     #S = make_field(sol2, :S, "salinity", "‰"))
neighbors = TMI.neighbors(y, γ).tracer

@time m̃ = massfractions(y);
Ã = watermassmatrix(m̃, γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(y.θ)

## reconstruct temperature
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ã,bθ,γ)

# compare to c.θ
Base.maximum(y.θ.tracer - θ̃.tracer)
Base.minimum(y.θ.tracer - θ̃.tracer)

y.θ.tracer[1, 1, :]
θ̃.tracer[1,1,:]
