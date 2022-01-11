#= placeholder for transient simulation
#how does a patch of tracer (concentration = 1) evolve throughout time?
=#

using Revise, TMI, Interpolations, PyPlot, NaNMath, DifferentialEquations, LinearAlgebra

#load TMI data
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
#here we will only use
#γ: corresponding lat/lon values for L 
#L: tendency matrix, essentially A (pathway matrix) scaled by residence time (tau) 
#each entry of L tells us how each point in the ocean is connected to its neighbors

#In MATLAB code: make_initial_condition - either rectangle or predefined region
#here we'll assume rectangle 
latbox = [50,60]
lonbox = [-50,0]
d = surfacepatch(lonbox, latbox, γ) #vector(same dims as γ, global)  describing surface patch in terms of γ

#γ.wet tells us, at each depth, which points are wet (binary mapping)
#we'll make the first level of γ.wet into a bit array, and access the first level at d at each index. this essentially flattens it into the same format as B 
dsfc =  d[:,:,1][γ.wet[:,:,1]]

#following make_initial_conditions.m
#initial conditions are the surface patch = 1, propagated down to the bottom of the mixed layer, which we get from B 
c0 = B * dsfc 

#Fixed euler timestep approximation
c = c0
Δt = 1e-3 #this becomes unstable if you go any lower
for tt = 1:20
    # forward Euler timestep
    global c += L*c*Δt
    println("Sum c = " *string(sum(c)) *", initial c = " *string(sum(c0)))
end

#make_boundary_conditions: can be "fixed" or "varying"
bc = "fixed" #"fixed" or "varying"
u0 = c0
du = similar(u0)
tspan = (0.0, 50.0)

#Solving differential equation
#NOTE: for DifferentialEquations.jl to work, follow naming conventions in docs
if bc == "fixed" 
    f(du,u,p,t) = mul!(du, L, u) #this avoids allocating a new array for each iteration

elseif bc == "varying"
    #example conditions
    #at t = 0, surface = 1. at t = 50, surface = 0.
    #i.e., surface values linearly decrease from 1 to 0 over the span of 50 years 
    tsfc = [0, 50]
    Csfc = zeros((2, length(dsfc)))
    Csfc[1, :] .= 1
    τ = 1 / 12 #monthly restoring timescale
    f(du, u, p, t) = varying(du, u, p, t, tsfc, Csfc, γ, τ, L, B)
else
    println("invalid boundary condition, must be 'fixed' or 'varying'")
end

#Solve diff eq 
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan)
println("Solving ode")
#solve using QNDF alg - tested against other alg and works fastest 
@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,calck=false)
println("ode solved")

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t), 90,45,33))

#stability check
stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
println("stable: " *string(stable))
println("Gain = "*string(NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])))
    
#everything in sol.u is size 74064 which is all points within the ocean
#we want to make a global map, so we need to use γ to get it back to a matrix
for i in 1:length(sol.t)
    sol_array[i, :, :, :] = vec2fld(sol.u[i], γ.I)
end 

#____PLOTTING____
#time plot
figure()
title("Time distribution")
plot(1:length(sol.t), sol.t, ".")
xlabel("time index [integer]")
ylabel("time [yrs]")

#longitudinal plots
lon_index = 85
dyeplot(γ.lat, γ.depth, sol_array[end, lon_index, :, :]', 0:0.05:1.05)
