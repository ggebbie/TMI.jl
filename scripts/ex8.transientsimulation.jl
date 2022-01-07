#= placeholder for transient simulation
#how does a patch of tracer (concentration = 1) evolve throughout time?
=#

using Revise, TMI, Interpolations, PyPlot, NaNMath, BenchmarkTools, DifferentialEquations

#load TMI data
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
#here we will only use
#γ: corresponding lat/lon values for L 
#L: tendency matrix, essentially A (pathway matrix) scaled by residence time (tau) 
#each entry of L tells us how each point in the ocean is connected to its neighbors

#define years to determine global tracer output for
years = range(1,3)
NY = length(years)

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

#make_boundary_conditions - assume fixed boundary conditions (for now) so we don't need to do anything

#Fixed euler timestep approximation
c = c0
Δt = 1e-3 #this becomes unstable if you go any lower
for tt = 1:20
    # forward Euler timestep
    global c += L*c*Δt
    println("Sum c = " *string(sum(c)) *", initial c = " *string(sum(c0)))
end

#Solving differential equation
#NOTE: for DifferentialEquations.jl to work, follow naming conventions
using LinearAlgebra
#f(u,p,t) = L*u
u0 = c0
du = similar(u0)
#mul! cuts time from 6m to 45s (245GiB->6.71GiB) for 0:5
f(du,u,p,t) = mul!(du, L, u) 
tspan = (0.0,5000.0)
func = ODEFunction(f, jac_prototype = L)
prob = ODEProblem(func, u0, tspan)
println("Solving ode")

@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,calck=false)
#solver choice
#TRBDF, Rodas5 ~ 180s
#QNDF ~ 130 seconds (on scale w
println("ode solved")

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t), 90,45,33))

#stability check
stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
println("stable: " *string(stable))
println("Gain = "*string(NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])))
    

#everything in sol.u is size 74064 which is all points within the ocean
#we want to make a global map, so we need to use γ to get it back 
for i in 1:length(sol.t)
    sol_array[i, :, :, :] = vec2fld(sol.u[i], γ.I)
end 

#time plot
figure()
plot(1:length(sol.t), sol.t, ".")
xlabel("time index [integer]")
ylabel("time [?yrs?]")

#surface plots
lev = 15
figure()
subplot(2,1,1)
cf = contourf(γ.lon,γ.lat, sol_array[begin, :, :, lev]', levels = 0:0.05:1)
contour(γ.lon, γ.lat, sol_array[begin, :, :, lev]',levels = 0:0.05:1, colors = "black", linewidths = 1)
colorbar(cf)
subplot(2,1,2)
cf = contourf(γ.lon,γ.lat, sol_array[end, :, :, lev]', levels = 0:0.05:1)
contour(γ.lon, γ.lat, sol_array[end, :, :, lev]', colors = "black", linewidths = 1, levels = 0:0.05:1)
colorbar(cf)

#longitudinal plots
lon_index = 85
figure()
subplot(2,1,1)
cf = contourf(γ.lat, γ.depth, sol_array[begin, lon_index, :, :]', levels = 0:0.05:1)
contour(γ.lat, γ.depth, sol_array[begin, lon_index, :, :]', levels = 0:0.05:1, linewidths = 1, colors = "black")
colorbar(cf)
ylim(maximum(γ.depth), minimum(γ.depth))
subplot(2,1,2)
cf = contourf(γ.lat, γ.depth, sol_array[end, lon_index, :, :]',levels = 0:0.05:1)
contour(γ.lat, γ.depth, sol_array[end, lon_index, :, :]', levels = 0:0.05:1, linewidths = 1, colors = "black")
colorbar(cf)
ylim(maximum(γ.depth), minimum(γ.depth))



