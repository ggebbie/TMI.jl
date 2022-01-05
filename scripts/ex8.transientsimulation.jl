#= placeholder for transient simulation
#how does a patch of tracer (concentration = 1) evolve throughout time?
=#

using Revise, TMI, Interpolations, DifferentialEquations

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
c = c0 
odefunc(A,t) = A*t
func = ODEFunction(odefunc, jac = L)
prob = ODEProblem(func, c0, (years[begin],years[end]))
sol = solve(prob)

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t), 90,45,33))

#everything in sol.u is size 74064 which is all points within the ocean
#we want to make a global map, so we need to use γ to get it back 
for i in 1:length(sol.t)
    sol_array[i, :, :, :] = vec2fld(sol.u[i], γ.I)
end 

using PyPlot, NaNMath

#surface plots - don't look very different 
figure()
subplot(2,1,1)
contourf(γ.lon,γ.lat, sol_array[begin, :, :, 1]')
subplot(2,1,2)
contourf(γ.lon,γ.lat, sol_array[end, :, :, 1]')

figure()
cf = contourf(γ.lon,γ.lat,sol_array[end, :, :, 1]'.-sol_array[begin, :, :, 1]', cmap = "Reds")
colorbar(cf)
title("END - BEGINNING: Gain") 
 

#longitudinal plots
lon_index = 85
figure()
subplot(2,1,1)
cf = contourf(γ.lat, γ.depth, sol_array[begin, lon_index, :, :]')
colorbar(cf)
ylim(maximum(γ.depth), minimum(γ.depth))
subplot(2,1,2)
cf = contourf(γ.lat, γ.depth, sol_array[end, lon_index, :, :]')
colorbar(cf)
ylim(maximum(γ.depth), minimum(γ.depth))

println("Gain = "*string(NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])))

figure()
cf = contourf(γ.lat, γ.depth,sol_array[end, lon_index, :, :]' .-sol_array[begin, lon_index, :, :]', cmap = "Reds")
colorbar(cf)
title("END - BEGINNING: Gain") 
ylim(maximum(γ.depth), minimum(γ.depth))
