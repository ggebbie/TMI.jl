#= 
Simulate CE circulation 
Based on ex8.transientsimulation.jl 
Data Input (both files in /data folder) 
  Uses 2° TMI grid
  Theta_anom_OPT-00015.nc downloaded from https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/gebbie2019/

=#

using Revise, TMI, Interpolations, PyPlot, NaNMath, DifferentialEquations, LinearAlgebra

#load TMI data
TMIversion = "modern_180x90x33_GH11_GH12"
bc_file = "/home/brynn/Code/TMI.jl/data/Theta_anom_OPT-0015.nc"

A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#load boundary condition data
years, bc = readopt(bc_file)
dsfc = bc[end, :, :]'[γ.wet[:,:,1]] 

#following make_initial_conditions.m
#initial conditions are the surface patch = 1, propagated down to the bottom of the mixed layer, which we get from B 
c0 = B * dsfc 
u0 = c0
du = similar(u0)
tspan = (years[begin], years[10])#tspan must occur within tsfc 

#define varying boundary conditions 
tsfc = years

Csfc = zeros((length(years), length(dsfc)))
for i in 1:length(years)
    Csfc[i, :] = bc[i, :, :]'[γ.wet[:, :, 1]]
end

τ = 1 / 12 #monthly restoring timescale
f(du, u, p, t) = varying(du, u, p, t, tsfc, Csfc, γ, τ, L, B)

#Solve diff eq 
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan)
println("Solving ODE")
#solve using QNDF alg - tested against other alg and works fastest 
@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,calck=false)
println("ODE solved")

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
