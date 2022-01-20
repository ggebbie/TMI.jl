#= 
Simulate CE circulation 
Based on varying case from ex8.transientsimulation.jl 
Data Input (both files in /data folder) 
  2° TMI grid
  Theta_anom_OPT-00015.nc downloaded from https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/gebbie2019/

Run instructions
 (a) Change bc_file variable to point to Theta_anom file on your machine
 (b) Put the 2degree TMI file into the "data" folder in TMI.jl folder
 (c) Change tspan to cover the period of interest
     If this covers whole time period, it will take ~2 hrs to run
     Change to "tspan = (tsfc[begin],tsfc[end])" to run a quick test
 (d) Output will save to a .nc file in data called ces_output.nc
     Make sure you rename or delete between runs so that it can write  
=#

using Revise, TMI, Interpolations, PyPlot, NaNMath, DifferentialEquations
using LinearAlgebra, ForwardDiff,OrdinaryDiffEq, PreallocationTools

#load TMI data
TMIversion = "modern_180x90x33_GH11_GH12"
bc_file = "/home/brynn/Code/TMI.jl/data/Theta_anom_OPT-0015.nc"
@time A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#load boundary condition data
tsfc, bc = readopt(bc_file,γ)
dsfc = bc[begin,:, :, 1][γ.wet[:,:,1]] 

#make initial conditions be output of ex1
#initial surface cond are 0 for anom file, so it just ends up being 0 everywhere  
c0 = tracerinit(γ.wet)
c0[γ.wet] = Alu\bc[begin, :, :, :][γ.wet]
u0 = c0[γ.wet]

#Timespan that diffeq solver will solve for, must be within tsfc 
tspan = (tsfc[begin], tsfc[2])

#Get surface boundary conditions from Theta_anom_OPT 
Csfc = zeros((length(tsfc), length(dsfc)))
[Csfc[i,:] = bc[i,:,:,1][γ.wet[:,:,1]] for i ∈ 1:length(tsfc)]

#more parameters for diffeq solver 
τ = 1 / 12 #monthly restoring timescale
li = LinearInterpolation(tsfc, 1:length(tsfc))

#Instantiate arrays that the diffeq solver will reallocate
LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
Cb = similar(Csfc[1,:])
surface_ind = findall(x-> x[3] == 1, γ.I) #Find which points in γ.I are on the surface

p = (Csfc,surface_ind,τ,L,B,li,LC,BF,Cb) #parameters
f(du, u, p, t) = varying!(du, u, p, t) #diffeq function to solve 
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan, p)

#solve using QNDF alg - tested against other alg and works fastest
#Solver will print out what time step it is on
println("Solving ODE")
@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,saveat=tsfc)
println("ODE solved")

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))
[sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

#stability check - no interior point should be less/greater than the surface extrema
stable = true ? NaNMath.maximum(sol_array) < NaNMath.maximum(Csfc) && NaNMath.minimum(sol_array) > NaNMath.minimum(Csfc) : false
println("stable: " *string(stable))

#____PLOTTING____
lon_index = 5
t_index = 2
absmax = NaNMath.maximum([NaNMath.maximum(abs.(sol_array)), NaNMath.maximum(abs.(sol_array))])
#longitudinal plots
dyeplot(γ.lat, γ.depth, sol_array[t_index, lon_index, :, :]',-absmax:0.05:absmax)

#surface plots at 3 depths 
figure()
suptitle("values at time = "*string(sol.t[t_index]))

subplot(3,1,1)
contourf(γ.lon, γ.lat, sol_array[t_index,:,:,1]',-absmax:0.05:absmax,cmap="coolwarm")
title("depth = "*string( γ.depth[1]))

subplot(3,1,2)
contourf(γ.lon, γ.lat, sol_array[t_index,:,:,15]',-absmax:0.05:absmax,cmap="coolwarm")
title("depth = "*string( γ.depth[15]))

subplot(3,1,3)
title("depth = "*string(γ.depth[25]))
cf = contourf(γ.lon, γ.lat, sol_array[t_index,:,:,25]',-absmax:0.05:absmax,cmap="coolwarm")
colorbar(cf,orientation="horizontal")

#write output to ncfile 
ces_ncwrite(γ,sol.t,sol_array)
