#= 
nSimulate CE circulation 
Based on ex8.transientsimulation.jl 
Data Input (both files in /data folder) 
  Uses 2° TMI grid
  Theta_anom_OPT-00015.nc downloaded from https://www.ncei.noaa.gov/pub/data/paleo/gcmoutput/gebbie2019/
=#

using Revise, TMI, Interpolations, PyPlot, NaNMath, DifferentialEquations, LinearAlgebra, ForwardDiff,OrdinaryDiffEq, PreallocationTools

#load TMI data
TMIversion = "modern_180x90x33_GH11_GH12"
bc_file = "/home/brynn/Code/TMI.jl/data/Theta_anom_OPT-0015.nc"
@time A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#load boundary condition data
tsfc, bc = readopt(bc_file,γ)
dsfc = bc[begin,:, :, 1][γ.wet[:,:,1]] 

#make initial conditions be output of ex1
#this ends up being silly - initial surface cond are 0 for anom file, so it just ends up being zero everywhere, but I left it in so I  know how to do it in the future 
c0 = tracerinit(γ.wet)
c0[γ.wet] = Alu\bc[begin, :, :, :][γ.wet]
u0 = c0[γ.wet]

#Timespan that diffeq solver will solve for, must be within tsfc 
tspan = (tsfc[begin], tsfc[2])

#Get surface boundary conditions from Theta_anom_OPT 
Csfc = zeros((length(tsfc), length(dsfc)))
for i in 1:length(tsfc)
    Csfc[i, :] = bc[i,:,:, 1][γ.wet[:, :, 1]]
end

#more parameters for diffeq solver 
τ = 1 / 12 #monthly restoring timescale
li = LinearInterpolation(tsfc, 1:length(tsfc))

#Instantiate arrays that the diffeq solver will reallocate
LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
Cb = similar(Csfc[1,:])

#this is silly - but I just need to know which in γ.I are surface layer 
surface_ind = []
for i in 1:size(γ.I)[1]
    if γ.I[i][3] == 1
        push!(surface_ind, i)
    end
end

p = (Csfc,surface_ind,τ,L,B,li,LC,BF,Cb) #parameters

f(du, u, p, t) = varying!(du, u, p, t) #diffeq function to solve 

#Solve diff eq 
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan, p)
println("Solving ODE")
#solve using QNDF alg - tested against other alg and works fastest
#Solver will print out what time step it is on 
@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,saveat=tsfc)
println("ODE solved")

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))

#everything in sol.u is size 74064 which is all points within the ocean
#we want to make a global map, so we need to use γ to get it back to a matrix
for i in 1:length(sol.t)
    sol_array[i, :, :, :] = vec2fld(sol.u[i], γ.I)
end

#stability check
stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
println("stable: " *string(stable))
println("Gain = "*string(NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])))

#____PLOTTING____
#longitudinal plots
lon_index = 5
t_index = 2
absmax = NaNMath.maximum([NaNMath.maximum(abs.(sol_array)), NaNMath.maximum(abs.(sol_array))])
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
