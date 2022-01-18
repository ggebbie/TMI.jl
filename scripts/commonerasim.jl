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
years, bc = readopt(bc_file,γ)
dsfc = bc[begin,:, :, 1][γ.wet[:,:,1]] 

#following make_initial_conditions.m
#initial conditions are the surface patch = 1, propagated down to the bottom of the mixed layer, which we get from B 
#c0 = B * dsfc 
#u0 = c0

#make initial conditions be output of ex1
c0 = tracerinit(γ.wet)
c0[γ.wet] = Alu\bc[begin, :, :, :][γ.wet]
u0 = c0[γ.wet]

tspan = (years[begin], years[end])#tspan must occur within tsfc 
#define varying boundary conditions 
tsfc = years

Csfc = zeros((length(years), length(dsfc)))
for i in 1:length(years)
    Csfc[i, :] = bc[i,:,:, 1][γ.wet[:, :, 1]]
end

τ = 1 / 12 #monthly restoring timescale
li= LinearInterpolation(tsfc, 1:length(tsfc))

#case where u and du are of type Vector{Float64}
du = similar(u0)
LC = DiffEqBase.dualcache(similar(du))
BF = DiffEqBase.dualcache(similar(du))
Cb = similar(Csfc[1,:])
γwet = @view γ.wet[:,:,1]
p = (Csfc,γwet,γ.I,τ,L,B,li,LC,BF,Cb) #parameters - "easier for compiler to handle local variables"
f(du, u, p, t) = varying!(du, u, p, t)

#Solve diff eq 
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan,p)
println("Solving ODE")
#solve using QNDF alg - tested against other alg and works fastest 
@time sol = solve(prob,QNDF(),abstol = 1e-2,reltol=1e-2,saveat=years)
println("ODE solved")
x = Vector{ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEq.OrdinaryDiffEqTag, Float64}, Float64, 12}}(undef, 1000)

x#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))

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
lon_index = 5
dyeplot(γ.lat, γ.depth, sol_array[end, lon_index, :, :]',NaNMath.minimum(sol_array[end, lon_index, :, :]):0.01:NaNMath.maximum(sol_array[end,lon_index, :, :]))

dyeplot(γ.lat, γ.depth, bc[5, lon_index, :, :]' ,NaNMath.minimum(sol_array[end, lon_index, :, :]):0.01:NaNMath.maximum(sol_array[end,lon_index, :, :]))

#surface plots at 3 depths 
figure()
suptitle("values at time = "*string(sol.t[end]))
subplot(3,1,1)
contourf(γ.lon, γ.lat, sol_array[end,:,:,1]',NaNMath.minimum(sol_array[end, :, :, :]):0.01:NaNMath.maximum(sol_array[end,:, :, :]),cmap="coolwarm")
title("depth = "*string( γ.depth[1]))
subplot(3,1,2)
contourf(γ.lon, γ.lat, sol_array[end,:,:,15]',NaNMath.minimum(sol_array[end, :, :, :]):0.01:NaNMath.maximum(sol_array[end,:, :, :]),cmap="coolwarm")
title("depth = "*string( γ.depth[15]))
subplot(3,1,3)
title("depth = "*string(γ.depth[25]))
cf = contourf(γ.lon, γ.lat, sol_array[end,:,:,25]',NaNMath.minimum(sol_array[end, :, :, :]):0.01:NaNMath.maximum(sol_array[end,:, :, :]),cmap="coolwarm")
colorbar(cf,orientation="horizontal")

figure()
contourf(γ.lon, γ.lat, bc[10,:,:,1]',NaNMath.minimum(sol_array[end, :, :, :]):0.01:NaNMath.maximum(sol_array[end,:, :, :]),cmap="coolwarm")

