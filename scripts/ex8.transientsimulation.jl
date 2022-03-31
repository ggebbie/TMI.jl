#=
 Example 8: Propagate initial conditions through time according to L-matrix 
 Steps: (a) Define a surface patch (d) of concentration 1 to propagate
        (b) Use a fixed euler timestep approximation to estimate how L
            will propagate d
        (c) Solve dc/dt = Lu for the fixed boundary condition case 
 Varying boundary condition: Solve for tracer concentration given an initial
 global surface tracer concentration of 1 that linearly decreases for 50 years
 until global surface tracer concentration is 0. Functionality can be turned on
 by toggling bc from "fixed" to "varying" 
 Solves dc/dt = Lu + Bf 
=#

using Revise, TMI, Interpolations, PyPlot
using NaNMath, DifferentialEquations, LinearAlgebra

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#In MATLAB code: make_initial_condition 
latbox = [50,60]
lonbox = [-50,0]
d = surfacepatch(lonbox, latbox, γ) #vector(same dims as γ, global)  describing surface patch in terms of γ

#flatten d and access at surface wet positions 
dsfc = d.tracer[d.wet]

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
bc = "varying" #"fixed" or "varying"
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
    li = LinearInterpolation(tsfc, 1:length(tsfc))

    #Instantiate arrays that the diffeq solver will reallocate
    LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
    BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
    Cb = similar(Csfc[1,:])
    surface_ind = findall(x-> x[3] == 1, γ.I)#which points in γ.I are on the surface  
    
    p = (Csfc, surface_ind, τ,L,B, li, LC,BF, Cb) #parameters 
    f(du, u, p, t) = varying!(du, u, p, t) #ODEfunction 
else
    println("invalid boundary condition, must be 'fixed' or 'varying'")
end

#Solve diff eq
operator = DiffEqArrayOperator(L)
#isconstant(operator) # not currently working
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan,p)
println("Solving ode")
#solve using QNDF alg - tested against other alg and works fastest 
@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,saveat =tspan[1]:tspan[2])
println("ode solved")

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))
[sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

#stability check
stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
println("stable: " *string(stable))

#____PLOTTING____
#longitudinal plots
#lon_index = 85
#dyeplot(γ.lat, γ.depth, sol_array[25, lon_index, :, :]', 0:0.05:1.05)

lon_section = 330; # only works if exact
lims = 0:0.05:3.0
snapshot = TMI.Field(sol_array[25,:,:,:],γ)
label = "Some quantity [units]"
sectionplot(snapshot,lon_section,lims,titlelabel=label)
