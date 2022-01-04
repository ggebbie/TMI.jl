#= placeholder for transient simulation

=#

#how does a patch of tracer (concentration = 1) evolve throughout time?
using Revise, TMI, Interpolations
#load TMI data
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#define years to determine global tracer output for
years = range(1,3)
NY = length(years)

#make_initial_condition - either rectangle or predefined region
#how is this already done? - I think there's a method for this
#following ex1.trackpathways.jl
#trackpathways->surfacepatch
latbox = [50,60]
lonbox = [-50,0]
d = surfacepatch(lonbox, latbox, γ) #vector describing surface patch in terms of γ

dsfc =  d[:,:,1][γ.wet[:,:,1]]

N = size(A)[1]

#following make_initial_conditions.m
N = size(A)[1]
c0 = zeros(N,1)
c0 = B * dsfc

#need to make indices of c0 equal to 1 where the box is
#make_boundary_conditions - assume fixed boundary conditions (for now) so we don't need to do anything
#
#dcdt = L*c0
c = c0
Δt = 1e-6
for tt = 1:1000
   
    # forward Euler timestep
    c += L*c*Δt
    println(sum(c))
    
end

