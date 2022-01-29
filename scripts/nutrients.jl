#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Given surface boundary conditions and interior sources of phosphate, produce 3D fields of various tracers.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#

using Revise, TMI, GoogleDrive
using PyPlot, PyCall, Test
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
qPO₄ = readtracer(TMIfile,"qPO₄")

# PO₄ inversion
PO₄obs = readtracer(TMIfile,"PO₄")
PO₄sfc = PO₄obs[:,:,1] # nutrient surface boundary condition
r = 1.0
PO₄model = steady_inversion(Alu,PO₄sfc,qPO₄,r,γ.wet)
@test isapprox(PO₄model[50,30,10],PO₄obs[50,30,10])

# NO₃ inversion
NO₃obs = readtracer(TMIfile,"NO₃")
NO₃sfc = NO₃obs[:,:,1] # nutrient surface boundary condition
r = 15.5
NO₃model = steady_inversion(Alu,NO₃sfc,qPO₄,r,γ.wet)
@test isapprox(NO₃model[50,30,10],NO₃obs[50,30,10])

# O₂ inversion
O₂obs = readtracer(TMIfile,"O₂")
O₂sfc = O₂obs[:,:,1] # nutrient surface boundary condition
r = -170.0
O₂model = steady_inversion(Alu,O₂sfc,qPO₄,r,γ.wet)
@test isapprox(O₂model[50,30,10],O₂obs[50,30,10])

# δ¹³C inversion
δ¹³Cobs = readtracer(TMIfile,"δ¹³C")
δ¹³Csfc = δ¹³Cobs[:,:,1] # nutrient surface boundary condition
r = -1.1
δ¹³Cmodel = steady_inversion(Alu,δ¹³Csfc,qPO₄,r,γ.wet)
@test isapprox(δ¹³Cmodel[50,30,10],δ¹³Cobs[50,30,10])
