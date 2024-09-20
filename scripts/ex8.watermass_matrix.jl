#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Given a set of five conservative, perfectly-observed tracers, invert for
the water-mass matrix, A. 

Algorithm: local water-mass inversion from Gebbie & Huybers (2010), J. Phys. Oceanogr.

Details of algorithm:
1. first guess of isotropic mixing and circulation
2. Use least-squares to find smallest deviation that perfectly fits tracers and conserves mass.
3. If non-negativity is violated, use quadratic programming to invert locally.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate("..")

using Revise
using LinearAlgebra
using TMI
using GeoPythonPlot # will load optional extension

TMIversion = versionlist()[6] # G14 has no remote mass fractions
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

# get observations at surface
# set them as surface boundary condition
y = (θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
    P★ = preformedphosphate(TMIversion,Alu,γ),
    δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

w = (θ =  0.01,
    S = 0.001,
    δ¹⁸O = 0.05,
    P★ = 0.05,
    δ¹³C★ = 0.05
)

@time m̃ = massfractions(y,w);
Ã = watermassmatrix(m̃, γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(y.θ)

## reconstruct temperature
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ãlu,bθ,γ)

# compare to c.θ
Base.maximum(y.θ - θ̃)
Base.minimum(y.θ - θ̃)

## Plot a plan view
# view the surface
cntrs = -2:2:16

# what model depth level?
level = 15
depth = γ.depth[level]

# Set up the figure display
#GeoPythonPlot.pygui(true) # to help plots appear on screen using Python GUI
#GeoPythonPlot.pygui() # to help plots appear on screen using Python GUI
# alternatively push to Julia backend (VS Code)
GeoPythonPlot.display(GeoPythonPlot.gcf())

!isdir(TMI.pkgplotsdir()) && mkpath(TMI.pkgplotsdir())
label1 = "θ reconstructed, depth = "*string(depth)*" m"
planviewplot(θ̃, depth, cntrs,
    titlelabel=label1,
    fname = TMI.pkgplotsdir("T_reconstructed_"*string(depth)*".pdf"))

label2 = "θ original, depth = "*string(depth)*" m"
planviewplot(y.θ, depth, cntrs,
    titlelabel=label2,
    fname = TMI.pkgplotsdir("T_obs_"*string(depth)*".pdf"))

colorlabels = [-.3, -.2, -.1, -0.05, -0.02, 0.02, 0.05, .1, .2, .3]
label3 = "θ error [°C], depth = "*string(depth)*" m"
planviewplot(y.θ-θ̃, depth, colorlabels,
    titlelabel=label3,
    fname = TMI.pkgplotsdir("T_error_"*string(depth)*".pdf"))

#### water mass analysis
list = TMI.regionlist()
region = list[4] # sample region

# what model depth level?
level = 26
depth = γ.depth[level]

# do numerical analysis
g̃ = watermassdistribution(TMIversion,Ãlu,region,γ);
g = watermassdistribution(TMIversion,Alu,region,γ);
colorlabels = 0:10:100
glabel = "NATL percentage, depth = "*string(depth)*" m"

planviewplot(100g̃, depth, colorlabels,
    titlelabel=glabel,
    fname = TMI.pkgplotsdir("gNATL_reconstructed_"*string(depth)*"m.pdf")
)

planviewplot(100g, depth, colorlabels,
    titlelabel=glabel,
    fname = TMI.pkgplotsdir("gNATL_obs_"*string(depth)*"m.pdf")
)

dcolorlabels = -10:10
dglabel = "NATL difference %, depth = "*string(depth)*" m"
planviewplot(100(g-g̃), depth, dcolorlabels,
    titlelabel=glabel,
    fname = TMI.pkgplotsdir("gNATL_difference_"*string(depth)*"m.pdf")
)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330 # only works if exact
lims = vcat(0,5,10:10:60,100)
tlabel = region * " water-mass [%], "*string(360-lon_section)*"°W"
sectionplot(100g̃,lon_section,lims,titlelabel = tlabel,fname=TMI.pkgplotsdir("gNATL_"*string(360-lon_section)*"W.pdf")) # a

lon_section = 330 # only works if matches grid exactly
colorlabels = vcat(-.3,-0.1,-0.05:0.01:0.05,0.1,.3)
#tlabel = region * " water-mass fraction [%]"
label3 = "θ error [°C], 30°W"
sectionplot(θ̃-y.θ,lon_section,colorlabels,
    titlelabel = label3,
    fname = TMI.pkgplotsdir("T_error_"*string(360-lon_section)*"W.pdf")
)

# compare against the "truth" as solved by TMI
m_true = (north = TMI.massfractions_north(A,γ),
    east   = TMI.massfractions_east(A,γ),
    south  = TMI.massfractions_south(A,γ),
    west   = TMI.massfractions_west(A,γ),
    up     = TMI.massfractions_up(A,γ),
    down   = TMI.massfractions_down(A,γ))

# compare m̃ and m_true (some big values in locations that are nearly homogeneous)
for nn in keys(m_true)
   Δ = m_true[nn].fraction - m̃[nn].fraction
   println(maximum(Δ[m_true[nn].γ.wet]))
   println(minimum(Δ[m_true[nn].γ.wet]))
end
