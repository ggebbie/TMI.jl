#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1. make a structure which contains m.south, m.north, m.east, m.west, m.up, m.down, m.remote
2. get masks right for each component? (Or compute on fly, just keep tracer mask.)
3. Read "A" 2014 -> construct m::MassFraction
4. Evaluate m*c and see if it equals zero.
5. Make "C" matrix
6. Invert for "m" given C_obs? How about multiple C_obs?
7. Let m have the form of a logistic curve (bounded by 0 and 1)
8. Define u_m control vector which again uses the logistic curve.
9. u_m control part of cost function.
10. Adjoint of TMI core equation wrt m.
11. Optimize for m.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

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

@time m̃ = massfractions(y);
Ã = watermassmatrix(m̃, γ)

# compare m̃ and m_true (Δ̄ = 1e-14 in my case)
for nn in keys(m̃)
    Δ = m̃[nn]
    println(maximum(m̃[nn]))
    println(minimum(m̃[nn]))
    
end

# compare against the "truth" as solved by TMI
m_true = (north = TMI.massfractions_north(A,γ),
    east   = TMI.massfractions_east(A,γ),
    south  = TMI.massfractions_south(A,γ),
    west   = TMI.massfractions_west(A,γ),
    up     = TMI.massfractions_up(A,γ),
    down   = TMI.massfractions_down(A,γ))

# compare m̃ and m_true (Δ̄ = 1e-14 in my case)
for nn in keys(m_true)
   Δ = m_true[nn].fraction - m̃[nn].fraction
   println(maximum(Δ[m_true[nn].γ.wet]))
   println(minimum(Δ[m_true[nn].γ.wet]))
end

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(y.θ)

## reconstruct temperature
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ã,bθ,γ)

# compare to c.θ
Base.maximum(y.θ - θ̃)
Base.minimum(y.θ - θ̃)

## Plot a plan view
# view the surface
cntrs = -2:2:16

# what model depth level?
level = 27
depth = γ.depth[level]

# Help: needs work with continents and labels
GeoPythonPlot.pygui(true) # to help plots appear on screen using Python GUI
GeoPythonPlot.pygui() # to help plots appear on screen using Python GUI

!isdir(TMI.pkgplotsdir()) && mkpath(TMI.pkgplotsdir())
label1 = "θ reconstructed, depth = "*string(depth)*" m"
planviewplot(θ̃, depth, cntrs,
    titlelabel=label1,
    fname = TMI.pkgplotsdir("T_reconstructed_"*string(depth)*".pdf"))

label2 = "θ original, depth = "*string(depth)*" m"
planviewplot(y.θ, depth, cntrs, titlelabel=label2)

colorlabels = [-.3, -.2, -.1, -0.05, -0.02, 0.02, 0.05, .1, .2, .3]
label3 = "θ error [°C], depth = "*string(depth)*" m"
planviewplot(y.θ-θ̃, depth, colorlabels,
    titlelabel=label3,
    fname = TMI.pkgplotsdir("T_error_"*string(depth)*".pdf"))

#### water mass analysis
list = TMI.regionlist()
region = list[4] # sample region

# do numerical analysis
g̃ = watermassdistribution(TMIversion,Ãlu,region,γ);
colorlabels = 0:10:100
glabel = "NATL percentage"
planviewplot(100g̃, depth, colorlabels, titlelabel=label3, fname = TMI.pkgplotsdir("gNATL_"*string(depth)*".pdf"))

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330 # only works if exact
lims = vcat(0,5,10:10:60,100)
tlabel = region * " water-mass [%], "*string(360-lon_section)*"°W"
sectionplot(100g̃,lon_section,lims,titlelabel = tlabel,fname=TMI.pkgplotsdir("gNATL_"*string(360-lon_section)*"W.pdf")) # a

lon_section = 330 # only works if exact
colorlabels = vcat(-.3,-0.1,-0.05:0.01:0.05,0.1,.3)
#tlabel = region * " water-mass fraction [%]"
label3 = "θ error [°C], 30°W"
sectionplot(θ̃-y.θ,lon_section,colorlabels,titlelabel = label3, fname = "T_error_30W.pdf") # a


###################################################
#### EXTRA CODE SNIPPETS ##########################
sum(m̃.down.γ.wet) +
sum(m̃.up.γ.wet) +
sum(m̃.west.γ.wet) +
sum(m̃.east.γ.wet) +
sum(m̃.south.γ.wet) +
sum(m̃.north.γ.wet)

sum(m_true.down.γ.wet) +
sum(m_true.up.γ.wet) +
sum(m_true.west.γ.wet) +
sum(m_true.east.γ.wet) +
sum(m_true.south.γ.wet) +
sum(m_true.north.γ.wet)



mc_test = TMI.tracer_contribution(yθ,m_true.north) # a test
mc_true = TMI.tracer_contribution(yθ,m_true)
maximum(mc_true) < 1e-10
minimum(mc_true) > -1e-10


# alternatively push to Julia backend (VS Code)
# GeoPythonPlot.display(GeoPythonPlot.gcf())

## Plot a lat-depth section
lon_section = 330; # only works if exact
lims = 0:0.1:3.0
sectionplot(PO₄total,lon_section,lims)

## oxygen distribution, just be sure it runs
yO₂ = readfield(TMIfile,"O₂",γ)
bO₂ = getsurfaceboundary(yO₂)
O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)

# Plan view of oxygen at same depth as phosphate
# but different contours
cntrs = 0:20:400 # μmol/kg
label = "oxygen [μmol/kg], depth = "*string(depth)*" m"
planviewplot(O₂, depth, cntrs, titlelabel=label) 

# Section view of oxygen on same phosphate section.
sectionplot(O₂,lon_section,cntrs)