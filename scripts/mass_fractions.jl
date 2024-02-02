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
import Pkg; Pkg.activate("..")

using Revise
using LinearAlgebra
using TMI
using SparseArrays
#using GeoPythonPlot # will load optional extension

#TMIversion = versionlist()[6] # G14 has no remote mass fractions
#A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);


#CSJ wants to load the MITgcm grid
TMIversion = "simple_MITgcm"


pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

TMIfile = pkgdatadir(""*TMIversion*".nc")
println(TMIfile)
γ = Grid(TMIfile, "hFacC", "XC", "YC", "Z")
θtrue = readfield(TMIfile,"THETA",γ)
u =  Field(-0.01.*ones(size(γ.wet)),γ,θtrue.name,θtrue.longname,θtrue.units)

# get observations at surface
# set them as surface boundary condition
c = (θ =  readfield(TMIfile, "THETA", γ),
    S = readfield(TMIfile, "SALT", γ),
#    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
#    P★ = preformedphosphate(TMIversion,Alu,γ),
#    δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

@time m̃ = TMI.local_solve(c) 

Ã = watermassmatrix(m̃, γ);

# compare against the "truth" as solved by TMI
#m_true = (north = TMI.massfractions_north(A,γ),
#    east   = TMI.massfractions_east(A,γ),
#    south  = TMI.massfractions_south(A,γ),
#    west   = TMI.massfractions_west(A,γ),
#    up     = TMI.massfractions_up(A,γ),
#    down   = TMI.massfractions_down(A,γ))

# compare m̃ and m_true (Δ̄ = 1e-14 in my case)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(c.θ)
uvec = vec(c.θ)
## reconstruct temperature
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ã,bθ,γ)

Adiag = diag(Ã)

Diag_a = spdiagm(Adiag)

A_no_diag = Ã - Diag_a

maxAnodiag = maximum(A_no_diag)
minAnodiag = minimum(A_no_diag)

println("max A no diagonal:$maxAnodiag")
println("min A no diagonal:$minAnodiag")

count_gt_0 = count(x->x>0, A_no_diag)

println("Agt0:$count_gt_0")

maxAdiag = maximum(Adiag)
minAdiag = minimum(Adiag)

println("max A diagonal:$maxAdiag")
println("min A diagonal:$minAdiag")

maxθ̃ = maximum(θ̃)
minθ̃ = minimum(θ̃)

println("max theta:$maxθ̃")
println("min theta:$minθ̃")

surfind = surfaceindex(γ.I)
#northind = findindex(γ.I,latindex,128)
northind = findall(latindex(γ.I) .>= 100)
#println(intersect(surfind,northind))
nadw_ind = intersect(surfind,northind)
nadwsurf = zeros(size(uvec))
nadwsurf[nadw_ind] .= 1

nadw = Ã \ nadwsurf
nadw_to_file = unvec(u,nadw)

writefield("theta_jake.nc",θ̃)
writefield("nadw_jake.nc",nadw_to_file)
# compare to c.θ

## Plot a plan view
# view the surface
#cntrs = -2:2:16

# what model depth level?
#level = 15
#depth = γ.depth[level]

# Help: needs work with continents and labels
#GeoPythonPlot.pygui(true) # to help plots appear on screen using Python GUI

#label1 = "θ reconstructed, depth = "*string(depth)*" m"
#planviewplot(θ̃, depth, cntrs, titlelabel=label1)

#label2 = "θ original, depth = "*string(depth)*" m"
#planviewplot(c.θ, depth, cntrs, titlelabel=label2)

#label3 = "θ difference, depth = "*string(depth)*" m"
#planviewplot(c.θ-θ̃, depth, -.2:.01:.2, titlelabel=label3)

###################################################
#### EXTRA CODE SNIPPETS ##########################
#sum(m̃.down.γ.wet) +
#sum(m̃.up.γ.wet) +
#sum(m̃.west.γ.wet) +
#sum(m̃.east.γ.wet) +
#sum(m̃.south.γ.wet) +
#sum(m̃.north.γ.wet)

#sum(m_true.down.γ.wet) +
#sum(m_true.up.γ.wet) +
#sum(m_true.west.γ.wet) +
#sum(m_true.east.γ.wet) +
#sum(m_true.south.γ.wet) +
#sum(m_true.north.γ.wet)



#mc_test = TMI.tracer_contribution(yθ,m_true.north) # a test
#mc_true = TMI.tracer_contribution(yθ,m_true)
#maximum(mc_true) < 1e-10
#minimum(mc_true) > -1e-10


# alternatively push to Julia backend (VS Code)
# GeoPythonPlot.display(GeoPythonPlot.gcf())

## Plot a lat-depth section
#lon_section = 330; # only works if exact
#lims = 0:0.1:3.0
#sectionplot(PO₄total,lon_section,lims)

## oxygen distribution, just be sure it runs
#yO₂ = readfield(TMIfile,"O₂",γ)
#bO₂ = getsurfaceboundary(yO₂)
#O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)

# Plan view of oxygen at same depth as phosphate
# but different contours
#cntrs = 0:20:400 # μmol/kg
#label = "oxygen [μmol/kg], depth = "*string(depth)*" m"
#planviewplot(O₂, depth, cntrs, titlelabel=label) 

# Section view of oxygen on same phosphate section.
#sectionplot(O₂,lon_section,cntrs)
