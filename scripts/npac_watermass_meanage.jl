#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate("..")

using Revise
using LinearAlgebra
using TMI
using GeoPythonPlot # will load optional extension
using SparseArrays

TMIversion = versionlist()[1] # G14 has no remote mass fractions
A, Alu, γ, TMIfile, L, B = config(TMIversion);

vol= cellvolume(γ)

## big control volume to stabilize results
xlim = (125,240); ylim = (30,65); zlim = (2000,3000)
#xlim = (209,211); ylim = (34,36); zlim = (2400,2600)
#xlim = (125,240); ylim = (20,40); zlim = (1500,2500)

# test one: get south water-mass fraction
b_south = ones(2,30,γ,:bc_south,"South Boundary","nondim")
b_up = ones(3,25,γ,:bc_upper,"Upper Boundary","nondim")
b_lo = ones(3,29,γ,:bc_lower,"Lower Boundary","nondim")
b_surface = ones(3,1,γ,:bc_surface,"Surface","nondim")

# fix A matrix
# find all locations where boundary set to one
b_ones = (surface = b_surface,
          south = b_south,
          upper = b_up,
          lower = b_lo)

bmask_south = TMI.boundarymask(b_south, γ)
bmask_up = TMI.boundarymask(b_up, γ)
bmask = TMI.boundarymask(b_ones, γ)
maximum(bmask)
sum(bmask)

γb = Grid(b_ones, γ)

# update water-mass matrix to be consistent with boundary points and grid 
Abc = watermassmatrix(A, γb)

# case 1: get m_south 
b_s  = (surface = 0.0*b_surface,
          south = true*b_south,
          upper = false*b_up,
          lower = false*b_lo)
m_south = steadyinversion(lu(Abc),b_s,γb)
maximum(m_south)
minimum(m_south)

## get average over cube
cube = TMI.incube(xlim,ylim,zlim, γ)
function cube_mean(θ::Field, cube, vol, γ)
    # get cube volume averaged temperature and total volume
    cubevec = cube[γ.wet]
    θcube = 0.0
    cumvol = 0.0
    for i in eachindex(cubevec)
        if cubevec[i]
            θcube += θ.tracer[γ.wet][i] * vol[i]
            cumvol += vol[i]
        end
    end
    return θcube /= cumvol
end

# capture fixed variables
cube_mean(θ) = cube_mean(θ, cube, vol, γb)
cube_mean(m_south) # 0.17

## lower boundary
b_l  = (surface = 0.0*b_surface,
          south = 0.0*b_south,
          upper = 0.0*b_up,
          lower = 1.0*b_lo)
m_lower = steadyinversion(lu(Abc),b_l,γb)
maximum(m_lower)
minimum(m_lower)
cube_mean(m_lower) # 0.51

## upper boundary
b_u  = (surface = 0.0*b_surface,
          south = 0.0*b_south,
          upper = 1.0*b_up,
          lower = 0.0*b_lo)
m_upper = steadyinversion(lu(Abc),b_u,γb)
cube_mean(m_upper) # 0.32

##surface boundary (non-local)
b_srf  = (surface = 1.0*b_surface,
          south = 0.0*b_south,
          upper = 0.0*b_up,
          lower = 0.0*b_lo)
m_srf = steadyinversion(lu(Abc),b_srf,γb)
cube_mean(m_srf) # 0.002

maximum(m_lower)
minimum(m_lower)

msum = cube_mean(m_srf) + cube_mean(m_south) + cube_mean(m_lower) + cube_mean(m_upper)

# get the control volume age (i.e., residence time)
setboundarycondition!(amask, a1_boundary)

    F₀ = readfield(TMIfile,"F₀",γ)
    qPO₄ = readsource(TMIfile,"qPO₄",γ) # use this to define mixed-layer

a_boundary = (surface =zeros(3,1,γb,:bc_surface,"Surface","nondim"), 
              south = zeros(2,30,γb,:bc_south,"South Boundary","nondim"),
              upper = zeros(3,25,γb,:bc_upper,"Upper Boundary","nondim"),
              lower = zeros(3,29,γb,:bc_lower,"Lower Boundary","nondim"))

qa = TMI.local_residence_time(TMIfile, Abc, γb)

a_npac = meanage(TMIfile, Abc, a_boundary, γb)
τ_npac = cube_mean( a_npac) # 38.3 years
V_npac = sum(vol.tracer[cube])
F_npac = V_npac / (τ_npac * 3.1e13) # 24.9 Sv (total flux)

minimum(a_npac)
minimum(qa)

maximum(a_npac)
maximum(qa)

# plot a section at 330 east longitude (i.e., 30 west)
GeoPythonPlot.display(GeoPythonPlot.gcf())
lon_section = 202.0 # only works if exact
lon_section = 330 # only works if exact
lims = -5:0
lims = vcat(0,5,10:10:60,100)
tlabel = "NPAC age [yr], "*string(360-lon_section)*"°W"
sectionplot(a_npac,lon_section,lims,titlelabel = tlabel,fname=TMI.pkgplotsdir("aNPAC_"*string(360-lon_section)*"W.pdf")) # a


