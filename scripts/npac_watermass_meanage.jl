#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate("..")

using Revise
using LinearAlgebra
using TMI
using GeoPythonPlot # will load optional extension
using SparseArrays
using GH19
using NCDatasets

versionno = 2
if versionno == 1
    TMIversion = versionlist()[versionno] # G14 has no remote mass fractions
    A, Alu, γ, TMIfile, L, B = config(TMIversion);
    vol= cellvolume(γ)
 

    # test one: get south water-mass fraction
    b_south = ones(2,30,γ,:bc_south,"South Boundary","nondim")
    b_up = ones(3,25,γ,:bc_upper,"Upper Boundary","nondim")
    b_lo = ones(3,29,γ,:bc_lower,"Lower Boundary","nondim")
    b_surface = ones(3,1,γ,:bc_surface,"Surface","nondim")

elseif versionno == 2
    TMIversion = versionlist()[versionno] # G14 has no remote mass fractions
    A, Alu, γ, TMIfile, L, B = config(TMIversion);
    vol= cellvolume(γ)

    # test one: get south water-mass fraction
    #b_south = ones(2,60,γ,:bc_south,"South Boundary","nondim")
    b_south = ones(2,45,γ,:bc_south,"South Boundary","nondim")
    b_up = ones(3,25,γ,:bc_upper,"Upper Boundary","nondim")
    b_lo = ones(3,29,γ,:bc_lower,"Lower Boundary","nondim")
    b_surface = ones(3,1,γ,:bc_surface,"Surface","nondim")
end

# fix A matrix
# find all locations where boundary set to one
b_ones = (surface = b_surface,
          south = b_south,
          upper = b_up,
          lower = b_lo)

# bmask_south = TMI.boundarymask(b_south, γ)
# bmask_up = TMI.boundarymask(b_up, γ)
# bmask = TMI.boundarymask(b_ones, γ)
# maximum(bmask)
# sum(bmask)

# update grid to be consistent with boundary conditions
γb = Grid(b_ones, γ)

# update water-mass matrix to be consistent with boundary points and grid 
Abc = watermassmatrix(A, γb)
Abc_lu = Abc_lu

## get average over cube
## big control volume to stabilize results
xlim = (125,240); ylim = (30,65); zlim = (2000,3000)
xlim = (125,240); ylim = (0,65); zlim = (2000,3000)
#xlim = (209,211); ylim = (34,36); zlim = (2400,2600)
#xlim = (125,240); ylim = (20,40); zlim = (1500,2500)
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

# case 1: get m_south 
b_s  = (surface = 0.0*b_surface,
          south = true*b_south,
          upper = false*b_up,
          lower = false*b_lo)
m_south = steadyinversion(Abc_lu,b_s,γb)
maximum(m_south)
minimum(m_south)

# capture fixed variables
cube_mean(θ) = cube_mean(θ, cube, vol, γb)
mbar_south = cube_mean(m_south) # 0.17/0.17

## lower boundary
b_l  = (surface = 0.0*b_surface,
          south = 0.0*b_south,
          upper = 0.0*b_up,
          lower = 1.0*b_lo)
m_lower = steadyinversion(Abc_lu,b_l,γb)
maximum(m_lower)
minimum(m_lower)
mbar_lower = cube_mean(m_lower) # 0.51/0.51

## upper boundary 
b_u  = (surface = 0.0*b_surface,
          south = 0.0*b_south,
          upper = 1.0*b_up,
          lower = 0.0*b_lo)
m_upper = steadyinversion(Abc_lu,b_u,γb)
mbar_upper = cube_mean(m_upper) # 0.32/0.32

## surface boundary (non-local)
b_srf  = (surface = 1.0*b_surface,
          south = 0.0*b_south,
          upper = 0.0*b_up,
          lower = 0.0*b_lo)
m_srf = steadyinversion(Abc_lu,b_srf,γb)
mbar_srf = cube_mean(m_srf) # 0.002/0.0003

msum = mbar_srf + mbar_south + mbar_lower + mbar_upper

a_boundary = (surface =zeros(3,1,γb,:bc_surface,"Surface","nondim"), 
              south = zeros(2,30,γb,:bc_south,"South Boundary","nondim"),
              upper = zeros(3,25,γb,:bc_upper,"Upper Boundary","nondim"),
              lower = zeros(3,29,γb,:bc_lower,"Lower Boundary","nondim"))

#qa = TMI.local_residence_time(TMIfile, Abc, γb)

a_npac = meanage(TMIfile, Abc, a_boundary, γb)
τ_npac = cube_mean( a_npac) # 38.3 years/ 33.7 
V_npac = sum(vol.tracer[cube])
F_npac = V_npac / (τ_npac * 3.1e13) # 24.9 Sv (total flux)

minimum(a_npac)
maximum(a_npac)

# plot a section at 330 east longitude (i.e., 30 west)
GeoPythonPlot.display(GeoPythonPlot.gcf())
lon_section = 202.0 # only works if exact
lon_section = 330 # only works if exact
lims = -5:0
lims = vcat(0,5,10:10:60,100)
tlabel = "NPAC age [yr], "*string(360-lon_section)*"°W"
sectionplot(a_npac,lon_section,lims,titlelabel = tlabel,fname=TMI.pkgplotsdir("aNPAC_"*string(360-lon_section)*"W.pdf")) # a

# compute endmembers
θraw = readfield(TMIfile,"θ",γb)

# get steady state potential temperature
θb = getsurfaceboundary(θraw)

## preformed phosphate
θ = steadyinversion(Alu,θb,γ)
θ_cube = cube_mean(θ)

# get effective endmember
θ_south = TMI.effective_endmember(Abc_lu,θ,b_s,γb, control_volume= cube)
θ_upper = TMI.effective_endmember(Abc_lu,θ,b_u,γb, control_volume= cube)  
θ_lower = TMI.effective_endmember(Abc_lu,θ,b_l,γb, control_volume= cube)
θ_srf = TMI.effective_endmember(Abc_lu,θ,b_srf,γb, control_volume= cube)

θprime_south = θ_south - θ_cube
θprime_upper = θ_upper - θ_cube
θprime_lower = θ_lower - θ_cube
θprime_srf = θ_srf - θ_cube

Fθ_south = 10000*mbar_south*θprime_south / τ_npac # cK/century 
Fθ_upper = 10000*mbar_upper*θprime_upper / τ_npac # cK/century 
Fθ_lower = 10000*mbar_lower*θprime_lower / τ_npac # cK/century 
Fθ_srf = 10000*mbar_srf*θprime_srf / τ_npac # cK/century 

# try with GH19 anomaly fields
anomaly = true
gh19filename = download_exp("OPT-0015",anomaly)
ds = NCDataset(gh19filename)
θgh19 = GH19.readfield_snapshot(gh19filename,"theta", 11, γ)

# test what the actual temperature change is
Δθgh19 = 100*cube_mean( GH19.readfield_snapshot(gh19filename,"theta", 1, γ) -
                    GH19.readfield_snapshot(gh19filename,"theta", 21, γ))

θgh19_cube = cube_mean(θgh19) # 0.165
θgh19_south = TMI.effective_endmember(Abc_lu,θgh19,b_s,γb, control_volume= cube) # 0.168
θgh19_upper = TMI.effective_endmember(Abc_lu,θgh19,b_u,γb, control_volume= cube)  # 0.157
θgh19_lower = TMI.effective_endmember(Abc_lu,θgh19,b_l,γb, control_volume= cube)  # 0.157
θgh19_srf = TMI.effective_endmember(Abc_lu,θgh19,b_srf,γb, control_volume= cube) # 0.30

θgh19prime_south = θgh19_south - θgh19_cube # 0.003 
θgh19prime_upper = θgh19_upper - θgh19_cube # -0.007
θgh19prime_lower = θgh19_lower - θgh19_cube # -0.008
θgh19prime_srf = θgh19_srf - θgh19_cube # 0.24

Fθgh19_south = 10000*mbar_south*θgh19prime_south / τ_npac # 0.15 cK/century 
Fθgh19_upper = 10000*mbar_upper*θgh19prime_upper / τ_npac # -0.64 cK/century 
Fθgh19_lower = 10000*mbar_lower*θgh19prime_lower / τ_npac # -1.18 cK/century 
Fθgh19_srf = 10000*mbar_srf*θgh19prime_srf / τ_npac # 0.013 cK/century

Fθgh19_south + Fθgh19_upper + Fθgh19_lower + Fθgh19_srf
# compare to Δθgh19
percent_error = 100* (Δθgh19 - (Fθgh19_south + Fθgh19_upper + Fθgh19_lower + Fθgh19_srf)) /
                (Δθgh19 + (Fθgh19_south + Fθgh19_upper + Fθgh19_lower + Fθgh19_srf)) 

