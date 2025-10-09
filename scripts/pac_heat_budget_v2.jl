#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using GeoPythonPlot # will load optional extension
using SparseArrays

TMIversion = versionlist()[1] # G14 has no remote mass fractions
A, Alu, γ, TMIfile, L, B = config(TMIversion);

# get observations at surface
# set them as surface boundary condition
y = (θ =  readfield(TMIfile, "θ", γ),
#    S = readfield(TMIfile, "Sp", γ),
    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
)

## get steady version of tracers₄
θb = getsurfaceboundary(y.θ)

## preformed phosphate
θ = steadyinversion(Alu,θb,γ)

## heat budget at one point
Inorthpac = nearestneighbor((210,35,2500), γ)
row_northpac  = γ.R[Inorthpac]
yprime = vec(θ)  .- θ.tracer[Inorthpac]
Tflux = (L[row_northpac,:] .* yprime) * 1e4  # cK/century
col_northpac = Tflux.nzind
γ.I[Tflux.nzind]

sum(L[row_northpac,:] .* vec(θ))

## big control volume to stabilize results
xlim = (125,240); ylim = (30,65); zlim = (2000,3000)
#xlim = (209,211); ylim = (34,36); zlim = (2400,2600)
#xlim = (125,240); ylim = (20,40); zlim = (1500,2500)

# control volume
cube = TMI.incube(xlim,ylim,zlim, γ)
cube_up = TMI.above_cube(xlim,ylim,zlim, γ)
cube_down = TMI.below_cube(xlim,ylim,zlim, γ)
cube_east = TMI.east_of_cube(xlim,ylim,zlim, γ)
cube_west = TMI.west_of_cube(xlim,ylim,zlim, γ)
cube_north = TMI.north_of_cube(xlim,ylim,zlim, γ)
cube_south = TMI.south_of_cube(xlim,ylim,zlim, γ)

τ_vec = -diag(L).nzval  #1/ years
τ = similar(γ.wet, Float64)  # or however you originally constructed τ
τ[γ.wet] .= τ_vec

volume_average(x, mask, volume, γ) = sum(x .* mask .* volume .* γ.wet) / sum(mask .* volume .* γ.wet)

L_pt = -diag(L)[16200, :].nzval

vol= cellvolume( γ)

θcube = volume_average(θ.tracer, cube, vol.tracer, γ)
τcube = volume_average(τ, cube, vol.tracer, γ) #average residence in years
Lcube = A .* τcube

# get A into m3 / s
scaler = Diagonal(vol ./ 3.154e+7 )
Lm3 = scaler * Lcube

# this version gets the local heat budget and then sums itdTdt = 0.0
dTdt = 0.0
dTdt_down = 0.0
dTdt_up = 0.0
dTdt_east = 0.0
dTdt_west = 0.0
dTdt_north = 0.0
dTdt_south = 0.0
dTdt_int = 0.0

cubevec = cube[γ.wet]
cumvol = 0.0
for i in eachindex(cubevec)
    if cubevec[i]
        local ypr = θ.tracer[γ.wet]  .- θcube
        global dTdt += sum(dropzeros(Lm3[i,:] .* ypr))
        global dTdt_down += sum(Lm3[i,:] .* (ypr.*cube_down[γ.wet]))
        global dTdt_up += sum(Lm3[i,:] .* (ypr.*cube_up[γ.wet]))
        global dTdt_east += sum(Lm3[i,:] .* (ypr.*cube_east[γ.wet]))
        global dTdt_west += sum(Lm3[i,:] .* (ypr.*cube_west[γ.wet])) 
        global dTdt_north += sum(Lm3[i,:] .* (ypr.*cube_north[γ.wet]))
        global dTdt_south += sum(Lm3[i,:] .* (ypr.*cube_south[γ.wet]))
        global dTdt_int += sum(Lm3[i,:] .* (ypr.*cube[γ.wet]))
        global cumvol += vol[i]
    end
end

fac = cumvol / (3.154e+7 * 1e4)
dTdt /= fac
dTdt_down /= fac
dTdt_up /= fac
dTdt_east /= fac
dTdt_west /= fac
dTdt_north /= fac
dTdt_south /= fac
dTdt_int /= fac

Inorthpac = nearestneighbor((210,35,2500), γ)
Inorthpac_b = nearestneighbor((210,35,3000), γ)
row_northpac  = γ.R[Inorthpac]
row_northpacb  = γ.R[Inorthpac_b]

(L * θ).tracer[row_northpac, :]
A[row_northpac, :] .* vol.tracer[γ.wet][row_northpac]
A[row_northpacb, :] .* vol.tracer[γ.wet][row_northpacb]

dTdt
dTdt_north + dTdt_south + dTdt_east + dTdt_west + dTdt_up + dTdt_down + dTdt_int

# check mass conservation at one point
Inorthpac = nearestneighbor((210,35,2500), γ)
row_northpac  = γ.R[Inorthpac]
yprime = vec(θ)  .- θ.tracer[Inorthpac]
Tflux = (L[row_northpac,:] .* yprime) * 1e4  # cK/century
col_northpac = Tflux.nzind

Fd = min.(Lsv[row_northpac,:].nzval,Lsv[:,row_northpac].nzval)
Fa = Lsv[row_northpac,:].nzval - Fd 

# split entire matrix into 2 components
Ld = deepcopy(Lsv)
La = deepcopy(Lsv)
I, J, V = findnz(Lsv)

for k in 1:length(Lsv.rowval)
    i = I[k]
    j = J[k]
    if Lsv[i,j] > Lsv[j,i]
        global La[i, j] = Lsv[i,j] - Lsv[j,i]
        global Ld[i,j] = Lsv[j,i]
    elseif Lsv[i,j] ≤ Lsv[j,i]
        global La[i, j] = 0.0
        # Ld[i, j] = Lsv[i, j] # doesn't need any change 
    end
end

# make diff and advection individually conserve mass
for r in 1:size(La,1)
    La[r,r] -= sum(La[r,:])
    Ld[r,r] -= sum(Ld[r,:])
end




## end new work #########################################

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
