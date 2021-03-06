module TMI

using Revise
using LinearAlgebra, SparseArrays, NetCDF, Downloads,
    GoogleDrive, Distances, GibbsSeaWater,  
    PyPlot, PyCall, Distributions, Optim,
    Interpolations, LineSearches, MAT, NCDatasets

export config, config_from_mat, config_from_nc,
    vec2fld, fld2vec, surfaceindex,
    lonindex, latindex, depthindex,
    surfacepatch, 
    layerthickness, cellarea, cellvolume,
    planview, planviewplot,
    section, sectionplot,
    plotextent, tracerinit,
    watermassmatrix, watermassdistribution,
    circulationmatrix, boundarymatrix,
    linearindex, nearestneighbor, updatelinearindex,
    nearestneighbormask, horizontaldistance,
    readtracer, readfield,
    cartesianindex, Γ,
    costfunction_gridded_obs, costfunction_gridded_obs!,
    costfunction_point_obs, costfunction_point_obs!,
    trackpathways, regeneratedphosphate, meanage,
    volumefilled, surfaceorigin, synthetic_observations,
    observe,
    steadyclimatology, steadyinversion,
    interpweights, interpindex,
    wetlocation, iswet,
    control2state, control2state!,
    surfacecontrol2field, surfacecontrol2field!,
    sparsedatamap, config2nc, gridprops,
    matrix_zyx2xyz,
    surface_oxygensaturation, oxygen, location_obs,
    getsurfaceboundary, zerosurfaceboundary,
    onesurfaceboundary, setboundarycondition!,
    adjustboundarycondition, adjustboundarycondition!,
    gsetboundarycondition, setsource!,
    zeros, ones, maximum, minimum, (+), (-), (*), dot,
    Grid, Field, BoundaryCondition, vec, unvec!, unvec

import Base: zeros, ones, maximum, minimum, (\)
import Base: (+), (-), (*), vec
import LinearAlgebra: dot

#Python packages - initialize them to null globally
#const patch = PyNULL()
#const ccrs = PyNULL()

# following example at ClimatePlots.jl
const mpl = PyNULL()
const plt = PyNULL()
const cmocean = PyNULL()
const cartopy = PyNULL()

#Initialize all Python packages - install with conda through Julia
function __init__()
    #copy!(patch, pyimport_conda("matplotlib.patches", "patches"))
    #copy!(ccrs, pyimport_conda("cartopy.crs", "ccrs"))

    # following ClimatePlots.jl
    copy!(mpl, pyimport_conda("matplotlib", "matplotlib", "conda-forge"))
    #copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge"))
    #copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))
    copy!(cartopy, pyimport_conda("cartopy", "cartopy", "conda-forge"))

    println("Python libraries installed")
 end

"""
    struct Grid

    TMI grid with accounting for wet/dry points
"""
struct Grid
    lon::Vector{Float64}
    lat::Vector{Float64}
    depth::Vector{Float64}
    I::Vector{CartesianIndex{3}} # index
    R::Array{Int,3}
#    R::LinearIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}} 
    wet::BitArray{3}
end

"""
    struct Field

    This structure permits the grid to be 
    automatically passed to functions with
    the tracer field.

    This structure assumes the Tracer type to be 
    three-dimensional.
"""
struct Field{T}
    tracer::Array{T,3}
    γ::Grid
end

"""
    struct BoundaryCondition

    a plane defined at `dim=dimval`

# Attributes
    `tracer::Array{T,2}`: values on plane
    `i::Vector{T}`: coordinate values on local x-plane
    `j::Vector{T}`: coordinate values on local y-plane
    `k::T`: fixed coordinate value on local z-plane that defines the Boundary Condition plane
    `dim::Int64`: dimension (1,2, or 3) along which the plane's index is fixed
    `dimval::Int64`: plane defined at dim = dimval where dimval is a 1-based index number
    `wet::BitArray{2}`: ocean mask for boundary condition
"""
struct BoundaryCondition{T}
    tracer::Array{T,2}
    i::Vector{T}
    j::Vector{T}
    k::T
    dim::Int64
    dimval::Int64
    wet::BitArray{2}
end

# Credit to DrWatson.jl for these functions
# Didn't want to add dependency for these small functions
#projectdir() = dirname(Base.active_project())

# find packagedir even if TMI is not the active project
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

pkgsrcdir() = joinpath(pkgdir(),"src")
pkgsrcdir(args...) = joinpath(pkgsrcdir(), args...)

include(pkgsrcdir("TMIconfig.jl"))

""" 
    function trackpathways(TMIversion,latbox,lonbox)
    Track the pathways of a user-defined water mass.
     Steps: (a) define the water mass by a rectangular surface patch dyed with passive tracer concentration of         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.
     See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `latbox`: min and max latitude of box
- `lonbox`: min and max longitude of box
- `γ`: TMI grid
# Output
- `c`: fraction of water from surface source
"""
function trackpathways(Alu,latbox,lonbox,γ)

    b = surfacepatch(lonbox,latbox,γ)
    c = steadyinversion(Alu,b,γ)

    return c
end

""" 
    function watermassdistribution(TMIversion,latbox,lonbox)
    Track the pathways of a user-defined water mass.
     Steps: (a) define the water mass by an oceanographically-relevant surface patch dyed with passive tracer concentration of one
         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.
     See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `region`: name of pre-defined surface region
- `γ`: TMI grid
# Output
- `g`: water-mass fraction
"""
function watermassdistribution(TMIversion,Alu,region,γ)

    
    b = surfaceregion(TMIversion,region,γ)
    g = steadyinversion(Alu,b,γ)

    # # do matrix inversion to get quantity of dyed water throughout ocean:
    # g = tracerinit(γ.wet); # pre-allocate c

    # # make methods that make the "wet" index unnecessary
    # g[γ.wet] = Alu\d[γ.wet] # equivalent but faster than `c = A\d`

    return g
end


""" 
    function regeneratedphosphate(TMIversion,Alu,γ)
    Regenerated (i.e., accumulated, remineralized) phosphate
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `PO₄ᴿ`: regenerated phosphate
"""
function regeneratedphosphate(TMIversion,Alu,γ)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    ## read phosphate source
    qPO₄ = readfield(TMIfile,"qPO₄",γ)

    # zero boundary condition
    b₀ = zerosurfaceboundary(γ)
    PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)

    return PO₄ᴿ
end

""" 
    function meanage(TMIversion,Alu,γ)
    Mean or ideal age
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `a`: mean age [yr]
"""
function meanage(TMIversion,Alu,γ)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    if TMIfile[end-1:end] == "nc"

        #F = ncread(file,"F")
        ## read age source
        F₀ = readfield(TMIfile,"F₀",γ)
        qPO₄ = readfield(TMIfile,"qPO₄",γ) # use this to define mixed-layer

        # better to define a Source type
        Iq = findall(x -> x > 0,qPO₄.tracer)

        qa = zeros(γ)
        qa.tracer[Iq] = 1 ./ F₀.tracer[Iq]
        # zero boundary condition
        b₀ = zerosurfaceboundary(γ)
        a = steadyinversion(Alu,b₀,γ,q=qa)

    else
        
        error("not implemented for mat input file")
    end
        
    return a
end

""" 
    function volumefilled(TMIversion)
    Find the ocean volume that has originated from each surface box.
     This is equivalent to solving a sensitivity problem:
     The total volume is V = vᵀ c , where v is the volume of each box 
     and c is the fraction of volume from a given source which
     satisfies the equation A c = d.                     
     Next, dV/d(d) = A⁻ᵀ v, and dV/d(d) is exactly the volume originating from each source.

     See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011. 
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI.grid
# Output
- `volume`: log10 of global ocean volume filled by a surface region, exists at surface, therefore given BoundaryCondition type
"""
function volumefilled(TMIversion,Alu,γ)::BoundaryCondition

    v = cellvolume(γ)
    area = cellarea(γ)
    
    # effectively take inverse of transpose A matrix.
    dVdd = zeros(γ.wet); # pre-allocate array
    dVdd[γ.wet] = Alu'\v[γ.wet]

    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    #volume = zeros(Float64,length(γ.lon),length(γ.lat))
    volume = zeros(γ.wet[:,:,1])
    # this step could use a function with γ.I argument
    [volume[I[ii][1],I[ii][2]] = dVdd[I[ii]] / area[I[ii][1],I[ii][2]] for ii ∈ eachindex(I) if I[ii][3] == 1]

    volume = log10.(volume)

    ∂V∂b  = BoundaryCondition(volume,γ.lon,γ.lat,γ.depth[1],3,1,γ.wet[:,:,1])
    
    return  ∂V∂b 
end

""" 
    function surfaceorigin(TMIversion,loc)
     Find the surface origin of water for some interior box 
     This is equivalent to solving a sensitivity problem:
     The mass fraction at a location `loc` of interest is 
    `c[loc] = δᵀ c`, where `δ` samples the location of the global mass-fraction variable, c.
    Then the sensitivity of `c[loc]` is: d(c[loc])/d(d) = A⁻ᵀ δ.
    The derivative is solved using the constraint: Ac = d.
    The sensitivity is exactly the mass fraction originating from each source.      
    This problem is mathematically similar to determining how the ocean is filled.
# Arguments
- `loc`: location (lon,lat,depth) of location of interest
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `origin`: surface map of fraction of source water for a given location, log10 of effective depth, in terms of a BoundaryCondition
"""
function surfaceorigin(loc,Alu,γ::Grid)::BoundaryCondition

    #A, Alu, γ = config(TMIversion)
    #ctmp = tracerinit(γ.wet)
    δ = interpweights(loc,γ)
    
    # Find nearest neighbor on grid
    # set δ = 1 at grid cell of interest
    #δ = nearestneighbormask(loc,γ)
    # Note: ctrue[γ.wet]'*δ[γ.wet] returns interpolated value

    dvlocdd = zeros(γ.wet); # pre-allocate c
    dvlocdd[γ.wet] = Alu'\δ[γ.wet]

    # origin is defined at sea surface
    #origin = view(dvlocdd,:,:,1)
    dvlocdd = log10.(dvlocdd[:,:,1])
    origin = BoundaryCondition(dvlocdd,γ.lon,γ.lat,γ.depth[1],3,1,γ.wet[:,:,1])
    
    return origin
end

"""
function steadyclimatology(u₀,fg!,iterations)
     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u₀
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)
    subject to Aỹ = d₀ + Γ u₀.                 
    W⁻ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u₀`:
- `fg!`: compute cost function and gradient in place
- `iterations`: number of optimization iterations
"""
function steadyclimatology(u₀,fg!,iterations)
#function steadyclimatology(u₀,Alu,d₀,y,W⁻,fg!,γ)

    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    #out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(),Optim.Options(show_trace=true, iterations = iterations))

    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(linesearch = LineSearches.BackTracking()),Optim.Options(show_trace=true, iterations = iterations))

    return out    
end

"""
    function sparsedatamap(u₀::Vector{T},Alu,b::BoundaryCondition{T},u::BoundaryCondition{T},y::Vector{T},W⁻,wis,locs,Q⁻,γ::Grid;iterations=10) where T <: Real

     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u₀
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)
    subject to Aỹ = d₀ + Γ u₀.                 
    W⁻ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u₀`:
- `Alu`:
- `b`: first guess of boundary conditions and interior sources
- `y`: observations on 3D grid
- `W⁻`: weighting matrix best chosen as inverse error covariance matrix
- `fg!`: compute cost function and gradient in place
- `γ`: grid
"""
function sparsedatamap(u₀::Vector{T},Alu,b::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Vector{T},W⁻,wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}}},locs,Q⁻,γ::Grid,iterations=10) where {N1, N2, T <: Real}

     fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ)
    
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(linesearch = LineSearches.BackTracking()),Optim.Options(show_trace=true, iterations = iterations))

    return out    
end

"""
    function cartesianindex(wet)
    Read and assemble the grid coordinates
    according to a 3D tracer in x,y,z order
# Arguments
- `wet`: BitArray logical mask for wet points
# Output
- `I`: 3D Cartesian indices
"""
cartesianindex(wet::BitArray{3}) = findall(wet)

"""
    function linearindex(wet)
    Read and assemble the grid coordinates.
# Arguments
- `wet`: 3D mask for wet points
# Output
- `R`: array of linear indices, but not a LinearIndices type
"""
function linearindex(wet)
    R = Array{Int64,3}(undef,size(wet))
    fill!(R,0)
    # R = Array{Union{Int64,Nothing},3}(nothing,size(wet))
    R[wet]=1:sum(wet)
    # R = LinearIndices((1:maximum(it),1:maximum(jt),1:maximum(kt)));
    # R = LinearIndices((it,jt,kt));
    #Rwet = R[γ.wet]
    return R
end
                                                  
"""
    function readtracer(file,tracername)
    Read a tracer field from NetCDF.
# Arguments
- `file`: TMI NetCDF file name
- `tracername`: name of tracer
# Output
- `c`: 3D tracer field
"""
function readtracer(file,tracername)
    c = ncread(file,tracername)
    return c
end

"""
    function readfield(file,tracername,γ)
    Read a tracer field from NetCDF but return it 
    as a Field. 
# Arguments
- `file`: TMI NetCDF file name
- `tracername`: name of tracer
- `γ::Grid`, TMI grid specification
# Output
- `c`::Field
"""
function readfield(file,tracername,γ::Grid)
    tracer = ncread(file,tracername)

    # perform a check of file compatibility
    # with grid
    if sum(isnan.(tracer[γ.wet])) > 0
        println("readfield warning: NaN on grid")
    end
    # check for non NaN or nonzero off grid
    # Need to rethink how to do this.
    # if sum( !iszero(tracer[.!γ.wet])) > 0 
    #     println("readfield warning: nonzero value off grid")
    # end
           
    c = Field(tracer,γ)
    return c
end

"""
    function depthindex(I) 
    Get the k-index (depth level) from the Cartesian index
"""
function depthindex(I)
    T = eltype(I[1])
    k = Vector{T}(undef,length(I))
    [k[n]=I[n][3] for n ∈ eachindex(I)]
    return k
end

"""
    function lonindex(I) 
    Get the i-index (lon index) from the Cartesian index
"""
function lonindex(I)
    T = eltype(I[1])
    i = Vector{T}(undef,length(I))
    [i[n]=I[n][1] for n ∈ eachindex(I)]
    return i
end

"""
    function latindex(I) 
    Get the j-index (latitude index) from the Cartesian index
"""
function latindex(I)
    T = eltype(I[1])
    j = Vector{T}(undef,length(I))
    [j[n]=I[n][2] for n ∈ eachindex(I)]
    return j
end

function findindex(I,indexfunc,indexnumber)
    Ifound = findall(indexfunc(I) .== indexnumber)
    return Ifound
end
    
"""
    function surfaceindex(I) 
    Get the vector-index where depth level == 1 and it is ocean.
"""
surfaceindex(I) = findindex(I,depthindex,1)

"""
    function southindex(I) 
    Get the vector-index on the southern open boundary
"""
southindex(I) = findindex(I,latindex,1)

"""
    function northindex(I) 
    Get the vector index on the northern open boundary
"""
northindex(I) = findindex(I,latindex,maximum(latindex(I)))

"""
    function westindex(I) 
    Get the vector index on the western open boundary
"""
westindex(I) = findindex(I,lonindex,1)

"""
    function eastindex(I) 
    Get the vector index on the northern open boundary
"""
eastindex(I) = findindex(I,lonindex,maximum(lonindex(I)))

"""
    Horizontal area of grid cell
"""
function cellarea(γ)
    dx = zonalgriddist(γ)
    dy = haversine((γ.lon[1],γ.lat[1])
                  ,(γ.lon[1],γ.lat[2]))

    area = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(area,0.0)

    # to calculate area everywhere
    #[area[i,j] = dx[j] * dy for i ∈ eachindex(γ.lon) for j ∈ eachindex(γ.lat)]

    # to calculate sea surface area
    I = γ.I
    [area[I[ii][1],I[ii][2]] = dx[I[ii][2]] * dy for ii ∈ eachindex(I) if I[ii][3] == 1]

    return area
end

"""
    Volume of each grid cell.
"""
function cellvolume(γ)
    dz = layerthickness(γ)
    area = cellarea(γ)
    volume = Array{Float64,3}(undef,length(γ.lon),length(γ.lat),length(γ.depth))
    fill!(volume,0.0)

    # for volume everywhere
    # [volume[i,j,k] = area[i,j] * dz[k] for i ∈ eachindex(γ.lon) for j ∈ eachindex(γ.lat) for k ∈ eachindex(γ.depth)]

    # for ocean volume only
    I = γ.I
    [volume[I[ii]] = area[I[ii][1],I[ii][2]] * dz[I[ii][3]] for ii ∈ eachindex(I)]
    return volume
end

function layerthickness(γ::Grid)
    zface= (γ.depth[1:end-1].+γ.depth[2:end])./2;
    dz = ([zface[1] ; diff(zface); 500]);
    return dz
end

function zonalgriddist(γ::Grid)
    dx = similar(γ.lat)
    for j in eachindex(γ.lat)
        dx[j] = haversine((γ.lon[1],γ.lat[j])
                         ,(γ.lon[2],γ.lat[j]))
    end
    return dx
end

"""
    function surfacepatch
    Make a surface boundary condition
    with a rectangular patch
# Arguments
- `lonbox`: longitudes of box edges
- `latbox`: latitudes of box edges
- `γ`: TMI.grid
# Output
- `d`: vector that describes surface patch
"""
function surfacepatch(lonbox,latbox,γ::Grid)::BoundaryCondition

    # ternary operator to handle longitudinal wraparound
    lonbox[1] ≤ 0 ? lonbox[1] += 360 : nothing
    lonbox[2] ≤ 0 ? lonbox[2] += 360 : nothing

    # preallocate
    patch = zeros(γ.wet[:,:,1])

    # can you add a logical to a Float64? yes, it's 1.0
    [patch[i,j] +=  latbox[1] ≤ γ.lat[j] ≤ latbox[2] && lonbox[1] ≤ γ.lon[i] ≤ lonbox[2] for i in eachindex(γ.lon) for j in eachindex(γ.lat)] 

    # 3,1 to identify surface
    b = BoundaryCondition(patch,γ.lon,γ.lat,γ.depth[1],3,1,γ.wet[:,:,1])
    
    return b
end

"""
    function nearestneighbormask
    Make a 3D tracer field that is 1 at location 
    of nearest neighbor, 0 elsewhere
# Arguments
- `loc`: location in a 3-tuple (lon,lat,depth)
- `γ`: TMI.grid
# Output
- `δ`: nearest neighbor mask 3D field
"""
function nearestneighbormask(loc,γ::Grid,N=1)

    Inn, Rnn = nearestneighbor(loc,γ,N)

    # preallocate
    δ = falses(size(γ.wet))
    #Array{BitArray,3}(undef,size(γ.wet))
    #fill!(δ,zero(Bool))
    δ[Inn] = 1 
    return δ
end

"""
    function nearestneighbor(loc,γ)
    return the Cartesian index and linear index 
    of the nearest N neighbors
# Arguments
- `loc`: 3-tuple of lon,lat,depth location
- `γ`: TMI.grid
# Output
- `Inn`: Cartesian indices of nearest neighbor
#- `Rnn`: linear indices of nearest neighbor, Removed from code
"""
function nearestneighbor(loc,γ,N=1)

    xydist = horizontaldistance(loc[1:2],γ)

    if N==1
        ijdist,ijmin = findmin(xydist[γ.wet[:,:,1]])
        kdist,kmin = findmin(abs.(loc[3] .- γ.depth))
    elseif N > 1
        ijmin = sortperm(xydist[γ.wet[:,:,1]])
        kmin = sortperm(abs.(loc[3] .- γ.depth))
    end

    if N == 1
        Inn = CartesianIndex.(γ.I[ijmin][1],γ.I[ijmin][2],kmin)
    elseif N > 1
        Inn = Vector{CartesianIndex}(undef,N)
        cN2 = ceil(Integer,N/2)
        for ii in 1:cN2
            # translate ijmin into imin, jmin
            Inn[ii] = CartesianIndex.(γ.I[ijmin[ii]][1],γ.I[ijmin[ii]][2],kmin[1])
        end
        for ii in 1:floor(Integer,N/2)
            Inn[cN2+ii] = CartesianIndex.(γ.I[ijmin[ii]][1],γ.I[ijmin[ii]][2],kmin[2])
        end        
    end
    
    return Inn
end

"""
    function horizontaldistance(loc,γ)
    return the Cartesian index and linear index 
    of the nearest N neighbors
# Arguments
- `loc`: 3-tuple of lon,lat,depth location
- `γ`: TMI.grid
# Output
- `hordist`: horizontal distance to nearest tracer grid points
"""
function horizontaldistance(loc,γ::Grid)

    # hordist will have same type as lon,lat,depth
    T = eltype(γ.lon)
    
    # pre-allocate horizontal distance
    hordist = Matrix{T}(undef,length(γ.lon),length(γ.lat))
    # will give NaN with appropriate precision
    fill!(hordist,zero(T)/zero(T))
    
    # calculate haversine horizontal distance on sphere
    [hordist[γ.I[ii]] = haversine((loc[1],loc[2]),                  (γ.lon[γ.I[ii][1]],γ.lat[γ.I[ii][2]]))
       for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
    return hordist
end


"""
function interpindex(loc,γ)
    Weights for linear interpolation.
    The derivative of linear interpolation is needed in sensitivity studies.
    ReverseDiff.jl could find this quantity automatically.
    Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.
# Arguments
- `c`: a temporary tracer field, would be nice to make it unnecessary
- `loc`: (lon,lat,depth) tuple of a location of interest
- `γ`: TMI grid
# Output
- `δ`: weights on a 3D tracer field grid
"""
function interpindex(loc,γ)

    # Handle a grid mismatch.
    loc_on_grid = shiftloc(loc,γ)
    
    # Handle longitudinal periodic condition (i.e., wraparound)
    lon = vcat(copy(γ.lon),γ.lon[1]+360.)
    list = vcat(1:length(γ.lon),1)
    nodes = (lon,γ.lat,γ.depth)

    # eliminate need to pass tracer value
    wis = Interpolations.weightedindexes((Interpolations.value_weights,),((Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))),nodes, loc_on_grid)

    # issue, some of weighted points may be NaNs in tracer field
    # handle this in the Interpolations.jl routines
    # may involve chaging Gridded(Linear()) above
    return wis
end

"""
function  shiftloc(loc)

    sometimes loc longitudes are outside of grid due to different conventions
    assumption: 360° shift is enough to get back on grid
"""
function shiftloc(loc,γ)
    # accounts for a half grid cell of overlap space
    newlon = loc[1]
    westlon = (3/2)*γ.lon[1] - (1/2)*γ.lon[2]
    eastlon =   (3/2)*γ.lon[end] - (1/2)*γ.lon[end-1]
    while newlon <  westlon
        newlon += 360.0
    end
    
    while newlon > eastlon
        newlon -= 360.0
    end

    if newlon > eastlon || newlon < westlon
        error("location not on grid")
    end

    return (newlon, loc[2], loc[3])    
end

"""
function interpweights(loc,γ)
    Weights for linear interpolation.
    The derivative of linear interpolation is needed in sensitivity studies.
    ReverseDiff.jl could find this quantity automatically.
    Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.
# Arguments
- `loc`: (lon,lat,depth) tuple of a location of interest
- `γ`: TMI grid
# Output
- `δ`: weights on a 3D tracer field grid
"""
function interpweights(loc,γ)

    # handle wraparound
    # repeated (unnecessarily?) in interpindex
    lon = vcat(copy(γ.lon),γ.lon[1]+360)
    list = vcat(1:length(γ.lon),1)

    wis = interpindex(loc,γ)

    # translate to weights via
    #http://juliamath.github.io/Interpolations.jl/latest/devdocs/
    δ = zeros(γ.wet)

    # changes in δwrap i=91 are translated back to δ i=1
    δwrap = view(δ,list,:,:)
    for ii = 1:2
        for jj = 1:2
            for kk = 1:2
                δwrap[wis[1].istart+ii-1,wis[2].istart+jj-1,wis[3].istart+kk-1] +=
                wis[1].weights[ii]*wis[2].weights[jj]*wis[3].weights[kk]
            end
        end
    end

    # if some adjacent points are dry, then re-normalize to keep this interpolation as an average.
    # The hope is that the interpolation is stable with this approach, but other side effects are likely.
    # Note that this should be handled earlier, like in the interpindex section. For this reason, there could be an inconsistency in the global map function.
    if iszero(sum(filter(!isnan,δ)))
        δ = nothing
    elseif sum(filter(!isnan,δ)) < 1.0
        δ ./= sum(filter(!isnan,δ))
    end
    
    return δ
end

""" 
    function zeros(γ::Grid)
      initialize tracer field on TMI grid
      using a Field struct and constructor
# Arguments
- `γ`::TMI.Grid
# Output
- `d`::Field,  3d tracer field with NaN on dry points
"""
function zeros(γ::Grid)::Field

    # use depth (could have been lon, lat)
    # to get element type
    T = eltype(γ.depth)
    
    # preallocate
    tracer = Array{T}(undef,size(γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer[γ.wet] .= zero(T)
    tracer[.!γ.wet] .= zero(T)/zero(T) # NaNs with right type

    d = Field(tracer,γ)

    return d
end

""" 
    function zeros(wet,ltype=Float64)
    initialize tracer field on TMI grid
    This version will give an array
# Arguments
- `wet`::BitArray mask of ocean points
- `ltype`:: optional type argument, default=Float64
# Output
- `d`:: 3d tracer field with NaN on dry points
"""
function zeros(wet,ltype=Float64)
    # preallocate
    d = Array{ltype}(undef,size(wet))

    # set ocean to zero, land to NaN
    d[wet] .= zero(ltype)
    d[.!wet] .= zero(ltype)/zero(ltype) # NaNs with right type
    return d
end

"""
    function boundaryconditionatts(dim::Int64,dimval::Int64,γ::Grid)

       Help initialize boundary condition by getting some attributes
"""
function boundaryconditionatts(dim::Int64,dimval::Int64,γ::Grid)

    dimsize = size(γ.wet)
    # dumb way to do it
    if dim == 1
        wet2d = γ.wet[dimval,:,:]
        i = γ.lat
        j = γ.depth
        k = γ.lon[dimval]
    elseif dim == 2
        wet2d = γ.wet[:,dimval,:]
        i = γ.lon
        j = γ.depth
        k = γ.lat[dimval]
    elseif dim == 3
        wet2d = γ.wet[:,:,dimval]
        i = γ.lon
        j = γ.lat
        k = γ.depth[dimval]
    else
        error("boundary condition not implemented in 4+ dimensions")
    end
    return i,j,k,wet2d
    
end

"""
    function zeros(dim::Int64,dimval::Int64,γ::Grid)::BoundaryCondition

       Initialize boundary condition with zeroes
"""
function zeros(dim::Int64,dimval::Int64,γ::Grid)::BoundaryCondition

    i,j,k,wet = boundaryconditionatts(dim,dimval,γ)

    tracer = Array{Float64}(undef,size(wet))
    tracer[wet] .= zero(Float64)
    tracer[.!wet] .= zero(Float64)/zero(Float64)
    b = BoundaryCondition(tracer,i,j,k,dim,dimval,wet)

end

"""
    function ones(dim::Int64,dimval::Int64,γ::Grid)::BoundaryCondition

       Initialize boundary condition with ones
"""
function ones(dim::Int64,dimval::Int64,γ::Grid)::BoundaryCondition

    i,j,k,wet = boundaryconditionatts(dim,dimval,γ)

    tracer = Array{Float64}(undef,size(wet))
    tracer[wet] .= ones(Float64)
    tracer[.!wet] .= zero(Float64)/zero(Float64)
    b = BoundaryCondition(tracer,i,j,k,dim,dimval,wet)

end

"""
   Get boundary condition by extracting from 3D tracer
"""
function getboundarycondition(tracer3d,dim,dimval,γ::Grid)::BoundaryCondition

    dimsize = size(γ.wet)
    # dumb way to do it
    if dim == 1
        wet2d = γ.wet[dimval,:,:]
        tracer2d = tracer3d[dimval,:,:]
        i = γ.lat
        j = γ.depth
        k = γ.lon[dimval]
    elseif dim == 2
        wet2d = γ.wet[:,dimval,:]
        tracer2d = tracer3d[:,dimval,:]
        i = γ.lon
        j = γ.depth
        k = γ.lat[dimval]
    elseif dim == 3
        wet2d = γ.wet[:,:,dimval]
        tracer2d = tracer3d[:,:,dimval]
        i = γ.lon
        j = γ.lat
        k = γ.depth[dimval]
    else
        error("boundary condition not implemented in 4+ dimensions")
    end
    
    b = BoundaryCondition(tracer2d,i,j,k,dim,dimval,wet2d)

end

# Define maximum for Field to not include NaNs
maximum(c::Field) = maximum(c.tracer[c.γ.wet])
minimum(c::Field) = minimum(c.tracer[c.γ.wet])
#mean(x::Field) = mean(x.tracer[x.γ.wet])

# Define max/min for BoundaryCondition
maximum(b::BoundaryCondition) = maximum(b.tracer[b.wet])
minimum(b::BoundaryCondition) = minimum(b.tracer[b.wet])

"""
    `function \\(A,d::Field)::Field`
    Define left division for Fields
    Need two slashes to prevent invalid escape
"""
function \(A,d::Field{T})::Field{T} where T <: Real
    # initialize output
    c = zeros(d.γ)
    c.tracer[c.γ.wet] = A\d.tracer[d.γ.wet]
    return c
end

"""
    `function +(c::BoundaryCondition,d::BoundaryCondition)::BoundaryCondition`
    Define addition for Fields
"""
function +(c::BoundaryCondition{T},d::BoundaryCondition{T})::BoundaryCondition{T} where T <: Real

    if c.wet != d.wet # check conformability
        error("BoundaryCondition's not conformable for addition")
    end
    array = zeros(c.wet)
    e = BoundaryCondition(array,c.i,c.j,c.k,c.dim,c.dimval,c.wet)
    
    # a strange formulation to get
    # return e to be correct
    e.tracer[e.wet] += c.tracer[c.wet]
    e.tracer[e.wet] += d.tracer[d.wet]
    return e
end

"""
    `function +(c::Field,d::Field)::Field`
    Define addition for Fields
"""
function +(c::Field{T},d::Field{T})::Field{T} where T <: Real
    # initialize output
    if c.γ.wet != d.γ.wet # check conformability
        error("Fields not conformable for addition")
    end
    e = zeros(d.γ)

    # a strange formulation to get
    # return e to be correct
    e.tracer[e.γ.wet] += c.tracer[c.γ.wet]
    e.tracer[e.γ.wet] += d.tracer[d.γ.wet]
    return e
end

"""
    `function -(c::Field,d::Field)::Field`
    Define addition for Fields
"""
function -(c::Field{T},d::Field{T})::Field{T} where T <: Real
    # initialize output
    if c.γ.wet !== d.γ.wet # check conformability
        error("Fields not conformable for addition")
    end
    e = zeros(d.γ)
    e.tracer[e.γ.wet] += c.tracer[c.γ.wet]
    e.tracer[e.γ.wet] -= d.tracer[d.γ.wet]
    return e
end

"""
    `function *(C,d::Field)::Field`
    Define scalar or matrix multiplication for fields
"""
function *(C,d::Field{T})::Field{T} where T <: Real

    e = zeros(d.γ)
    e.tracer[e.γ.wet] += C*d.tracer[d.γ.wet]
    return e
end

"""
    `function *(C,d::BoundaryCondition)::BoundaryCondition`
    Define scalar or matrix multiplication for BoundaryCondition`s
"""
function *(C,d::BoundaryCondition{T})::BoundaryCondition{T} where T <: Real
    array = zeros(d.wet)
    e = BoundaryCondition(array,d.i,d.j,d.k,d.dim,d.dimval,d.wet)
    e.tracer[e.wet] += C*d.tracer[d.wet]
    return e
end

"""
    `function *(c::Field,d::Field)::Field`
    Field by field multiplication is element-by-element.
"""
function *(c::Field{T},d::Field{T})::Field{T} where T <: Real

    e = zeros(d.γ)
    e.tracer[e.γ.wet] += c.tracer[c.γ.wet] .* d.tracer[d.γ.wet]
    return e
end


"""
    `function *(c::Field,d::Field)::Field`
    Field by field multiplication is element-by-element.
"""
function dot(c::Field{T},d::Field{T})::T where T <: Real

    e =  c.tracer[c.γ.wet]' * d.tracer[d.γ.wet]
    return e
end


""" 
    function setsource!(d::Field,q::Field,r::Number)
    apply interior source q to the equation constraints d
# Arguments
- `d`::Field, equation constraints (i.e., right hand side)
- `q`::Field, interior source
- `r`::Number, default = 1.0, stoichiometric ratio
"""
function setsource!(d::Field{T},q::Field{T},r=1.0) where T<: Real

    # warning d.γ.wet and q.γ.wet must match
    d.tracer[d.γ.wet] -= r * q.tracer[q.γ.wet]

end

"""
    function vec(u)

    Turn a collection of controls into a vector
    for use with Optim.jl. 
    An in-place version of this function would be handy.
"""
function vec(u::BoundaryCondition{T}) where T <: Real
    uvec = u.tracer[u.wet]
    return uvec
end

"""
    function vec(u)

    Turn a collection of controls into a vector
    for use with Optim.jl. 
    An in-place version of this function would be handy.
"""
function vec(u::NamedTuple{<:Any,NTuple{N,BoundaryCondition{T}}}) where {N, T <: Real}

    uvec = Vector{T}(undef,0)
    for v in u
        append!(uvec,v.tracer[v.wet])
    end
    return uvec
end

"""
    function unvec!(u,uvec)

    Undo the operations by vec(u)
    Needs to update u because attributes of 
    u need to be known at runtime.
"""
function unvec!(u::NamedTuple{<:Any,NTuple{N,BoundaryCondition{T}}},uvec::Vector{T}) where {N, T <: Real}

    counter = 0
    for v in u
        n = sum(v.wet)
        v.tracer[v.wet] = uvec[counter+1:counter+n]
        counter += n
    end
end

"""
    function unvec(utemplate,uvec)

    Undo the operations by vec(u)
    Needs to update u because attributes of 
    u need to be known at runtime.
"""
function unvec(utemplate::NamedTuple{<:Any,NTuple{N,BoundaryCondition{T}}},uvec::Vector{T}) where {N, T <: Real}

    counter = 0
    vals = Vector{BoundaryCondition}(undef,N)
    for (ii,vv) in enumerate(utemplate)
        n = sum(vv.wet)
        vals[ii] = unvec(vv, uvec[counter+1:counter+n])
        counter += n
    end
    u = (;zip(keys(utemplate), vals)...)
    return u
end

"""
    function unvec!(u,uvec)

    Undo the operations by vec(u)
    Needs to update u because attributes of 
    u need to be known at runtime.
"""
function unvec!(u::BoundaryCondition{T},uvec::Vector{T}) where T <: Real

    I = findall(u.wet)
    counter = 0
    for i in I
        counter +=1
    # doesn't work to pass u back
        #u.tracer[u.wet] = uvec
        u.tracer[i] = uvec[counter]
    end
    
end

"""
    function unvec(u,uvec)

    Replace u with new u
    Undo the operations by vec(u)
    Needs to update u because attributes of 
    u need to be known at runtime.
"""
function unvec(utemplate::BoundaryCondition{T},uvec::Vector{T}) where T <: Real

    tracer = zeros(utemplate.wet)
    tracer[utemplate.wet] = uvec
    u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,utemplate.wet)
    return u
end

"""
Surface oxygen saturation value and fraction of saturation value in field 
"""
function surface_oxygensaturation(file)
    # read temperature and o2.
    θ = readtracer(file,"θ")
    θsurface = view(θ,:,:,1)

    S = readtracer(file,"Sp")
    Ssurface = view(S,:,:,1)

    # GibbsSeaWater.jl for saturation value
    O₂sol = gsw_o2sol_sp_pt.(Ssurface, θsurface)

    O₂ = readtracer(file,"O₂")
    O₂surface = view(O₂,:,:,1)
    O₂fraction = O₂surface./O₂sol

    return O₂sol, O₂fraction 
end

"""
Reconstruct dissolved oxygen (that doesn't exist in TMI product)
by assuming same oxygen saturation fraction as modern
"""
function oxygen(version,O₂fraction)

    A, Alu, γ, file = config_from_nc(version)

    o2po4ratio = 170
    
    # read temperature and o2.
    θ = readtracer(file,"θ")
    θsurface = view(θ,:,:,1)

    S = readtracer(file,"Sp")
    Ssurface = view(S,:,:,1)

    # GibbsSeaWater.jl for saturation value
    O₂sol = gsw_o2sol_sp_pt.(Ssurface, θsurface)

    O₂surface = O₂sol.*O₂fraction

    # invert with stoichiometric ratio
    O₂ = tracerinit(γ.wet)
    qPO₄ = readtracer(file,"qPO₄")
    d = o2po4ratio*qPO₄

    #d = qO₂lgm
    d[:,:,1] = O₂surface;

    O₂[γ.wet] =  Alu\d[γ.wet]

    return O₂
end


""" 
    function synthetic_observations(TMIversion,variable)
    Synthetic observations that are a contaminated version of real observations
    This version: gridded observations
# Arguments
- `TMIversion::String`: version of TMI water-mass/circulation model
- `variable::String`: variable name to use as template
# Output
- `y`: contaminated observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `θtrue`: real observations, 3D field
"""
function synthetic_observations(TMIversion,variable,γ)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    θtrue = readfield(TMIfile,variable,γ)
    σθ = readfield(TMIfile,"σ"*variable,γ)

    #ntrue = zeros(γ)
    #ntrue = zeros(γ.wet)
    #ntrue += rand(Normal(),length(σθ[γ.wet])) .* σθ[γ.wet]
    ntrue = Field(rand(Normal(),size(γ.wet)),γ)
    ntrue *= σθ

    y = θtrue + ntrue
    
    #y = θtrue .+ ntrue

    # get cost function (J) based on model misfit
    # here the data-model misfit is weighted by the expected error

    # weighting matrix
    #W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    W⁻ = (1/sum(γ.wet)) .* Diagonal(1 ./σθ.tracer[γ.wet].^2)
    return y, W⁻, θtrue
end
 
""" 
    function synthetic_observations(TMIversion,variable,locs)
    Synthetic observations that are a contaminated version of real observations
    This version: observations with random (uniform) spatial sampling
# Arguments
- `TMIversion::String`: version of TMI water-mass/circulation model
- `variable::String`: variable name to use as template
- `N`: number of observations
# Output
- `y`: contaminated observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `ytrue`: uncontaminated observations, 3D field
- `locs`: 3-tuples of locations for observations
- `wis`: weighted indices for interpolation to locs sites
"""
function synthetic_observations(TMIversion,variable,γ,N,σ=nothing)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    
    θtrue = readfield(TMIfile,variable,γ)
    replace!(θtrue.tracer,NaN=>0.0)

    if isnothing(σ)
        σθ = readfield(TMIfile,"σ"*variable,γ)
        replace!(σθ.tracer,NaN=>0.0)
    end

    # get random locations that are wet (ocean)
    locs = Vector{Tuple{Float64,Float64,Float64}}(undef,N)
    [locs[i] = wetlocation(γ) for i in eachindex(locs)]

    # get weighted interpolation indices
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    ytrue = observe(θtrue,wis,γ)
    if isnothing(σ)
        σtrue = observe(σθ,wis,γ)
    else
        σtrue = σ * ones(N)
    end

    #ntrue = rand(Normal.(zeros(N),σtrue),N)# .* σtrue
    ntrue = rand(Normal(),N).*σtrue
    y = ytrue .+ ntrue

    # weighting matrix
    #W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    W⁻ = (1/N) .* Diagonal(1 ./σtrue.^2)
    return y, W⁻, θtrue, ytrue, locs, wis
end

"""
    function observe
    Take a observation at location given by weights wis
"""
function observe(c::Field{T},wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}}},γ::Grid)::Vector{T} where T <: Real

        # look at total weight, < 1 if there are land points
    # later make sure total weight = 1 for proper average
    sumwis = Vector{Float64}(undef,length(wis))
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)
    [sumwis[i] = wetwrap[wis[i]...] for i in eachindex(wis)]

    # sample the true field at these random locations
    y = Vector{Float64}(undef,length(wis))
    replace!(c.tracer,NaN=>0.0)
    ywrap = view(c.tracer,list,:,:)
    [y[i] = ywrap[wis[i]...]/sumwis[i] for i in eachindex(wis)]

    return y
end

"""
    function gobserve(gy::Vector{T},c::Field{T},wis,γ) where T <: Real

    ADJOINT Take a observation at location given by weights wis
    Arguments not symmetric with `observe` due to splat operator
"""
function gobserve(gy::Vector{T},c::Field{T},locs) where T <: Real

    #initialize gc this sneaky way
    gc = 0.0 * c
    for ii in eachindex(gy)
        # interpweights repeats some calculations
        gc.tracer[c.γ.wet] .+= gy[ii] * interpweights(locs[ii],c.γ)[c.γ.wet]
    end

    return gc
end

function location_obs(field, locs, γ)

    tlength = length(size(field)) > 3 ? size(field)[1] : 1
    #determine weights for locations 
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    #
    sumwis = Vector{Float64}(undef,length(wis))
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)
    [sumwis[i] = wetwrap[wis[i]...] for i in eachindex(wis)]
    
    field_sample = Matrix{Float64}(undef,(length(wis),tlength))
    replace!(field,NaN=>0.0)
    if tlength > 1
        for t in 1:tlength
            field_wrap = view(field[t,:,:,:], list, :, :)
            [field_sample[i,t] = field_wrap[wis[i]...]/sumwis[i] for i in eachindex(wis)]
        end

        #below code is an attempt to replace for loop 
        
        #field_wrap = view(field, :, list, :, :) #handle wraparound
        #[field_sample[i,:] = field_wrap[:,wis[i][1],wis[i][2],wis[i][3]]/sumwis[i] for i in eachindex(wis)]

    else
    
    field_wrap = view(field,list,:,:)
    [field_sample[i] = field_wrap[wis[i]...]/sumwis[i] for i in eachindex(wis)]
    end
    
    return field_sample
end

""" 
    function costfunction_gridded_obs(uvec::Vector{T},Alu,b::BoundaryCondition{T},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real

    squared model-data misfit for gridded data
    controls are a vector input for Optim.jl
# Arguments
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
- `u`: controls, field format
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary conditions
- `y`: observations on grid
- `Wⁱ`: inverse of W weighting matrix for observations
- `γ`: grid
"""
function costfunction_gridded_obs(uvec,Alu,b₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}

    # turn uvec into a boundary condition
    u = unvec(u₀,uvec)

    b = adjustboundarycondition(b₀,u) #b += u # easy case where u and b are on the same boundary
    n = steadyinversion(Alu,b,γ) - y  # gives the misfit
    J = n ⋅ (Wⁱ * n) # dot product

    # adjoint equations
    gy = -2Wⁱ * n
    gb = gsteadyinversion( gy, Alu, b, γ)
    gu = gadjustboundarycondition(gb,u)
    guvec = vec(gu)

    return J, guvec
end

"""
    function costfunction_gridded_obs!(J,guvec,uvec::Vector{T},Alu,b₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}
"""
function costfunction_gridded_obs!(J,guvec,uvec::Vector{T},Alu,b₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}

    # turn uvec into a boundary condition
    u = unvec(u₀,uvec)

    b = adjustboundarycondition(b₀,u) #b += u # easy c
    #adjustboundarycondition!(b,u) # easy case where u and b are on the same boundary
    y -= steadyinversion(Alu,b,γ)  # gives the misfit

    if guvec != nothing
        # adjoint equations
        gy = -2Wⁱ * y
        gb = gsteadyinversion( gy, Alu, b, γ)
        gu = gadjustboundarycondition(gb,u)
        #guvec = vec(gu)

        # next block just to modify the contents
        tmp = vec(gu)
        for (ii,vv) in enumerate(tmp)
            guvec[ii] = vv
        end
    end
    
    if J !=nothing
        return  y ⋅ (Wⁱ * y) # dot product
    end
end

""" 
    function costfunction_point_obs(uvec::Vector{T},Alu,b₀::BoundaryCondition{T},u₀::BoundaryCondition{T},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where T <: Real

    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue #1: couldn't figure out how to nest with costfunction_obs!
    
# Arguments
- `uvec`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary condition
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `Q⁻`: weights for control vector
- `γ`: grid
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
function costfunction_point_obs(uvec::Vector{T},Alu,b₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where {N1, N2, T <: Real}

    # control penalty and gradient
    Jcontrol = uvec'*(Q⁻*uvec)
    guvec = 2*(Q⁻*uvec)

    u = unvec(u₀,uvec)
    b = adjustboundarycondition(b₀,u) # combine b₀, u

    c = steadyinversion(Alu,b,γ)  # gives the misfit

    # observe at right spots
    ỹ = observe(c,wis,γ)
    n = ỹ - y
    
    Jdata = n ⋅ (Wⁱ * n) # dot product
    J = Jdata + Jcontrol

    gn = 2Wⁱ * n
    gỹ = gn
    
    gc = gobserve(gỹ,c,locs)
    gb = gsteadyinversion(gc, Alu, b, γ)
    gu = gadjustboundarycondition(gb,u)
    guvec += vec(gu)
    
    return J, guvec
end

""" 
    function costfunction_point_obs!(J,guvec,uvec::Vector{T},Alu,b₀::BoundaryCondition{T},u₀::BoundaryCondition{T},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where T <: Real

    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue #1: couldn't figure out how to nest with costfunction_obs!
    
# Arguments
- `J`: cost function of sum of squared misfits
- `guvec`: derivative of cost function wrt to controls
- `uvec`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary condition
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `Q⁻`: weights for control vector
- `γ`: grid
"""
function costfunction_point_obs!(J,guvec,uvec::Vector{T},Alu,b₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where {N1, N2, T <: Real}

    u = unvec(u₀,uvec)
    b = adjustboundarycondition(b₀,u) # combine b₀, u
    c = steadyinversion(Alu,b,γ)  # gives the misfit

    # observe at right spots
    ỹ = observe(c,wis,γ)
    n = ỹ - y

    if guvec != nothing
        # adjoint equations
        gtmp = 2*(Q⁻*uvec)
        gn = 2Wⁱ * n
        gỹ = gn
        gc = gobserve(gỹ,c,locs)
        gb = gsteadyinversion(gc, Alu, b, γ)
        gu = gadjustboundarycondition(gb,u)
        gtmp += vec(gu)
        for (ii,vv) in enumerate(gtmp)
            guvec[ii] = vv
        end
    end

    if J !=nothing
        # control penalty and gradient
        Jcontrol = uvec'*(Q⁻*uvec)
        Jdata = n ⋅ (Wⁱ * n) # dot product
        return Jdata + Jcontrol
    end
end

""" 
    function steadyinversion(Alu,b;q=nothing,r=1.0)
    invert for a steady-state tracer distribution
# Arguments
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary condition
- `γ`::Grid
# Optional Arguments
- `q`: interior sources/sinks of phosphate
- `r`: stochiometric ratio of tracer:phosphate
# Output
- `c`::Field, steady-state tracer distribution
"""
function steadyinversion(Alu,b::BoundaryCondition{T},γ::Grid;q=nothing,r=1.0)::Field{T} where T <: Real

    #println("running steady inversion")

    # preallocate Field for equation constraints
    d = zeros(γ)
    
    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end

    # define ldiv with fields
    c = zeros(d.γ)
    #c.tracer[c.γ.wet] =  Alu\(d.tracer[d.γ.wet])
    c = Alu \ d
    
    return c
end

""" 
    function gsteadyinversion(Alu,b;q=nothing,r=1.0)

    ADJOINT invert for a steady-state tracer distribution

# Arguments
- `Alu`: LU decomposition of water-mass matrix
- `b`: BoundaryCondition
- `γ`::Grid
# Optional Arguments
- `q`: interior sources/sinks of phosphate
- `r`: stochiometric ratio of tracer:phosphate
# Output
- `c`::Field, steady-state tracer distribution
"""
function gsteadyinversion(gc,Alu,b::BoundaryCondition{T},γ::Grid;q=nothing,r=1.0)::BoundaryCondition{T} where T <: Real

    #println("running adjoint steady inversion")

    gd = Alu' \ gc

    # still need to develop this
    # and to find a way to return the value
    if !isnothing(q)
        gq = gsetsource!(gd,d,q,r)
    end

    gb = gsetboundarycondition(gd,b)
    return gb
end

"""
    function steadyinversion(Alu,b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}},γ::Grid;q=nothing,r=1.0)::Field{T} where {N, T <: Real}

    steady inversion for b::NamedTuple
"""
function steadyinversion(Alu,b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}},γ::Grid;q=nothing,r=1.0)::Field{T} where {N, T <: Real}

    # preallocate Field for equation constraints
    d = zeros(γ)

    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end

    # define ldiv with fields
    c = zeros(d.γ)
    c.tracer[c.γ.wet] =  Alu\(d.tracer[d.γ.wet])
    #c = Alu \ d

    return c
end

"""
    function gsteadyinversion(gc::Field{T},Alu,b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}},γ::Grid;q=nothing,r=1.0)::Field{T} where {N, T <: Real}

    ADDJOINT steady inversion for b::NamedTuple
"""
function gsteadyinversion(gc::Field{T},Alu,b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}},γ::Grid;q=nothing,r=1.0) where {N, T <: Real}

    #println("running adjoint steady inversion")

    # sensitivity of Field of equation constraints
    gd = Alu' \ gc

    # still need to develop this
    # and to find a way to return the value
    if !isnothing(q)
        gq = gsetsource!(gd,d,q,r)
    end

    # update d with the boundary condition b
    gb = gsetboundarycondition(gd,b)

    return gb
end

"""
    function wetlocation(γ)
    Get (lon,lat,depth) tuples of wet locations.
    Allow a location to be wet if at least one out of 8 nearby gridpoints is wet.
    Certainly "wet" gridpoints could be defined more strictly.
# Arguments
- `γ`: TMI.grid
# Output
- `loc`: lon,lat,depth """
function wetlocation(γ)

    confirmwet = false
    neighbors  = 8
    while !confirmwet
        loc = (rand(minimum(γ.lon):0.1:maximum(γ.lon)),
               rand(minimum(γ.lat):0.1:maximum(γ.lat)),
               rand(minimum(γ.depth):1.0:maximum(γ.depth)))

        iswet(loc,γ) && return loc
        println("dry point, try again")
    end # if not, then start over.
end
    
function iswet(loc,γ,neighbors)
    # two approaches
    # approach 2
    # find 8 nearest neighbors
    Inn = nearestneighbor(loc,γ,neighbors)

    # are any of them wet?
    for ii = 1:neighbors
        if γ.wet[Inn[ii]]
            return true
        end
    end
    return false
end

function iswet(loc,γ)

    # wetness bounded by 0 and 1
    # should be an argument
    # 1 = very strict
    # 0 = all points
    wetness = 0.2
    
    wis = interpindex(loc,γ)

    # handle wraparound
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)

    # are any of them wet?
    # interpolate ones and zeros on to this loc.
    # if there is land nearby, the interpolated value
    # will be greater than 0.
    # this criterion only requires on land point nearby,
    # where nearby is one of the 8 corners of the cube that contains loc
    return wetwrap[wis...] > wetness
end

# define the correct dimension and index for each control plane
# maybe someday find a way to hide γ
zerosurfaceboundary(γ) = zeros(3,1,γ)::BoundaryCondition
zeronorthboundary(γ) = zeros(2,maximum(latindex(γ.I)),γ)::BoundaryCondition
zeroeastboundary(γ) = zeros(1,maximum(lonindex(γ.I)),γ)::BoundaryCondition
zerosouthboundary(γ) = zeros(2,1,γ)::BoundaryCondition
zerowestboundary(γ) = zeros(1,1,γ)::BoundaryCondition

onesurfaceboundary(γ) = ones(3,1,γ)::BoundaryCondition

getsurfaceboundary(c::Field) = getboundarycondition(c.tracer,3,1,c.γ)::BoundaryCondition
getnorthboundary(c::Field) = getboundarycondition(c.tracer,2,maximum(latindex(c.γ.I)),c.γ)::BoundaryCondition
geteastboundary(c::Field) = getboundarycondition(c.tracer,1,maximum(lonindex(c.γ.I)),c.γ)::BoundaryCondition
getsouthboundary(c::Field) = getboundarycondition(c.tracer,2,1,c.γ)::BoundaryCondition
getwestboundary(c::Field) = getboundarycondition(c.tracer,1,1,c.γ)::BoundaryCondition

""" 
    function setboundarycondition!(d::Field,b::BoundaryCondition)
    apply boundary condition to the equation constraints
# Arguments
- `d`::Field, equation constraints (i.e., right hand side)
- `b`::BoundaryCondition
"""
function setboundarycondition!(d::Field{T},b::BoundaryCondition{T}) where T<: Real

    if b.dim == 1
        d.tracer[b.dimval,:,:] += b.tracer
    elseif b.dim == 2
        d.tracer[:,b.dimval,:] += b.tracer
    elseif b.dim == 3
        d.tracer[:,:,b.dimval] += b.tracer
    else
        error("controls not implemented for 4+ dimensions")
    end
    return d
end

""" 
    function gsetboundarycondition(gd::Field{T},b::BoundaryCondition{T}) where T<: Real

    ADJOINT: apply boundary condition to the equation constraints
# Arguments
- `d`::Field, equation constraints (i.e., right hand side)
- `b`::BoundaryCondition
"""
function gsetboundarycondition(gd::Field{T},b::BoundaryCondition{T}) where T<: Real
    #gb = 0.0 * b # initialize to zero
    if b.dim == 1
        #gb.tracer = gd.tracer[b.dimval,:,:]
        gb = BoundaryCondition(gd.tracer[b.dimval,:,:],b.i,b.j,b.k,b.dim,b.dimval,b.wet)
    elseif b.dim == 2
        #gb.tracer = gd.tracer[:,b.dimval,:]
        gb = BoundaryCondition(gd.tracer[:,b.dimval,:],b.i,b.j,b.k,b.dim,b.dimval,b.wet)
    elseif b.dim == 3
        #gb.tracer .+= gd.tracer[:,:,b.dimval] 
        gb = BoundaryCondition(gd.tracer[:,:,b.dimval],b.i,b.j,b.k,b.dim,b.dimval,b.wet)
    else
        error("controls not implemented for 4+ dimensions")
    end
    return gb
end

"""
    function setboundarycondition!(d::Field{T},b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}}) where {N, T <: Real}

    set all boundary conditions in a Named Tuple
"""
function setboundarycondition!(d::Field{T},b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}}) where {N, T <: Real}
    for b1 in b
        setboundarycondition!(d,b1)
    end
end

""" 
    function gsetboundarycondition(gd::Field{T},b::BoundaryCondition{T}) where T<: Real

    ADJOINT: apply boundary condition to the equation constraints
# Arguments
- `d`::Field, equation constraints (i.e., right hand side)
- `b`::BoundaryCondition
"""
function gsetboundarycondition(gd::Field{T},b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}}) where {N, T<: Real}

    gb1 = Vector{BoundaryCondition{T}}(undef,length(keys(b)))
    for (ii,vv) in enumerate(b)
        gb1[ii] = gsetboundarycondition(gd,vv)
    end

    # https://discourse.julialang.org/t/construct-namedtuple-dynamically/15394/7
    gb = (;zip(keys(b), gb1)...)
                                       
    return gb
end

"""
    function adjustboundarycondition!(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real

    adjust the (one) boundary condition 
    problem: passes back a mutated b
"""
function adjustboundarycondition!(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real
    b.tracer[b.wet] += u.tracer[u.wet] # write it out so b changes when returned
end

"""
    function adjustboundarycondition(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real

    adjust the (one) boundary condition 
"""
function adjustboundarycondition(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real

    bnew = b + u
    return bnew

end

"""
    function gadjustboundarycondition!(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real

    adjust the (one) boundary condition 
    Just copy the variable.
    Keep this function so that calling functions can look alike.
    Could probably combine with lower function, use Union type
"""
function gadjustboundarycondition(gb::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real
    gu  = gb
    return gu
end

"""
    function adjustboundarycondition!(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where N1, N2, T <: Real

    adjust all boundary conditions b that are described in u
"""
function adjustboundarycondition!(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where {N1, N2, T <: Real}

    ukeys = keys(u)
    for ukey in keys(u)
        b[ukey].tracer[b[ukey].wet] += u[ukey].tracer[b[ukey].wet] 
    end
    
end

"""
    function adjustboundarycondition(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where N1, N2, T <: Real

    adjust all boundary conditions b that are described in u
"""
function adjustboundarycondition(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where {N1, N2, T <: Real}

    bnew = deepcopy(b)
    ukeys = keys(u)
    for ukey in keys(u)
        bnew[ukey].tracer[bnew[ukey].wet] += u[ukey].tracer[bnew[ukey].wet] 
    end
    return bnew
end

"""
    function gadjustboundarycondition!(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real

    ADJOINT CODE
    adjust the (one) boundary condition 
    Just copy the variable.
    Keep this function so that calling functions can look alike.
"""
function gadjustboundarycondition(gb::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where {N1, N2, T <: Real}
    gu = gb[keys(u)] # grab the parts of the named tuple corresponding to u
    return gu
end

"""
    function section
    View latitude-depth slice of field
# Arguments
- `c::Field`, 3D tracer field plus meta data
- `lon`: longitude of section
# Output
- `csection`: 2d slice of field
"""
function section(c::Field{T},lon)::Array{T,2} where T <: Real

    # handle longitudinal ambiguities
    if lon < minimum(c.γ.lon)
        lon += 360
    end
    
    isec = findall(==(lon),c.γ.lon)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    csection= dropdims(view(c.tracer,isec,:,:),dims=1)
    return csection
end

function planview(c::Field{T},depth)::Array{T,2} where T <: Real
 
    isec = findall(==(depth),c.γ.depth)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    cplan = dropdims(view(c.tracer,:,:,isec),dims=3)
    return cplan
end

"""
    function plotextent
    Generate image showing user-specified ROI
# Arguments
- `latbox`: in format [lat_start, lat_stop]
- `lonbox`: in format [lon_start, lon_stop]

"""
function plotextent(latbox, lonbox)
    
#    ccrs = pyimport("cartopy.crs")
    lower_left = [minimum(lonbox), minimum(latbox)] #array of lower left of box

    #calc width and height of box
    w = maximum(lonbox) - minimum(lonbox)
    h = maximum(latbox) - minimum(latbox)

    #init GeoAxes
    fig = figure()
    ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())

    #plot rectangle
    ax.add_patch(TMI.mpl.patches.Rectangle(xy=lower_left,
                                 width=w, height=h,
                                 facecolor="blue",
                                 alpha=0.2,
                                 transform=TMI.cartopy.crs.PlateCarree()))
    #define extent of figure
    pad = 10 #how many deg lat and lon to show outside of bbox
    pad_add = [-pad, pad] #add this to latbox and lonbox
    padded_lat = latbox + pad_add
    padded_lon = lonbox + pad_add
    ext = vcat(padded_lon, padded_lat) #make into one vector
    ax.set_extent(ext)

    # using cartopy 0.18 and NaturalEarth is missing
    ax.coastlines() #show coastlines

    #add gridlines
    gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
    gl.top_labels = false
    gl.right_labels = false

    ax.set_title("User-defined surface patch")
end

"""
    function sectionplot
    Plot of section (lat-depth) in ocean
# Arguments
- `field::Field`, 3d filed of values to be plotted
- `lon`: longitude of section
- `lims`: contour levels
- `titlelabel`: optional title labeln
"""
function sectionplot(field::Field{T}, lon, lims;titlelabel="section plot") where T <: Real

    Psection = section(field,lon)
    cmap_seismic = get_cmap("seismic")
    z = field.γ.depth/1000.0
    
    #calc fignum - based on current number of figures
    figure()
    contourf(field.γ.lat, z, Psection', lims, cmap=cmap_seismic)
    #fig, ax = plt.subplots()
    CS = gca().contour(field.γ.lat, z, Psection', lims,colors="k")
    gca().clabel(CS, CS.levels, inline=true, fontsize=10)
    xlabel("Latitude [°N]")
    ylabel("Depth [km]")
    gca().set_title(titlelabel)
    gca().invert_yaxis()
    colorbar(orientation="horizontal")
    
end

"""
    function planviewplot
    Plot of plan view (lon-lat) in ocean
# Arguments
- `field::Field`, 3d filed of values to be plotted
- `depth`: depth of plan view
- `lims`: contour levels
- `titlelabel`: optional title label
"""
function planviewplot(c::Field{T}, depth, lims;titlelabel="section plot") where T <: Real

    cplan = planview(c::Field{T},depth)

    cmap_seismic = get_cmap("seismic")
    
    #calc fignum - based on current number of figures
    figure()
    contourf(c.γ.lon,c.γ.lat, cplan', lims, cmap=cmap_seismic)
    #fig, ax = plt.subplots()
    CS = gca().contour(c.γ.lon,c.γ.lat, cplan', lims, cmap=cmap_seismic)
    gca().clabel(CS, CS.levels, inline=true, fontsize=10)
    ylabel("Latitude [°N]")
    xlabel("Longitude [°E]")
    gca().set_title(titlelabel)
    colorbar(orientation="vertical")
    
end

"""
    function planviewplot
    Plot of plan view (lon-lat) in ocean
# Arguments
- `field::BoundaryCondition`, 3d filed of values to be plotted
- `depth`: depth of plan view
- `lims`: contour levels
- `γ::Grid`, needed for lat, lon but not in BoundaryCondition! (could refactor)
- `titlelabel`: optional title label
"""
function planviewplot(b::BoundaryCondition{T}, lims,γ::Grid;titlelabel="surface plot") where T <: Real

    # is the boundary condition oriented correctly?
    if b.dim != 3
        error("boundary condition not horizontal")
    end
    
    cplan = b.tracer
    
    cmap_seismic = get_cmap("seismic")
    
    #calc fignum - based on current number of figures
    figure()
    contourf(γ.lon,γ.lat, cplan', lims, cmap=cmap_seismic)
    #fig, ax = plt.subplots()
    CS = gca().contour(γ.lon,γ.lat, cplan', lims, cmap=cmap_seismic)
    gca().clabel(CS, CS.levels, inline=true, fontsize=10)
    ylabel("Latitude [°N]")
    xlabel("Longitude [°E]")
    gca().set_title(titlelabel)
    colorbar(orientation="vertical")
    
end

""" 
n    function control2state(tracer2D,γ)
    turn 2D surface field into 3D field with zeroes below surface    
# Arguments
- `tracer2D`:: 2D surface tracer field
- `wet`::BitArray mask of ocean points
# Output
- `tracer3D`:: 3d tracer field with NaN on dry points
"""
function control2state(tracer2D::Matrix{T},wet) where T<: Real
    # preallocate
    tracer3D = Array{T}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer3D[wet] .= zero(T)
    tracer3D[.!wet] .= zero(T)/zero(T)
    tracer3D[:,:,1] = tracer2D
    return tracer3D
end

""" 
    function surfacecontrol2field(usfc,γ.wet)
    turn surface control vector into 3D field with zeroes below surface    
# Arguments
- `usfc`:: surface control vector
- `wet`::BitArray mask of ocean points
# Output
- `tracer3D`:: 3d tracer field with NaN on dry points
"""
function surfacecontrol2field(usfc::Vector{T},wet) where T<: Real
    # preallocate
    tracer3D = Array{T}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer3D[wet] .= zero(T)
    tracer3D[.!wet] .= zero(T)/zero(T)
    tracer3D[:,:,1][wet[:,:,1]] = usfc
    return tracer3D
end

"""
    function Γsfc 
    Γsfc anonymously calls surfacecontrol2field
"""
Γsfc = surfacecontrol2field

"""
    function surfacecontrol2field!(c,u,γ)
    Add surface control vector to existing 3D field 
# Arguments
- `c`:: state field, 3d tracer field with NaN on dry points, modified by function
- `usfc`:: surface control vector
- `wet`::BitArray mask of ocean points
"""
function surfacecontrol2field!(c::Array{T,3},usfc::Vector{T},γ) where T<: Real
    #c[:,:,1][wet[:,:,1]] .+= u # doesn't work
#    [c[γ.I[ii][1],γ.I[ii][2],γ.I[ii][3]] += u[ii] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
    list = surfaceindex(γ.I)
    [c[γ.I[ii]] += usfc[list[ii]] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
end

""" 
    function surfacecontrol2field!(c,u,γ)
    Add surface control vector to tracer vector
# Arguments
- `c`:: state field, 3d tracer field with NaN on dry points, modified by function
- `u`:: surface control vector
- `wet`::BitArray mask of ocean points
"""
function surfacecontrol2field!(c::Vector{T},u::Vector{T},γ) where T<: Real
    list = surfaceindex(γ.I)
    [c[ii] += u[list[ii]] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
end

"""
    function Γsfc! 
    Γsfc! anonymously calls surfacecontrol2field!
"""
Γsfc! = surfacecontrol2field!

function field2obs(cvec,wis,γ)
    # interpolate onto data points
    N = length(wis)
    sumwis = Vector{Float64}(undef,N)
    list = vcat(1:length(γ.lon),1)

    # perhaps the most clever line in TMI.jl?
    wetwrap = view(γ.wet,list,:,:)

    # some interpolation weights on land, oh no
    # sum up all weights in ocean
    [sumwis[i] = wetwrap[wis[i]...] for i in eachindex(wis)]

    # reconstruct the observations
    ỹ = Vector{Float64}(undef,N)
    c̃ = tracerinit(γ.wet)
    c̃[γ.wet] = cvec
    replace!(c̃,NaN=>0.0)
    cwrap = view(c̃,list,:,:)

    # divide by sum of all ocean weights so that this is still a true average
    [ỹ[i] = cwrap[wis[i]...]/sumwis[i] for i in eachindex(wis)]
    return ỹ
end

"""
    function E 
    E anonymously calls field2obs
"""
E = field2obs

""" 
    function tracerinit(wet,vec,I)
          initialize tracer field on TMI grid
        perhaps better to have a tracer struct and constructor
# Arguments
- `wet`:: BitArray mask of ocean points
- `vec`:: vector of values at wet points
- `I`:: Cartesian Index for vector
# Output
- `field`:: 3d tracer field with NaN on dry points
"""
function tracerinit(vec,I,wet)

    # preallocate
    T = eltype(vec)
    field = Array{T}(undef,size(wet))
    fill!(field,zero(T)/zero(T))    

    #- a comprehension
    [field[I[n]]=vec[n] for n ∈ eachindex(I)]
    return field
end

"""
    function vec2fld
    Transfer a vector to a 3D field with accounting for ocean bathymetry
# Arguments
- `vector`: field in vector form (no land points)
- `I`: cartesian indices of ocean points
# Output
- `field`: field in 3d form including land points (NaN)
"""
function vec2fld(vector::Vector{T},I::Vector{CartesianIndex{3}}) where T<:Real
    # choose NaN for now, zero better? nothing better?
    fillvalue = zero(T)/zero(T)
    
    nx = maximum(I)[1]
    ny = maximum(I)[2]
    nz = maximum(I)[3]

    # faster instead to allocate as undef and then fill! ?
    field = (NaN .* zero(T)) .* zeros(nx,ny,nz)

    # a comprehension
    [field[I[n]]=vector[n] for n ∈ eachindex(I)]
    return field
end

"""
    function fld2vec
    Transfer 3D field with accounting for ocean bathymetry to a vector without land points.
    This is done more easily with a BitArray mask, i.e., vector = field[mask].
    This function may be removed in the future.
# Arguments
- `field`: field in 3d form including land points (NaN)
- `I`: cartesian indices of ocean points
# Output
- `vector`: field in vector form (no land points)
"""
function fld2vec(field::Array{T,3},I::Vector{CartesianIndex{3}}) where T<:Real
    vector = Vector{T}(undef,length(I))
    #- a comprehension
     [vector[n] = field[I[n]] for n ∈ eachindex(I)];
     return vector
 end

end
