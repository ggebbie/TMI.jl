module TMI

using Revise
using LinearAlgebra, SparseArrays, NetCDF, Downloads,
    GoogleDrive, Distances, DrWatson, GibbsSeaWater,  
    PyPlot, PyCall, Distributions, Optim

export config, download,
    vec2fld, fld2vec, depthindex, surfaceindex,
    surfacepatch, section,
    layerthickness, cellarea, cellvolume,
    planview, dyeplot, plotextent, tracerinit,
    updateLinearindex,
    watermassmatrixXYZ, watermassmatrixZYX,
    linearindexXYZ, nearestneighbor,
    nearestneighbormask, horizontaldistance,
    readtracer, cartesianindexZYX, Γ,
    misfit_gridded_data, misfit_gridded_data!,
    trackpathways, regeneratedphosphate, volumefilled,
    surfaceorigin, sample_observations, filterdata

#export JULIA_SSL_NO_VERIFY_HOSTS:"naturalearth.s3.amazonaws.com"

#Python packages - initialize them to null globally
#const patch = PyNULL()
#const ccrs = PyNULL()

# from ClimatePlots.jl
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

struct grid
    lon::Vector{Real}
    lat::Vector{Real}
    depth::Vector{Real}
    I::Vector{CartesianIndex{3}} # index
    R::Array{Integer,3}
#        R::Array{Union{Integer,Nothing},3}
#    R::LinearIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}} 
    wet::BitArray{3}
end

"""
    function config(url,inputdir)
# Arguments
- `TMIversion`: TMI version for water-mass/circulation model
# Output
- `A`: TMI steady-state water-mass matrix
- `Alu`: LU decomposition of A
- `γ`: TMI grid properties
- `TMIfile`: TMI file name
"""
function config(TMIversion)

    #- `url`: Google Drive URL for data
    url = gdriveurl(TMIversion)
    
    TMIfile = datadir(TMIversion*".nc")
    !isdir(datadir()) ? mkpath(datadir()) : nothing
    !isfile(TMIfile) ? download(url,datadir()) : nothing

    ncdata = NetCDF.open(TMIfile)
    #println(ncdata)
    
    # move this to runtests.jl to see if it is read correctly
    # Azyx = watermassmatrixZYX(TMIfile)

    # make a sample field from zyx cartesian indices
    Izyx = cartesianindexZYX(TMIfile)

    # make a mask
    wet = BitArray{3}(undef,maximum(Izyx)[1],maximum(Izyx)[2],maximum(Izyx)[3])
    fill!(wet,0)
    wet[Izyx] .= 1

    # if a tracer is available, should be consistent with this definition
    #wet = .!isnan.(c)
    
    # need to write this function
    I = cartesianindexXYZ(wet)

    R = linearindexXYZ(wet)

    A = watermassmatrixXYZ(TMIfile,R)
    #R = R[ wet ] # eliminate land points

    # LU factorization for efficient matrix inversions
    Alu = lu(A)
    
    # get properties of grid
    lat,lon,depth = gridprops(TMIfile)

    γ = grid(lon,lat,depth,I,R,wet)

    return  A, Alu, γ, TMIfile

end

"""
    function cartesianindexZYX(file)
    Read and assemble the grid coordinates
    according to the legacy MATLAB code (z,y,x order).
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function cartesianindexZYX(file)
    # make the Cartesian tracer grid
    it = convert(Array{Int,1},ncread(file,"xgrid"))
    jt = convert(Array{Int,1},ncread(file,"ygrid"))
    kt = convert(Array{Int,1},ncread(file,"zgrid"))
    I = CartesianIndex.(it,jt,kt)
    return I
end

"""
    function cartesianindexXYZ(wet)
    Read and assemble the grid coordinates
    according to a 3D tracer in x,y,z order
# Arguments
- `wet`: BitArray logical mask for wet points
# Output
- `I`: 3D Cartesian indices
"""
cartesianindexXYZ(wet) = findall(wet)

"""
    function linearindexXYZ(file)
    Read and assemble the grid coordinates.
# Arguments
- `wet`: 3D mask for wet points
# Output
- `R`: array of linear indices, but not a LinearIndices type
"""
function linearindexXYZ(wet)
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
    function gridprops(file)
    Read and assemble the grid properties.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function gridprops(file)
    lat = ncread(file,"lat")
    lon = ncread(file,"lon")
    depth = ncread(file,"depth")
    return lat,lon,depth
end

"""
    function watermassmatrixZYX(file)
    Read and assemble the water-mass matrix.
    Legacy version from MATLAB.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `A`: water-mass matrix
"""
function watermassmatrixZYX(file)
    i = ncread(file,"i")
    j = ncread(file,"j")
    m = ncread(file,"m")
    A = sparse(i,j,m)
    return A
end

"""
        function watermassmatrixXYZ(file,R)
    Read and assemble the water-mass matrix from MATLAB.
    Transfer to updated x,y,z version
# Arguments
- `file`: TMI NetCDF file name
- `γ`: TMI grid
# Output
- `A`: water-mass matrix
"""
function watermassmatrixXYZ(file,R)

    # MATLAB accounting z,y,x
    izyx = convert(Vector{Int},ncread(file,"i"))
    jzyx = convert(Vector{Int},ncread(file,"j"))

    # Julia accounting x,y,z
    Izyx = cartesianindexZYX(file)
    ixyz = updatelinearindex(izyx,Izyx,R)
    jxyz = updatelinearindex(jzyx,Izyx,R)
    
    m = ncread(file,"m")

    # use grid indices to switch i,j values
    A = sparse(ixyz,jxyz,m)
    #A = A[wet,:];
    return A
end


"""
    function updatelinearindex(izyx,Izyx,R)
    Linear index translated from z,y,x to x,y,z accounting
# Arguments
- `izyx`: index of interest in z,y,x accounting
- `Izyx`: wet Cartesian Index for z,y,x
- `R`: Linear indices for x,y,z 
# Output
- `ixyz`: index of interest in x,y,z accounting
"""
function updatelinearindex(izyx,Izyx,R)
    # get Izyx Cartesian index stored from legacy MATLAB code
    ixyz = R[Izyx[izyx]]
    return ixyz
end

"""
    function watermassmatrix(file)
    Read and assemble the water-mass matrix.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `A`: water-mass matrix
"""
function watermassmatrix(file)

    m = massFractions(file)

    A = watermassmatrix(m)
    
    # assemble m into A
    # Atest = Array{SparseArrays.SparseMatrixCSC{Float64, Int64},3}

    
    # it = convert(Vector{Int},ncread(file,"xgrid"))
    # jt = convert(Vector{Int},ncread(file,"ygrid"))
    # kt = convert(Vector{Int},ncread(file,"zgrid"))
    # i = convert(Vector{Int},ncread(file,"i"))
    # j = convert(Vector{Int},ncread(file,"j"))
    # m = ncread(file,"m")

    # for gg = 1:length(i)

    #     list = findall(x -> x == gg, i)

    #     # get linear coordinates.
    #     [linear[q] = R[it[list[q]],jt[list[q]],kt[list[q]]] for q = 1:length(list)] 
        
    #     # can't store sparse 3D matrix. just 2D.
    #     # need to go from i,j,k to linear index
    #     Arowsparse = sparse(  m[gg])
    #     A[it[i[gg]],jt[i[gg]],kt[i[gg]]] = sparse([it[j[gg]],jt[j[gg]],kt[j[gg]]] = m[gg]
    # end
    
    #A = sparse(i,j,m)
    return A
end

"""
    function watermassmatrix(file)
    Assemble the water-mass matrix given mass fractions `m`
# Arguments
- `m`: mass fractions
# Output
- `A`: water-mass matrix
"""
# function watermassmatrix(m)
#     # Goal:assemble m into A

#     # preallocate A
#     A = Array{SparseArrays.SparseMatrixCSC{Float64, Int64},3}
                                                  
#     # loop over each equation
    
#     # get linear index for equation (destination)

#     # get linear indices for sources

#     # make a list of row (destination), column (destination), m value

#     # complete loop

#     # make sparse matrix
#     #A = sparse(i,j,m)
#     return A
# end

"""
    function massFractions(file)
    Read and assemble the water-mass fractions, `m`
# Arguments
- `file`: TMI NetCDF file name
# Output
- `m`: mass fractions
"""
# function massFractions(file)
#     Atest = Array{SparseArrays.SparseMatrixCSC{Float64, Int64},3}
    
#     it = convert(Vector{Int},ncread(file,"xgrid"))
#     jt = convert(Vector{Int},ncread(file,"ygrid"))
#     kt = convert(Vector{Int},ncread(file,"zgrid"))
#     i = convert(Vector{Int},ncread(file,"i"))
#     j = convert(Vector{Int},ncread(file,"j"))
#     mMat = ncread(file,"m") # MATLAB generated list

#     for gg = 1:length(i)

#         list = findall(x -> x == gg, i)

#         # get linear coordinates.
#         [linear[q] = R[it[list[q]],jt[list[q]],kt[list[q]]] for q = 1:length(list)] 
        
#         # can't store sparse 3D matrix. just 2D.
#         # need to go from i,j,k to linear index
#         Arowsparse = sparse(  m[gg])
#         A[it[i[gg]],jt[i[gg]],kt[i[gg]]] = sparse([it[j[gg]],jt[j[gg]],kt[j[gg]]] = m[gg]
#     end
    
#     return m
# end
                                                  
"""
    function readtracer(file,tracername)
    Read and assemble the water-mass matrix.
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

function layerthickness(γ::grid)
    zface= (γ.depth[1:end-1].+γ.depth[2:end])./2;
    dz = ([zface[1] ; diff(zface); 500]);
    return dz
end

function zonalgriddist(γ::grid)
    dx = similar(γ.lat)
    for j in eachindex(γ.lat)
        dx[j] = haversine((γ.lon[1],γ.lat[j])
                         ,(γ.lon[2],γ.lat[j]))
    end
    return dx
end

"""
    function download(url,inputdir)
    Read and assemble all TMI inputs.
# Arguments
- `url`: Google Drive location of TMI input
- `inputdir`: input directory location to store file
# Output
- none
"""
function download(url,inputdir)
    # later include options and move these settings to arguments.

    # make sure input dir exists
    !isdir(inputdir) ? mkdir(inputdir) : nothing

    # two ways to download
    # 1. Use `run` for a shell command (less portable). Also difficult for Google Drive. See downloadTMIfromGoogleDrive.sh.

    # 2. Use GoogleDrive.jl package
    google_download(url,inputdir)
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
function vec2fld(vector::Vector{Float64},I::Vector{CartesianIndex{3}})

    nx = maximum(I)[1]
    ny = maximum(I)[2]
    nz = maximum(I)[3]
    field = NaN .* zeros(nx,ny,nz)

    # a comprehension
    [field[I[n]]=vector[n] for n ∈ eachindex(I)]
    return field
end

"""
    function fld2vec
    Transfer 3D field with accounting for ocean bathymetry to a vector without land points
# Arguments
- `field`: field in 3d form including land points (NaN)
- `I`: cartesian indices of ocean points
# Output
- `vector`: field in vector form (no land points)
"""
function fld2vec(field::Array{Float64,3},I::Vector{CartesianIndex{3}})
    vector = Vector{Real}(undef,length(I))
    #- a comprehension
     [vector[n] = field[I[n]] for n ∈ eachindex(I)];
     return vector
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
function surfacepatch(lonbox::Vector{T},latbox::Vector{T},γ::grid)::Array{Float64} where T<:Real

    # ternary operator to handle longitudinal wraparound
    lonbox[1] ≤ 0 ? lonbox[1] += 360 : nothing
    lonbox[2] ≤ 0 ? lonbox[2] += 360 : nothing

    # define the surface boundary condition

    # preallocate
    d = tracerinit(γ.wet)

    [d[i,j,1] =  latbox[1] ≤ γ.lat[j] ≤ latbox[2] && lonbox[1] ≤ γ.lon[i] ≤ lonbox[2] for i in eachindex(γ.lon) for j in eachindex(γ.lat)] 
    d[.!γ.wet] .= NaN # double check that NaNs stay NaNs

    # old method for vectors
        #nfield = length(γ.I) # number of ocean points
    #d = zeros(Int,nfield) # preallocate
    #[d[n]=1 for n ∈ 1:nfield if γ.I[n][3]==1 && latbox[1] ≤ γ.lat[γ.I[n][2]] ≤ latbox[2]
    #     && lonbox[1] ≤ γ.lon[γ.I[n][1]] ≤ lonbox[2] ]
    return d
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
function nearestneighbormask(loc,γ::grid)

    Inn, Rnn = nearestneighbor(loc,γ)

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
- `Rnn`: linear indices of nearest neighbor
"""
function nearestneighbor(loc,γ)

    xydist = horizontaldistance(loc[1:2],γ)
    ijdist,ijmin = findmin(xydist[γ.wet[:,:,1]])
    
    kdist,kmin = findmin(abs.(loc[3] .- γ.depth))

    # translate ijmin into imin, jmin
    Inn = CartesianIndex.(γ.I[ijmin][1],γ.I[ijmin][2],kmin)
    Rnn = γ.R[Inn]

    # what if there is no Rnn? Outside of grid. Could move vertically by one.
    if iszero(Rnn)
        Inn = CartesianIndex.(γ.I[ijmin][1],γ.I[ijmin][2],kmin-1)
        Rnn = γ.R[Inn]
    end

    # if still not available, then give up.
    iszero(Rnn) && println("Warning: outside of grid")
    
    return Inn, Rnn
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
function horizontaldistance(loc,γ::grid)

    # pre-allocate horizontal distance
    hordist = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(hordist,NaN)
    
    # calculate haversine horizontal distance on sphere
    [hordist[γ.I[ii]] = haversine((loc[1],loc[2]),                  (γ.lon[γ.I[ii][1]],γ.lat[γ.I[ii][2]]))
       for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
    return hordist
end

"""
    function section
    View latitude-depth slice of field
# Arguments
- `c`: 3d tracer field
- `lon`: longitude of section
- `γ`: TMI.grid
# Output
- `csection`: 2d slice of field
"""
function section(c,lon,γ)

    isec = findall(==(lon),γ.lon)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    csection= dropdims(view(c,isec,:,:),dims=1)
    return csection
end

function planview(c,depth,γ)

    isec = findall(==(depth),γ.depth)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    cplan = dropdims(view(c,:,:,isec),dims=3)
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
    function dyeplot
    Plot of dye in ocean
# Arguments
- `lat`: latitude arrays
- `depth`: depth array
- `vals`: lat x depth value array
- `lims`: contour levels
"""
function dyeplot(lat, depth, vals, lims)
    #println("turned off due to matplotlib CI setup issue")
    #calc fignum - based on current number of figures
    figure()
    contourf(lat, depth, vals, lims) 
    gca().set_title("Meridional dye concentration")
end

function nearest_gridpoints(lon::Float64,lat::Float64,depth::Float64,γ::grid)

    # do each dimension separately
    lon = 4.4
    lonval, lonloc = findmin(abs.(γ.lon .- lon))

    gridvals = lonloc
    return gridvals
end

"""
    function depthindex(I) 
    
    Get the k-index (depth level) from the Cartesian index
"""
function depthindex(I) 
    k = Vector{Int64}(undef,length(I))
    [k[n]=I[n][3] for n ∈ 1:length(I)]
    return k
end

"""
    function surfaceindex(I) 
    
    Get the vector-index where depth level == 1 and it is ocean.
"""
function surfaceindex(I)
    Isfc = findall(depthindex(I) .==1)
    return Isfc
end

""" 
    function tracerinit(wet,ltype=Float64)
      initialize tracer field on TMI grid
    perhaps better to have a tracer struct and constructor
# Arguments
- `wet`::BitArray mask of ocean points
# Output
- `d`:: 3d tracer field with NaN on dry points
"""
function tracerinit(wet,ltype=Float64)
    # preallocate
    d = Array{ltype}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    d[wet] .= zero(ltype)
    d[.!wet] .= zero(ltype)/zero(ltype) # NaNs with right type
    return d
end

""" 
    function Γ(tracer2D,γ)
    turn 2D surface field into 3D field with zeroes below surface    
# Arguments
- `tracer2D`:: 2D surface tracer field
- `wet`::BitArray mask of ocean points
# Output
- `tracer3D`:: 3d tracer field with NaN on dry points
"""
function Γ(tracer2D::Matrix{T},wet) where T<: Real
    # preallocate
    tracer3D = Array{T}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer3D[wet] .= zero(T)
    tracer3D[.!wet] .= zero(T)/zero(T)
    tracer3D[:,:,1] = tracer2D
    return tracer3D
end

# """ 
#     function misfit_gridded_data(u,Alu,y,d,Wⁱ,Qⁱ,wet)
#     squared model-data misfit
# # Arguments
# - `u`: controls, 2D surface perturbation
# - `Alu`: LU decomposition of water-mass matrix
# - `y`: observations on grid
# - `d`: model constraints
# - `Wⁱ`: inverse of W weighting matrix for observations
# - `Qⁱ`: inverse of Q weighting matrix for controls
# - `wet`: BitArray ocean mask
# # Output
# - `J`: cost function value
# - `gJ`: derivative of cost function wrt to controls
# """
# function misfit_gridded_data(u::Matrix{T},Alu,y::Array{T,3},d::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},Qⁱ::T,wet::BitArray{3}) where T <: Real
#     # a first guess: observed surface boundary conditions are perfect.
#     # set surface boundary condition to the observations.
#     # below surface = 0 % no internal sinks or sources.
#     ỹ = tracerinit(wet,T)
#     n = tracerinit(wet,T)
#     dJdn = tracerinit(wet,T)
#     dJdd = tracerinit(wet,T)
    
#     gJ = tracerinit(wet[:,:,1],T)
    
#     # first-guess reconstruction of observations
#     Δd = d + Γ(u,wet)
#     ỹ[wet] =  Alu\Δd[wet]
#     n = y .- ỹ
    
#     J = n[wet]'* (Wⁱ * n[wet])
#     J += u[wet[:,:,1]]'* (Qⁱ * u[wet[:,:,1]])

#     dJdn[wet] = 2Wⁱ*n[wet]
#     dJdd[wet] = -( Alu'\dJdn[wet])
    
#     gJ[wet[:,:,1]] = dJdd[:,:,1][wet[:,:,1]]
#     gJ[wet[:,:,1]] += 2Qⁱ*u[wet[:,:,1]]

#     # any way to put J and gradient into one function?
#     return J, gJ
# end

# """ 
#     function misfit_gridded_data(u,Alu,y,d,Wⁱ,Qⁱ,wet)
#     squared model-data misfit
#     controls are a vector input for Optim.jl
# # Arguments
# - `u`: controls, vector format
# - `Alu`: LU decomposition of water-mass matrix
# - `y`: observations on grid
# - `d`: model constraints
# - `Wⁱ`: inverse of W weighting matrix for observations
# - `Qⁱ`: inverse of Q weighting matrix for controls
# - `wet`: BitArray ocean mask
# # Output
# - `J`: cost function of sum of squared misfits
# - `gJ`: derivative of cost function wrt to controls
# """
# function misfit_gridded_data(uvec::Vector{T},Alu,y::Array{T,3},d::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},wet::BitArray{3}) where T <: Real
#     # a first guess: observed surface boundary conditions are perfect.
#     # set surface boundary condition to the observations.
#     # below surface = 0 % no internal sinks or sources.
#     u = tracerinit(wet[:,:,1],T)
#     u[wet[:,:,1]] = uvec

#     ỹ = tracerinit(wet,T)
#     n = tracerinit(wet,T)
#     dJdn = tracerinit(wet,T)
#     dJdd = tracerinit(wet,T)
    
#     # first-guess reconstruction of observations
#     Δd = d + Γ(u,wet)
#     ỹ[wet] =  Alu\Δd[wet]
#     n = y .- ỹ
    
#     J = n[wet]'* (Wⁱ * n[wet])
#     #J += u[wet[:,:,1]]'* (Qⁱ * u[wet[:,:,1]])

#     dJdn[wet] = 2*(Wⁱ*n[wet])
#     dJdd[wet] = -( Alu'\dJdn[wet])
    
#     gJ = dJdd[:,:,1][wet[:,:,1]]
#     #gJ += 2Qⁱ*u[wet[:,:,1]]

#     return J, gJ
# end

# """ 
#     function misfit_gridded_data!(J,gJ,u,Alu,y,d,Wⁱ,Qⁱ,wet)
#     in-place version for computations, no allocation
#     squared model-data misfit
#     controls are a vector input for Optim.jl
# # Arguments
# - `u`: controls, vector format
# - `Alu`: LU decomposition of water-mass matrix
# - `y`: observations on grid
# - `d`: model constraints
# - `Wⁱ`: inverse of W weighting matrix for observations
# - `Qⁱ`: inverse of Q weighting matrix for controls
# - `wet`: BitArray ocean mask
# # Output
# - `J`: cost function of sum of squared misfits
# - `gJ`: derivative of cost function wrt to controls
# """
# function misfit_gridded_data!(J::T,gJ::Vector{T},uvec::Vector{T},Alu,y::Array{T,3},d::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},wet::BitArray{3}) where T <: Real
#     # a first guess: observed surface boundary conditions are perfect.
#     # set surface boundary condition to the observations.
#     # below surface = 0 % no internal sinks or sources.
#     u = tracerinit(wet[:,:,1],T)
#     u[wet[:,:,1]] = uvec

#     ỹ = tracerinit(wet,T)
#     n = tracerinit(wet,T)
    
#     dJdn = tracerinit(wet,T)
#     dJdd = tracerinit(wet,T)
    
#     # first-guess reconstruction of observations
#     Δd = d + Γ(u,wet)
#     ỹ[wet] =  Alu\Δd[wet]
#     n = y .- ỹ
    
#     J = n[wet]'* (Wⁱ * n[wet])
#     #J += u[wet[:,:,1]]'* (Qⁱ * u[wet[:,:,1]])

#     dJdn[wet] = 2*(Wⁱ*n[wet])
#     dJdd[wet] = -( Alu'\ dJdn[wet])
    
#     gJ = dJdd[:,:,1][wet[:,:,1]]
#     #gJ += 2Qⁱ*u[wet[:,:,1]]
    
#     return J, gJ
# end

""" 
    function trackpathways(TMIversion,latbox,lonbox)
    Fraction of water originating from surface box    
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `latbox`: min and max latitude of box
- `lonbox`: min and max longitude of box
# Output
- `c`: fraction of water from surface source
- `γ`: TMI grid
"""
function trackpathways(TMIversion,latbox,lonbox)

    A, Alu, γ = config(TMIversion)

    d = surfacepatch(lonbox,latbox,γ)

    # do matrix inversion to get quantity of dyed water throughout ocean:
    c = tracerinit(γ.wet); # pre-allocate c

    # make methods that make the "wet" index unnecessary
    c[γ.wet] = Alu\d[γ.wet] # presumably equivalent but faster than `c = A\d`

    return c, γ
end

""" 
    function gdriveurl(TMIversion)
    placeholder function to give location (URL) of Google Drive input
    in the future, consider a struct or Dict that describes all TMI versions.
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `url`: location (URL) for download
"""
function gdriveurl(TMIname)
    if TMIname == "TMI_2010_2012_4x4x33"
        #        url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
        url = "https://docs.google.com/uc?export=download&id=1Fn_cY-90_RDbBGh6kV0kpXmsvwdjp1Cd"
    else
        url = nothing
    end
end

""" 
    function regeneratedphosphate(TMIversion)
    Regenerated (i.e., accumulated, remineralized) phosphate
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `PO₄ᴿ`: regenerated phosphate
- `γ`: TMI grid
"""
function regeneratedphosphate(TMIversion)

    A, Alu, γ, inputfile = config(TMIversion)
    ΔPO₄ = readtracer(inputfile,"qpo4")

    # PO₄ᴿ = cumulative regenerated phosphate
    PO₄ᴿ = tracerinit(γ.wet); # pre-allocate 
    PO₄ᴿ[γ.wet] = -(Alu\ΔPO₄[γ.wet])
    return PO₄ᴿ, γ
end

""" 
    function volumefilled(TMIversion)

# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `volume`: global ocean volume filled by a surface region
"""
function volumefilled(TMIversion)

    A, Alu, γ = config(TMIversion)
    
    v = cellvolume(γ)
    area = cellarea(γ)
    
    # effectively take inverse of transpose A matrix.
    dVdd = tracerinit(γ.wet); # pre-allocate c
    dVdd[γ.wet] = Alu'\v[γ.wet]

    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    volume = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(volume,0.0)

    # this step could use a function with γ.I argument
    [volume[I[ii][1],I[ii][2]] = dVdd[I[ii]] / area[I[ii][1],I[ii][2]] for ii ∈ eachindex(I) if I[ii][3] == 1]

    return volume
end

""" 
    function surfaceorigin(TMIversion,loc)

# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `loc`: location (lon,lat,depth) of location of interest
# Output
- `origin`: surface map of fraction of source water for a given location
- `γ`: TMI grid
"""
function surfaceorigin(TMIversion,loc)

    A, Alu, γ = config(TMIversion)

    # Find nearest neighbor on grid
    # set δ = 1 at grid cell of interest
    δ = nearestneighbormask(loc,γ)

    # Include option: find grid coordinate by linear interpolation/extrapolation
    # try Interpolations.jl, need interpolation factors that add up to one
    # a false start to fix this issue 
    #v = cellVolume(γ)
    #itp = interpolate(v)

    #
    dVdδ = tracerinit(γ.wet); # pre-allocate c
    dVdδ[γ.wet] = Alu'\δ[γ.wet]

    # surfaceorigin only exists at sea surface
    origin = view(dVdδ,:,:,1)
    return origin, γ
end

function filterdata(u₀,Alu,y,d₀,W⁻,wet)
#using Distributions, PyPlot, PyCall,
 #   LinearAlgebra,  Zygote, ForwardDiff, Optim

    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    fg(x) = misfit_gridded_data(x,Alu,y,d₀,W⁻,wet)
    F0,G0 = fg(u₀)
    fg!(F,G,x) = misfit_gridded_data!(F,G,x,Alu,y,d₀,W⁻,wet)
    #F1,G1 = fg!(F0,G0,u₀)

    # optimize with Optim.jl
    #optimize(Optim.only_fg!(fg!), u₀, Optim.LBFGS())
    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(),Optim.Options(show_trace=true, iterations = 5))

    # no gradient
    #optimize(misfit,u₀,NelderMead())

    # with gradient
    # G = gmisfit(u₀vec)
    # function g!(G, x)
    #     G = gmisfit(x)
    #     return G
    # end
    # out = optimize(misfit,g!,u₀vec,LBFGS())
    return out    
end

function sample_observations(TMIversion,variable)
    
    A, Alu, γ, inputfile = config(TMIversion)

    # take synthetic observations
    # get observational uncertainty
    # could be better to read indices to form wet mask rather than a sample variable.
    θtrue = readtracer(inputfile,variable)

    #ΔPO₄ = readtracer(inputfile,"qpo4")
    σθ = readtracer(inputfile,"σ"*variable)

    ntrue = tracerinit(γ.wet)
    ntrue[γ.wet] = rand(Normal(),length(σθ[γ.wet])) .* σθ[γ.wet]

    y = θtrue .+ ntrue

    # get cost function (J) based on model misfit
    # here the data-model misfit is weighted by the expected error

    # weighting matrix
    W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    Wⁱ = (1/sum(γ.wet)) .* Diagonal(1 ./σθ[γ.wet].^2)
    return y, Wⁱ
end

""" 
    function misfit_gridded_data(u,Alu,y,d,Wⁱ,wet)
    squared model-data misfit
    controls are a vector input for Optim.jl
# Arguments
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `y`: observations on grid
- `d`: model constraints
- `Wⁱ`: inverse of W weighting matrix for observations
- `wet`: BitArray ocean mask
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
function misfit_gridded_data(uvec::Vector{T},Alu,y::Array{T,3},d::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},wet::BitArray{3}) where T <: Real
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    # below surface = 0 % no internal sinks or sources.
    u = tracerinit(wet[:,:,1],T)
    u[wet[:,:,1]] = uvec

    ỹ = tracerinit(wet,T)
    n = tracerinit(wet,T)
    dJdn = tracerinit(wet,T)
    dJdd = tracerinit(wet,T)
    
    # first-guess reconstruction of observations
    Δd = d + Γ(u,wet)
    ỹ[wet] =  Alu\Δd[wet]
    n = y .- ỹ
    
    J = n[wet]'* (Wⁱ * n[wet])
    #J += u[wet[:,:,1]]'* (Qⁱ * u[wet[:,:,1]])

    dJdn[wet] = 2Wⁱ*n[wet]
    dJdd[wet] = Alu'\dJdn[wet]
    
    gJ = -dJdd[:,:,1][wet[:,:,1]]
    #gJ += 2Qⁱ*u[wet[:,:,1]]

    return J, gJ
end

""" 
    function misfit_gridded_data!(J,gJ,u,Alu,y,d,Wⁱ,Qⁱ,wet)
    in-place version for computations, no allocation
    squared model-data misfit
    controls are a vector input for Optim.jl
# Arguments
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `y`: observations on grid
- `d`: model constraints
- `Wⁱ`: inverse of W weighting matrix for observations
- `wet`: BitArray ocean mask
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
#function misfit_gridded_data!(J,gJ,uvec::Vector{T},Alu,y::Array{T,3},d::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},Qⁱ::T,wet::BitArray{3}) where T <: Real
# eliminate type definitions to help with automatic differentiation
function misfit_gridded_data!(J,gJ,uvec,Alu,y,d,Wⁱ,wet) 
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    # below surface = 0 % no internal sinks or sources.
    T = eltype(uvec)
    u = tracerinit(wet[:,:,1],T)
    u[wet[:,:,1]] = uvec

    ỹ = tracerinit(wet,T)
    n = tracerinit(wet,T)
    
    dJdn = tracerinit(wet,T)
    dJdd = tracerinit(wet,T)
    
    # first-guess reconstruction of observations
    Δd = d + Γ(u,wet)
    
    ỹ[wet] =  Alu\Δd[wet]
    n = y .- ỹ

    if gJ != nothing
        dJdn[wet] = 2Wⁱ*n[wet]
        dJdd[wet] =  Alu' \ dJdn[wet]
    
        gJ .= -dJdd[:,:,1][wet[:,:,1]]
        #gJ .+= 2Qⁱ*u[wet[:,:,1]]
    end

    if J != nothing
        # don't understand why I can't compute J and then return J
        # when I tried that, J was not passed out of function
        # because J is a scalar not a vector?
        return n[wet]'* (Wⁱ * n[wet]) # +  u[wet[:,:,1]]'* (Qⁱ * u[wet[:,:,1]])
    end
    
end


end
