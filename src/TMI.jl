module TMI

using Revise
using LinearAlgebra, SparseArrays, NetCDF, Downloads,
    GoogleDrive, Distances
using PyPlot, PyCall

export config, download,
    vec2fld, fld2vec, kindex, surfaceIndex,
    surfacePatch, section,
    layerthickness, cellArea, cellVolume,
    planview, dyeplot, plotextent, tracerFieldInit,
    updateLinearIndex,
    watermassMatrixXYZ, watermassMatrixZYX,
    linearIndexXYZ

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
    R::Array{Union{Integer,Nothing},3}
#    R::LinearIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}} 
    wet::BitArray{3}
end

"""
    function config(url,inputdir)
# Arguments
- `ncfilename`: NetCDF file name
# Output
- `A`: TMI steady-state water-mass matrix
- `grid`: TMI grid coordinates
"""

function config(url,inputdir)

    TMIfile = inputdir * "/TMI_4deg_2010.nc"
    !isfile(TMIfile) ? download(url,inputdir) : nothing
    ncdata = NetCDF.open(TMIfile)

    # put together the sparse matrix, A
    Azyx = watermassMatrixZYX(TMIfile)
    
    c = readTracer(TMIfile,"θ")

    # make a mask
    wet = .!isnan.(c)
    
    # Do bookkeeping for TMI grid coordinates
    I = cartesianIndexZYX(TMIfile)

    R = linearIndexXYZ(wet)

    Axyz = watermassMatrixXYZ(TMIfile,I,R)
    #R = R[ wet ] # eliminate land points

    # LU factorization for efficient matrix inversions
    Alu = lu(Axyz)
    
    # get properties of grid
    lat,lon,depth = gridProps(TMIfile)

    γ = grid(lon,lat,depth,I,R,wet)

    return Azyx, Axyz, Alu, c, γ
end

"""
    function cartesianIndexZYX(file)
    Read and assemble the grid coordinates
    according to the legacy MATLAB code.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function cartesianIndexZYX(file)
    # make the Cartesian tracer grid
    it = convert(Array{Int,1},ncread(file,"xgrid"))
    jt = convert(Array{Int,1},ncread(file,"ygrid"))
    kt = convert(Array{Int,1},ncread(file,"zgrid"))
    I = CartesianIndex.(it,jt,kt)
    return I
end

"""
    function linearIndexXYZ(file)
    Read and assemble the grid coordinates.
# Arguments
- `wet`: 3D mask for wet points
# Output
- `R`: array of linear indices, but not a LinearIndices type
"""
function linearIndexXYZ(wet)
    R = Array{Union{Int64,Nothing},3}(nothing,size(wet))
    R[wet]=1:sum(wet)
    # R = LinearIndices((1:maximum(it),1:maximum(jt),1:maximum(kt)));
    # R = LinearIndices((it,jt,kt));
    #Rwet = R[γ.wet]
    return R
end

"""
    function gridProps(file)
    Read and assemble the grid properties.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function gridProps(file)
    lat = ncread(file,"lat")
    lon = ncread(file,"lon")
    depth = ncread(file,"depth")
    return lat,lon,depth
end

"""
    function watermassMatrixZYX(file)
    Read and assemble the water-mass matrix.
    Legacy version from MATLAB.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `A`: water-mass matrix
"""
function watermassMatrixZYX(file)
    i = ncread(file,"i")
    j = ncread(file,"j")
    m = ncread(file,"m")
    A = sparse(i,j,m)
    return A
end

"""
        function watermassMatrixXYZ(file,γ)
    Read and assemble the water-mass matrix from MATLAB.
    Transfer to updated x,y,z version
# Arguments
- `file`: TMI NetCDF file name
- `γ`: TMI grid
# Output
- `A`: water-mass matrix
"""
function watermassMatrixXYZ(file,I,R)

    # MATLAB accounting z,y,x
    izyx = convert(Vector{Int},ncread(file,"i"))
    jzyx = convert(Vector{Int},ncread(file,"j"))

    # Julia accounting x,y,z
    ixyz = updateLinearIndex(izyx,I,R)
    jxyz = updateLinearIndex(jzyx,I,R)
    
    m = ncread(file,"m")

    # use grid indices to switch i,j values
    A = sparse(ixyz,jxyz,m)
    #A = A[wet,:];
    return A
end


"""
    function updateLinearIndex(i,I,R)
    Linear index translated from z,y,x to x,y,z accounting
# Arguments
- `izyx`: index of interest in z,y,x accounting
- `I`: wet Cartesian Index for z,y,x
- `R`: Linear indices for x,y,z 
# Output
- `ixyz`: index of interest in x,y,z accounting
"""
function updateLinearIndex(izyx,I,R)
    #ixyz = similar(izyx)
    #[ixyz[vv] = R[I[izyx[vv]]] for vv in eachindex(izyx)]
    ixyz = R[I[izyx]]
    return ixyz
end

"""
    function watermassMatrix(file)
    Read and assemble the water-mass matrix.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `A`: water-mass matrix
"""
function watermassMatrix(file)

    m = massFractions(file)

    A = watermassMatrix(m)
    
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
    function watermassMatrix(file)
    Assemble the water-mass matrix given mass fractions `m`
# Arguments
- `m`: mass fractions
# Output
- `A`: water-mass matrix
"""
# function watermassMatrix(m)
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
    function readTracer(file,tracername)
    Read and assemble the water-mass matrix.
# Arguments
- `file`: TMI NetCDF file name
- `tracername`: name of tracer
# Output
- `c`: 3D tracer field
"""
function readTracer(file,tracername)
    c = ncread(file,tracername)
    return c
end

function cellArea(γ)
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

function cellVolume(γ)
    dz = layerthickness(γ)
    area = cellArea(γ)
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
    function surfacePatch
    Make a surface boundary condition
    with a rectangular patch
# Arguments
- `lonbox`: longitudes of box edges
- `latbox`: latitudes of box edges
- `γ`: TMI.grid
# Output
- `d`: vector that describes surface patch
"""
function surfacePatch(lonbox::Vector{T},latbox::Vector{T},γ::grid)::Array{Float64} where T<:Real

    # ternary operator to handle longitudinal wraparound
    lonbox[1] ≤ 0 ? lonbox[1] += 360 : nothing
    lonbox[2] ≤ 0 ? lonbox[2] += 360 : nothing

    # define the surface boundary condition

    # preallocate
    d = tracerFieldInit(γ.wet)

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
    function kindex(I) 
    
    Get the k-index (depth level) from the Cartesian index
"""
function kindex(I) 
    k = Vector{Int64}(undef,length(I))
    [k[n]=I[n][3] for n ∈ 1:length(I)]
    return k
end

"""
    function surfaceIndex(I) 
    
    Get the vector-index where depth level == 1 and it is ocean.
"""
function surfaceIndex(I)
    Isfc = findall(kindex(I) .==1)
    return Isfc
end

""" 
    function tracerFieldInit(wet)
      initialize tracer field on TMI grid
    perhaps better to have a tracer struct and constructor
"""
function tracerFieldInit(wet)
    # preallocate
    d = Array{Float64}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    d[wet] .= 0.0
    d[.!wet] .= NaN
    return d
end

end
