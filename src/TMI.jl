module TMI

using Revise
using LinearAlgebra, SparseArrays, NetCDF, Downloads,
    GoogleDrive, Distances
#using PyPlot, PyCall

export configTMI, downloadTMI, vec2fld, fld2vec, surfacepatch, section
export layerthickness, cellarea, cellvolume, planview
#export dyeplot, plotextent 

#Python packages - initialize them to null globally
#const patch = PyNULL()
#const ccrs = PyNULL()

#Initialize all Python packages - install with conda through Julia
# function __init__()
#     copy!(patch, pyimport_conda("matplotlib.patches", "patches"))
#     copy!(ccrs, pyimport_conda("cartopy.crs", "ccrs"))
#     print("Python libraries installed")
# end

struct grid
    lon::Vector{Real}
    lat::Vector{Real}
    depth::Vector{Real}
    I::Vector{CartesianIndex{3}} # index
    #nx::Int
    #ny::Int
    #nz::Int
end

"""
    function readTMI(ncfilename)
# Arguments
- `ncfilename`: NetCDF file name
# Output
- `A`: TMI steady-state water-mass matrix
- `grid`: TMI grid coordinates
"""
function configTMI(url,inputdir)

    TMIfile = inputdir * "/TMI_4deg_2010.nc"
    !isfile(TMIfile) ? downloadTMI(url,inputdir) : nothing
    ncdata = NetCDF.open(TMIfile)

    # put together the sparse matrix, A
    A = watermassmatrix(TMIfile)

    # LU factorization for efficient matrix inversions
    Alu = lu(A)

    # Do bookkeeping for TMI grid coordinates
    I = gridindex(TMIfile)

    # get properties of grid
    lat,lon,depth = gridprops(TMIfile)

    γ = grid(lon,lat,depth,I)

    return A, Alu, γ
end

"""
    function gridcoords(file)
    Read and assemble the grid coordinates.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function gridindex(file)
    # make the Cartesian tracer grid
    it = convert(Array{Int,1},ncread(file,"xgrid"))
    jt = convert(Array{Int,1},ncread(file,"ygrid"))
    kt = convert(Array{Int,1},ncread(file,"zgrid"))
    I = CartesianIndex.(it,jt,kt)
    return I
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
    function watermassmatrix(file)
    Read and assemble the water-mass matrix.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `A`: water-mass matrix
"""
function watermassmatrix(file)
    i = ncread(file,"i")
    j = ncread(file,"j")
    m = ncread(file,"m")
    A = sparse(i,j,m)
    return A
end

function cellarea(γ)
    dx = zonalgriddist(γ)
    dy = haversine((γ.lon[1],γ.lat[1])
                  ,(γ.lon[1],γ.lat[2]))

    area = Vector{Float64}(undef,length(γ.I))
    [area[v] = dx[γ.I[v][2]] for v ∈ eachindex(γ.I)]
    area *= dy
    return area
end

function cellvolume(γ)
    dz = layerthickness(γ)
    area = cellarea(γ)
    volume = similar(area)
    [volume[v] = dz[γ.I[v][3]] for v ∈ eachindex(γ.I)]
    volume .*= area
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
    function downloadTMI(url,inputdir)
    Read and assemble all TMI inputs.
# Arguments
- `url`: Google Drive location of TMI input
- `inputdir`: input directory location to store file
# Output
- none
"""
function downloadTMI(url,inputdir)
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
function surfacepatch(lonbox,latbox,γ)

    # ternary operator to handle longitudinal wraparound
    lonbox[1] ≤ 0 ? lonbox[1] += 360 : nothing
    lonbox[2] ≤ 0 ? lonbox[2] += 360 : nothing

    # define the surface boundary condition
    nfield = length(γ.I) # number of ocean points
    d = zeros(Int,nfield) # preallocate
    [d[n]=1 for n ∈ 1:nfield if γ.I[n][3]==1 && latbox[1] ≤ γ.lat[γ.I[n][2]] ≤ latbox[2]
         && lonbox[1] ≤ γ.lon[γ.I[n][1]] ≤ lonbox[2] ]
    return d
end

"""
    function section
    View latitude-depth slice of field
# Arguments
- `cfld`: 3d tracer field
- `lon`: longitude of section
- `γ`: TMI.grid
# Output
- `csection`: 2d slice of field
"""
function section(cfld,lon,γ)

    isec = findall(==(lon),γ.lon)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    csection= dropdims(view(cfld,isec,:,:),dims=1)
    return csection
end

function planview(cfld,depth,γ)

    isec = findall(==(depth),γ.depth)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    cplan = dropdims(view(cfld,:,:,isec),dims=3)
    return cplan
end


"""
    function plotextent
    Generate image showing user-specified ROI
# Arguments
- `latbox`: in format [lat_start, lat_stop]
- `lonbox`: in format [lon_start, lon_stop]

"""
# function plotextent(latbox, lonbox)
    
#     ccrs = pyimport("cartopy.crs")
#     lower_left = [minimum(lonbox), minimum(latbox)] #array of lower left of box

#     #calc width and height of box
#     w = maximum(lonbox) - minimum(lonbox)
#     h = maximum(latbox) - minimum(latbox)

#     #init GeoAxes
#     fig = figure()
#     ax = fig.add_subplot(projection = ccrs.PlateCarree())

#     #plot rectangle
#     ax.add_patch(patch.Rectangle(xy=lower_left,
#                                  width=w, height=h,
#                                  facecolor="blue",
#                                  alpha=0.2,
#                                  transform=ccrs.PlateCarree()))
#     #define extent of figure
#     pad = 10 #how many deg lat and lon to show outside of bbox
#     pad_add = [-pad, pad] #add this to latbox and lonbox
#     padded_lat = latbox + pad_add
#     padded_lon = lonbox + pad_add
#     ext = vcat(padded_lon, padded_lat) #make into one vector
#     ax.set_extent(ext)

#     # using cartopy 0.18 and NaturalEarth is missing
#     # ax.coastlines() #show coastlines

#     #add gridlines
#     gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
#     gl.top_labels = false
#     gl.right_labels = false

#     ax.set_title("User-defined surface patch")
# end

# """
#     function dyeplot
#     Plot of dye in ocean
# # Arguments
# - `lat`: latitude arrays
# - `depth`: depth array
# - `vals`: lat x depth value array
# - `lims`: contour levels
# """
# function dyeplot(lat, depth, vals, lims)
#     #calc fignum - based on current number of figures
#     figure()
#     contourf(lat, depth, vals, lims) 
#     gca().set_title("Meridional dye concentration")
# end

function nearest_gridpoints(lon::Float64,lat::Float64,depth::Float64,γ::grid)

    # do each dimension separately
    lon = 4.4
    lonval, lonloc = findmin(abs.(γ.lon .- lon))

    gridvals = lonloc
    return gridvals
end

end
