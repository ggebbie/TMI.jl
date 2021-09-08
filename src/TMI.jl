module TMI

using Revise
using LinearAlgebra
using SparseArrays
using NetCDF, Downloads
using GoogleDrive

export configTMI, downloadTMI, vec2fld, surfacepatch, section

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
    Transfer a vector to a 3D field
"""
function vec2fld(vector::Array{Float64,1},I::Array{CartesianIndex{3},1})

    nfield = length(vector)

    nx = maximum(I)[1]
    ny = maximum(I)[2]
    nz = maximum(I)[3]
    field = zeros(nx,ny,nz)

    # a comprehension
    [field[I[n]]=vector[n] for n ∈ eachindex(I)]
    return field
end

# function tmifld2vec!(vector::Array{Real,1},field::Array{Real,3},index::Array{CartesianIndex{3},1})

#     nv = length(index)
    
#     #- a comprehension
#     [vector[n] = field[indextmi[n]] for n ∈ 1:nv];
#     return vector
# end

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

end
