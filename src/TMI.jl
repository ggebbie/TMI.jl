module TMI

using Revise
using LinearAlgebra
using SparseArrays
using NetCDF, Downloads
using GoogleDrive

export configTMI, downloadTMI, vec2fld

struct TMIgrid
    lon::Vector{Real}
    lat::Vector{Real}
    depth::Vector{Real}
    coords::Vector{CartesianIndex{3}}                 
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
    coords = gridcoords(TMIfile)

    # get properties of grid
    lat,lon,depth = gridprops(TMIfile)

    γ = TMIgrid(lon,lat,depth,coords)

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
function gridcoords(file) 
    # make the Cartesian tracer grid
    it = convert(Array{Int,1},ncread(file,"xgrid"))
    jt = convert(Array{Int,1},ncread(file,"ygrid"))
    kt = convert(Array{Int,1},ncread(file,"zgrid"))
    coords = CartesianIndex.(it,jt,kt) # formerly Itmi
    return coords
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

function vec2fld(vector::Array{Float64,1},coords::Array{CartesianIndex{3},1})

    nfield = length(vector)
    
    i = zeros(Int,nfield)
    j = zeros(Int,nfield)
    k = zeros(Int,nfield)

    # fix this: no need to repeat this calculation
    [i[n] = coords[n][1] for n ∈ 1:nfield]
    [j[n] = coords[n][2] for n ∈ 1:nfield]
    [k[n] = coords[n][3] for n ∈ 1:nfield]

    nx = maximum(i)
    ny = maximum(j)
    nz = maximum(k)
    field = zeros(nx,ny,nz)

    #- a comprehension
    [field[coords[n]]=vector[n] for n ∈ 1:nfield]
    return field
end

# function tmifld2vec!(vector::Array{Real,1},field::Array{Real,3},index::Array{CartesianIndex{3},1})

#     nv = length(index)
    
#     #- a comprehension
#     [vector[n] = field[indextmi[n]] for n ∈ 1:nv];
#     return vector
# end

end
