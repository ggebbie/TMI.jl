module TMI

using Revise
using LinearAlgebra, SparseArrays, NetCDF, Downloads,
    GoogleDrive, Distances, GibbsSeaWater,  
    PyPlot, PyCall, Distributions, Optim,
    Interpolations, LineSearches, MAT, NCDatasets,
    OrdinaryDiffEq, PreallocationTools

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
    costfunction, costfunction!,
    trackpathways, regeneratedphosphate, meanage,
    volumefilled, surfaceorigin, synthetic_observations,
    observe,
    steadyclimatology, steadyinversion,
    interpweights, interpindex,
    wetlocation, iswet,
    control2state, control2state!,
    surfacecontrol2field, surfacecontrol2field!,
    sparsedatamap, config2nc, gridprops,
    matrix_zyx2xyz, varying!, readopt, ces_ncwrite,
    surface_oxygensaturation, oxygen, location_obs,
    getsurfaceboundary, zerosurfaceboundary,
    onesurfaceboundary, setboundarycondition!,
    adjustboundarycondition, adjustboundarycondition!,
    gsetboundarycondition, setsource!,
    zeros, ones, maximum, minimum, (+), (-), (*), dot,
    Grid, Field, BoundaryCondition, vec, unvec!, unvec

    #pkgdir, pkgdatadir, pkgsrcdir, not needed?

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
    struct Tracer

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
    Can array have other element types?
    Are indices needed?
"""
struct BoundaryCondition{T}
    tracer::Array{T,2}
    #I::Vector{CartesianIndex{3}} # index
    #R::Array{Int,3}
    dim::Int64
    dimval::Int64
    wet::BitArray{2}
    #T::DataType
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

"""
    function config_from_nc(TMIversion)
    Configure TMI environment from NetCDF input file.
# Arguments
- `TMIversion`: TMI version for water-mass/circulation model
# Output
- `A`: TMI steady-state water-mass matrix
- `Alu`: LU decomposition of A
- `γ`: TMI grid properties
- `TMIfile`: TMI file name
"""
function config_from_nc(TMIversion)

    #make datdir() if it doesn't exist 
    !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 
    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    #if TMIfile doesn't exist, get GDrive url and download 
    if !isfile(TMIfile)

        # add a workaround for large files
        if TMIversion == "modern_180x90x33_GH11_GH12"
            println("workaround for 2° x 2°")
            shellscript = pkgsrcdir("read_nc_modern_180x90x33_GH11_GH12.sh")
            run(`sh $shellscript`)
            mv(joinpath(pwd(),"TMI_"*TMIversion*".nc"),TMIfile)
        else
            println("read via GoogleDrive.jl")
            #- `url`: Google Drive URL for data
            url = ncurl(TMIversion)
            google_download(url,pkgdatadir())
        end

    end 
    
    # get properties of grid
    lon,lat,depth = gridprops(TMIfile)

    # read Cartesian Index from file.
    I = cartesianindex(TMIfile)

    # make a mask
    # first choice: smaller but inconsistent with input grid
    #wet = falses(maximum(I)[1],maximum(I)[2],maximum(I)[3])
    wet = falses(length(lon),length(lat),length(depth))
    wet[I] .= 1

    R = linearindex(wet)

    println("A")
    @time A = watermassmatrix(TMIfile)

    # LU factorization for efficient matrix inversions
    println("Alu")
    @time Alu = lu(A)
    
    γ = Grid(lon,lat,depth,I,R,wet)

    # would be good to make this optional
    println("L=")
    @time L = circulationmatrix(TMIfile,A,γ)

    println("B=")
    @time B = boundarymatrix(TMIfile,γ)
    
    return  A, Alu, γ, TMIfile, L, B

end

"""
Configure TMI environment from original MATLAB output
"""
function config_from_mat(TMIversion)

    #- `url`: Google Drive URL for data
    url = maturl(TMIversion)
    TMIfile = pkgdatadir("TMI_"*TMIversion*".mat")
    TMIfilegz = TMIfile*".gz"
    println(TMIfile)
    !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 
    #    !isfile(TMIfilegz) & !isfile(TMIfile) ? google_download(url,pkgdatadir()) : nothing

    if TMIversion == "modern_180x90x33_GH11_GH12"
        println("workaround for 2° x 2°")
        shellscript = pkgsrcdir("read_mat_modern_180x90x33_GH11_GH12.sh")
        run(`sh $shellscript`)
        mv(joinpath(pwd(),"TMI_"*TMIversion*".mat.gz"),TMIfilegz,force=true)
    else
        !isfile(TMIfilegz) & !isfile(TMIfile) && google_download(url,pkgdatadir())
    end
    
    # cloak mat file in gz to get Google Drive spam filter to shut down
    isfile(TMIfilegz) & !isfile(TMIfile) && run(`gunzip $TMIfilegz`) 
    
    # # make a sample field from zyx cartesian indices
    Izyx = cartesianindex(TMIfile)

    # # make a mask
    wet = falses(maximum(Izyx)[1],maximum(Izyx)[2],maximum(Izyx)[3])
    wet[Izyx] .= 1

    I = cartesianindex(wet)

    R = linearindex(wet)

    Azyx = watermassmatrix(TMIfile)
    A = matrix_zyx2xyz(TMIfile,Azyx,R)

    # # LU factorization for efficient matrix inversions
    Alu = lu(A)

    # get properties of grid
    lat,lon,depth = gridprops(TMIfile)

    γ = Grid(lon,lat,depth,I,R,wet)

    # need to make this optional
    L = circulationmatrix(TMIfile,γ)
    
    B = boundarymatrix(TMIfile,γ)
    
    # consider re-ordering this.
    # some output should be optional
    # return Izyx or I or neither?
    #return  A, Alu, γ, TMIfile, I, L, B
    return  A, Alu, γ, TMIfile, L, B
end

"""
    function cartesianindex(file)
    Read and assemble the grid coordinates
    according to the legacy MATLAB code (z,y,x order).
# Arguments
- `file`: TMI NetCDF file name
# Output
- `I`: TMI Cartesian index for wet points
"""
function cartesianindex(file::String)
    # make the Cartesian tracer grid
    if file[end-1:end] == "nc"

        it = convert(Vector{Int},ncread(file,"i"))
        jt = convert(Vector{Int},ncread(file,"j"))
        kt = convert(Vector{Int},ncread(file,"k"))

        I = CartesianIndex.(it,jt,kt)

    elseif file[end-2:end] == "mat"

        matobj = matopen(file)
        haskey(matobj,"it") ? it=convert(Vector{Integer},vec(read(matobj,"it"))) : it = convert(Vector{Integer},vec(read(matobj,"i")))
        haskey(matobj,"jt") ? jt=convert(Vector{Integer},vec(read(matobj,"jt"))) : jt=convert(Vector{Integer},vec(read(matobj,"j")))
        haskey(matobj,"kt") ? kt=convert(Vector{Integer},vec(read(matobj,"kt"))) : kt=convert(Vector{Integer},vec(read(matobj,"k")))
        close(matobj)
        I = CartesianIndex.(it,jt,kt) 
    end
    return I
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
    function gridprops(file)
    Read and assemble the grid properties.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function gridprops(file)
    if file[end-1:end] == "nc"
        
        lon = convert(Vector{Float64},ncread(file,"lon"))
        lat = convert(Vector{Float64},ncread(file,"lat"))
        depth = convert(Vector{Float64},ncread(file,"depth"))

    elseif file[end-2:end] == "mat"
        
        matobj = matopen(file)
        lon=convert(Vector{Float64},vec(read(matobj,"LON")))
        lat=convert(Vector{Float64},vec(read(matobj,"LAT")))
        depth=convert(Vector{Float64},vec(read(matobj,"DEPTH")))
        close(matobj)

    end
    
    return lon,lat,depth
end

"""
    function watermassmatrix(file)
    Read and assemble the water-mass matrix.
# Arguments
- `file`: TMI NetCDF or MATLAB file name
# Output
- `A`: water-mass matrix
"""
function watermassmatrix(file)

    # consider adding a catch if A doesn't exist in file.

    if file[end-1:end] == "nc"
        # Int or Integer?
        i = convert(Vector{Int},ncread(file,"Arow"))
        j = convert(Vector{Int},ncread(file,"Acol"))
        m = ncread(file,"m")
        A = sparse(i,j,m)
    elseif file[end-2:end] == "mat"
        matobj = matopen(file)
        if haskey(matobj,"A")
            A=read(matobj,"A")
            close(matobj)
        else
            close(matobj)
            return nothing
        end

        # But MATLAB had zyx format and we need xyz format.
        # linearindices R not available so will do conversion in higher scope

    end
    return A
end

"""
        function matrix_zyx2xyz(TMIfile,Azyx,γ)
   
    Transfer zyx format water-mass matrix A to xyz format
# Arguments
- `Azyx`: water-mass matrix in zyx format
- `γ`: TMI grid
# Output
- `Axyz`: water-mass matrix in xyz format
"""
function matrix_zyx2xyz(file,Azyx,R)

    izyx, jzyx, mzyx = findnz(Azyx)
    Izyx = cartesianindex(file)
        
    # Julia accounting x,y,z
    ixyz = updatelinearindex(izyx,Izyx,R)
    jxyz = updatelinearindex(jzyx,Izyx,R)
    
    # use grid indices to switch i,j values
    Axyz = sparse(ixyz,jxyz,mzyx)
    return Axyz
end

"""
    function circulationmatrix(file,γ)
    Read and assemble the circulation matrix from MATLAB.
    Transfer to updated x,y,z version
# Arguments
- `file`: TMI MATLAB file name
- `γ`: TMI grid
# Output
- `L`: circulation matrix in xyz format
"""
function circulationmatrix(file,γ)

    if file[end-2:end] == "mat" 

        matobj = matopen(file)
        if haskey(matobj,"L")
            # Matlab output in zyx format
            Lzyx=read(matobj,"L")
            close(matobj)

            Izyx = cartesianindex(file)
            izyx, jzyx, Fzyx = findnz(Lzyx)
            # Julia accounting x,y,z
            ixyz = updatelinearindex(izyx,Izyx,γ.R)
            jxyz = updatelinearindex(jzyx,Izyx,γ.R)
            L = sparse(ixyz,jxyz,Fzyx)

        else
            close(matobj)
            return nothing
        end

    elseif file[end-1:end] == "nc"

        # based on function arguments, read from inefficient storage of L matrix.
        if haskey(NCDataset(file),"F")

            i = convert(Vector{Int},ncread(file,"Lrow"))
            j = convert(Vector{Int},ncread(file,"Lcol"))
            F = ncread(file,"F")
            L = sparse(i,j,F)
        else
            return nothing
        end
    end
    
    return L
end

"""
    function circulationmatrix(file,A,γ)
    Read and assemble the circulation matrix from the efficient storage of A and F₀ variables. 
# Arguments
- `file`: TMI MATLAB file name
- `A`: TMI water-mass matrix
- `γ`: TMI grid
# Output
- `L`: circulation matrix in xyz format
"""
function circulationmatrix(file,A,γ)

    file[end-1:end] !== "nc" && error("not a NetCDF file")

    # based on function arguments, read F₀ to efficiently reproduce L matrix.

    if haskey(NCDataset(file),"F₀")
        F₀ = ncread(file,"F₀")
        F₀vec = F₀[γ.wet]
    
        # For each row of A, multiply by F₀
        i, j, F = findnz(A)

        # careful, this loop can be really slow
        for nn in eachindex(i)
            F[nn] *= F₀vec[i[nn]]
        end

        L = sparse(i,j,F)
        return L
    else
        return nothing
    end
    
end

"""
        function boundarymatrix(file,γ)
    Read and assemble the boundary matrix from MATLAB.
    Transfer to updated x,y,z version
# Arguments
- `file`: TMI MATLAB file name
- `γ`: TMI grid
# Output
- `B`: boundary condition matrix
"""
function boundarymatrix(file,γ)

    if file[end-2:end] == "mat"

        matobj = matopen(file)
        if haskey(matobj,"B")
            Bzyx=read(matobj,"B")
            close(matobj)

            # matlab in zyx format.
            # consider using Azyx2xyz here.
            Izyx = cartesianindex(file)
            izyx, jzyx, Fzyx = findnz(Bzyx)
            # for B, rows are 3D grid space, columns are for the surface index. 
            # Julia accounting x,y,z
            Isfc = surfaceindex(Izyx)
            ixyz = updatelinearindex(izyx,Izyx,γ.R)
            jxyz = updatelinearindex(Isfc[jzyx],Izyx,γ.R)

            # assume surface at k = 1 (revisit for LGM problem)
            # give the full dimension of sparse matrix
            B = sparse(ixyz,jxyz,Fzyx,sum(γ.wet),sum(γ.wet[:,:,1]))
        else
            close(matobj)
            return nothing
        end

    elseif file[end-1:end] == "nc"

        # based on function arguments, read from inefficient storage of L matrix.
        if haskey(NCDataset(file),"b")

            i = convert(Vector{Int},ncread(file,"Brow"))
            j = convert(Vector{Int},ncread(file,"Bcol"))
            b = ncread(file,"b")
            B = sparse(i,j,b,sum(γ.wet),sum(γ.wet[:,:,1]))

        else
            return nothing
        end
    end
    return B
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
    b = BoundaryCondition(patch,3,1,γ.wet[:,:,1])
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

# """
# Find the locations of the control adjustments
# Is this necessary?
# """
# function controlindex(plane)
#     i = lonindex(I)
#     j = latindex(I)
#     Isfc = findall(depthindex(I) .==1 .|| i.==maximum(i) .|| i.==minimum(i) .|| j.==maximum(j) .|| j.==minimum(j))
#     return Isfc
# end

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
   Initialize boundary condition with zeroes
"""
function zeros(dim::Int64,dimval::Int64,wet::BitArray{3})::BoundaryCondition

    dimsize = size(wet)
    # dumb way to do it
    if dim == 1
        wet2d = wet[dimval,:,:]
    elseif dim == 2
        wet2d = wet[:,dimval,:]
    elseif dim == 3
        wet2d = wet[:,:,dimval]
    else
        error("boundary condition not implemented in 4+ dimensions")
    end
    
    tracer = Array{Float64}(undef,size(wet2d))
    tracer[wet2d] .= zero(Float64)
    tracer[.!wet2d] .= zero(Float64)/zero(Float64)
    
    b = BoundaryCondition(tracer,dim,dimval,wet2d)

end

"""
   Initialize boundary condition with ones
"""
function ones(dim::Int64,dimval::Int64,wet::BitArray{3})::BoundaryCondition

    dimsize = size(wet)
    # dumb way to do it
    if dim == 1
        wet2d = wet[dimval,:,:]
    elseif dim == 2
        wet2d = wet[:,dimval,:]
    elseif dim == 3
        wet2d = wet[:,:,dimval]
    else
        error("boundary condition not implemented in 4+ dimensions")
    end
    
    tracer = Array{Float64}(undef,size(wet2d))
    tracer[wet2d] .= ones(Float64)
    tracer[.!wet2d] .= zero(Float64)/zero(Float64)
    
    b = BoundaryCondition(tracer,dim,dimval,wet2d)

end

"""
   Get boundary condition by extracting from 3D tracer
"""
function getboundarycondition(tracer3d,dim,dimval,wet)::BoundaryCondition

    dimsize = size(wet)
    # dumb way to do it
    if dim == 1
        wet2d = wet[dimval,:,:]
        tracer2d = tracer3d[dimval,:,:]
    elseif dim == 2
        wet2d = wet[:,dimval,:]
        tracer2d = tracer3d[:,dimval,:]
    elseif dim == 3
        wet2d = wet[:,:,dimval]
        tracer2d = tracer3d[:,:,dimval]
    else
        error("boundary condition not implemented in 4+ dimensions")
    end
    
    b = BoundaryCondition(tracer2d,dim,dimval,wet2d)

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
    e = BoundaryCondition(array,c.dim,c.dimval,c.wet)
    
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
    e = BoundaryCondition(array,d.dim,d.dimval,d.wet)
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


# define the correct dimension and index for each control plane
# maybe someday find a way to hide γ
zerosurfaceboundary(γ) = zeros(3,1,γ.wet)::BoundaryCondition
zeronorthboundary(γ) = zeros(2,maximum(latindex(γ.I)),γ.wet)::BoundaryCondition
zeroeastboundary(γ) = zeros(1,maximum(lonindex(γ.I)),γ.wet)::BoundaryCondition
zerosouthboundary(γ) = zeros(2,1,γ.wet)::BoundaryCondition
zerowestboundary(γ) = zeros(1,1,γ.wet)::BoundaryCondition

onesurfaceboundary(γ) = ones(3,1,γ.wet)::BoundaryCondition

getsurfaceboundary(c::Field) = getboundarycondition(c.tracer,3,1,c.γ.wet)::BoundaryCondition
getnorthboundary(c::Field) = getboundarycondition(c.tracer,2,maximum(latindex(c.γ.I)),c.γ.wet)::BoundaryCondition
geteastboundary(c::Field) = getboundarycondition(c.tracer,1,maximum(lonindex(c.γ.I)),c.γ.wet)::BoundaryCondition
getsouthboundary(c::Field) = getboundarycondition(c.tracer,2,1,c.γ.wet)::BoundaryCondition
getwestboundary(c::Field) = getboundarycondition(c.tracer,1,1,c.γ.wet)::BoundaryCondition

### TURN THIS INTO SET BOUNDARY CONDITION COMMAND.
# """ 
#     function constraint(b::BoundaryCondition,γ)
#     turn control adjustment into adjustment of constraint for all 3d grid points    
# # Arguments
# - `b`:: BoundaryCondition
# - `γ`:: TMI.Grid
# # Output
# - `d`:: right hand side adjustment
# """
# function constraint(b::BoundaryCondition,γ::Grid) 
#     # preallocate
#     T = eltype(b.tracer)
#     d = Array{T}(undef,size(γ.wet))

#     # set ocean to zero, land to NaN
#     # consider whether land should be nothing or missing
#     d[γ.wet]   .= zero(T)
#     d[.!γ.wet] .= zero(T)/zero(T)
#     if b.dim == 1
#         d[b.dimval,:,:] = b.tracer
#     elseif b.dim == 2
#         d[:,b.dimval,:] = b.tracer
#     elseif b.dim == 3
#         d[:,:,b.dimval] = b.tracer
#     else
#         error("controls not implemented for 4+ dimensions")
#     end
#     return d
# end

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
    gb = 0.0 * b # initialize to zero
    if b.dim == 1
        gb.tracer = gd.tracer[b.dimval,:,:]
    elseif b.dim == 2
        gb.tracer = gd.tracer[:,b.dimval,:]
    elseif b.dim == 3
        gb.tracer .+= gd.tracer[:,:,b.dimval] 
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
    u = BoundaryCondition(tracer,utemplate.dim,utemplate.dimval,utemplate.wet)
    return u
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
    #u.tracer[b.wet] += u.tracer[u.wet] # write it out so b changes when returned
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
    function control2state(tracer2D,γ)
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
    function ncurl(TMIversion)
    placeholder function to give location (URL) of NetCDF Google Drive input
    in the future, consider a struct or Dict that describes all TMI versions.
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `url`: location (URL) for download
"""
function ncurl(TMIname)
    if TMIname == "modern_90x45x33_GH10_GH12"
        url = "https://docs.google.com/uc?export=download&id=1Fn_cY-90_RDbBGh6kV0kpXmsvwdjp1Cd"
    elseif TMIname == "modern_180x90x33_GH11_GH12"
        url = "https://docs.google.com/uc?export=download&id=1-YEkB_YeQGqPRH6kauhBb2bi_BjVGt9b"
    elseif TMIname == "modern_90x45x33_unpub12"
        url = "https://docs.google.com/uc?export=download&id=1Kw_Mr7fiKqan0nx0dKvGHnSInP0hQ7AV"
    elseif TMIname == "modern_90x45x33_G14"
        url = "https://docs.google.com/uc?export=download&id=1aeE7EXA-vy3Cm_drt4qCFw4AlpYrdudk"
    elseif TMIname == "modern_90x45x33_G14_v2"
        url = "https://docs.google.com/uc?export=download&id=1Mwhv70soBX6-pYijU0ElNl0TZw0vSbXN"
    elseif TMIname == "LGM_90x45x33_G14"
        url = "https://docs.google.com/uc?export=download&id=1yoDi7_foBt3TVULCstlWnNLHFc2G47Fz"  
    elseif TMIname == "LGM_90x45x33_G14A"
        url = "https://docs.google.com/uc?export=download&id=1ADkDI3Fc3z4Vm75K5u6hx0Yu1P0iVnW1"
    elseif TMIname == "LGM_90x45x33_GPLS1"
        url = "https://docs.google.com/uc?export=download&id=1VOrZGUsO7lp21qw6Yw0dBqGlIy_ImdRW"
    elseif TMIname == "LGM_90x45x33_GPLS2"
        url = "https://docs.google.com/uc?export=download&id=1cOCrty9kvA2s3NoD1QZjnNehlVbP0rHP"
    elseif TMIname == "LGM_90x45x33_OG18"
        url = "https://docs.google.com/uc?export=download&id=19zccG1BSdspD9rti2OttsF2Dm4P2OLjt"
    else
        url = nothing
    end
end

""" 
    function maturl(TMIversion)
    Find *mat file here.
    placeholder function to give location (URL) of Google Drive input
    in the future, consider a struct or Dict that describes all TMI versions.
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `url`: location (URL) for download
"""
function maturl(TMIname)
    if TMIname == "modern_90x45x33_GH10_GH12"
        url = "https://docs.google.com/uc?export=download&id=1Z2knDctAmZHO2lcWTBCdR8zjkmbcyCGg"
    elseif TMIname == "modern_180x90x33_GH11_GH12"
        url = "https://docs.google.com/uc?export=download&id=11zD1nOfT6V7G0qIHdjK2pDGHFk-ExXwU"
    elseif TMIname == "modern_90x45x33_unpub12"
        url = "https://docs.google.com/uc?export=download&id=1sqkjFCPxZT_2Bm9rsp0acyxxkBri9YAT"
    elseif TMIname == "modern_90x45x33_G14"
        url = "https://docs.google.com/uc?export=download&id=1dCrDe5VXrsXiOf04mbuHID7xc5Ymm8-z"
    elseif TMIname == "modern_90x45x33_G14_v2"
        url = "https://docs.google.com/uc?export=download&id=1Axaqn88HZv3i4rYn0g5dUztVXlx_TmGt"
    elseif TMIname == "LGM_90x45x33_G14"
                url = "https://docs.google.com/uc?export=download&id=1qWnL9SrcjFRt4TUQ0k7v0jEFrKtUYpev"
    elseif TMIname == "LGM_90x45x33_G14A"
                url = "https://docs.google.com/uc?export=download&id=1E7qYWdAnoz4YXFvJ0AkVcqpHDq8aFP4R"
    elseif TMIname == "LGM_90x45x33_GPLS1"
                url = "https://docs.google.com/uc?export=download&id=1nF5KjWE3n--NkOTgIymIhP6pkafJ7jl3"
    elseif TMIname == "LGM_90x45x33_GPLS2"
                url = "https://docs.google.com/uc?export=download&id=1CHiR2J60HPRrxXTvJupyWEHmggMkgnxV"
    elseif TMIname == "LGM_90x45x33_OG18"
                url = "https://docs.google.com/uc?export=download&id=1WKAWxDYDFls6C5YzO1sMLDlUhC8_tKV2"
    else
        url = nothing
    end
    return url
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

    ∂V∂b  = BoundaryCondition(volume,3,1,γ.wet[:,:,1])
    
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
    origin = BoundaryCondition(dvlocdd,3,1,γ.wet[:,:,1])
    
    return origin
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

    # Handle longitudinal periodic condition (i.e., wraparound)

    lon = vcat(copy(γ.lon),γ.lon[1]+360.)
    list = vcat(1:length(γ.lon),1)
    nodes = (lon,γ.lat,γ.depth)

    # eliminate need to pass tracer value
    wis = Interpolations.weightedindexes((Interpolations.value_weights,),((Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))),nodes, loc)

    # issue, some of weighted points may be NaNs in tracer field
    # handle this in the Interpolations.jl routines
    # may involve chaging Gridded(Linear()) above
    return wis
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
    function sparsedatamap(u₀::Vector{T},Alu,b::BoundaryCondition{T},y::Vector{T},W⁻,wis,locs,Q⁻,γ::Grid;iterations=10) where T <: Real

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
function sparsedatamap(u₀::Vector{T},Alu,b::BoundaryCondition{T},y::Vector{T},W⁻,wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}}},locs,Q⁻,γ::Grid,iterations=10) where T <: Real

     fg!(F,G,x) = costfunction!(F,G,x,Alu,b,y,W⁻,wis,locs,Q⁻,γ)
    
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(linesearch = LineSearches.BackTracking()),Optim.Options(show_trace=true, iterations = iterations))

    return out    
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
function synthetic_observations(TMIversion,variable,γ,N)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    
    θtrue = readfield(TMIfile,variable,γ)
    replace!(θtrue.tracer,NaN=>0.0)
    
    σθ = readfield(TMIfile,"σ"*variable,γ)
    replace!(σθ.tracer,NaN=>0.0)

    # get random locations that are wet (ocean)
    locs = Vector{Tuple{Float64,Float64,Float64}}(undef,N)
    [locs[i] = wetlocation(γ) for i in eachindex(locs)]

    # get weighted interpolation indices
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    ytrue = observe(θtrue,wis,γ)
    σtrue = observe(σθ,wis,γ)

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
function costfunction_gridded_obs(uvec,Alu,b::BoundaryCondition{T},utemplate::BoundaryCondition{T},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {T <: Real}

    # turn uvec into a boundary condition
    u = unvec(utemplate,uvec)

    bnew = adjustboundarycondition(b,u) #b += u # easy case where u and b are on the same boundary
    println("max b ",maximum(b))
    n = steadyinversion(Alu,bnew,γ) - y  # gives the misfit
    J = n ⋅ (Wⁱ * n) # dot product

    # adjoint equations
    gy = -2Wⁱ * n
    gb = gsteadyinversion( gy, Alu, b, γ)
    gu = gadjustboundarycondition(gb,u)
    guvec = vec(gu)

    return J, guvec
end

function costfunction_gridded_obs!(J,guvec,uvec::Vector{T},Alu,b::BoundaryCondition{T},utemplate::BoundaryCondition{T},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real

    # turn uvec into a boundary condition
    u = unvec(utemplate,uvec)

    bnew = adjustboundarycondition(b,u) #b += u # easy c
    #adjustboundarycondition!(b,u) # easy case where u and b are on the same boundary
    y -= steadyinversion(Alu,bnew,γ)  # gives the misfit

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
    function costfunction_gridded_obs(uvec,Alu,b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}
"""
function costfunction_gridded_obs(uvec,Alu,b₀::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},utemplate::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}

    # turn uvec into a boundary condition
    u = unvec(utemplate,uvec)

    b = adjustboundarycondition(b₀,u) #b += u # easy case where u and b are on the same boundary
    y -= steadyinversion(Alu,b,γ)  # gives the misfit
    J = y ⋅ (Wⁱ * y) # dot product

    # adjoint equations
    gy = -2Wⁱ * y
    gb = gsteadyinversion( gy, Alu, b, γ)
    gu = gadjustboundarycondition(gb,u)
    guvec = vec(gu)

    return J, guvec
end

function costfunction_gridded_obs!(J,guvec,uvec::Vector{T},Alu,b₀::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u₀::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}

    # turn uvec into a boundary condition
    u = unvec(u₀,uvec)
    
    b = adjustboundarycondition(b₀,u) #b += u # easy case where u and b are on the same boundary
    y -= steadyinversion(Alu,b,γ)  # gives the misfit

    if guvec != nothing
        # adjoint equations
        gy = -2Wⁱ * y
        gb = gsteadyinversion( gy, Alu, b, γ)
        gu = gadjustboundarycondition(gb,u)
        tmp = vec(gu)
        for (ii,vv) in enumerate(tmp)
            guvec[ii] = vv
        end
    end
    
    if J !=nothing
        return  y ⋅ (Wⁱ * y) # dot product
    end
end


# """
# function costfunction_obs(uvec,Alu,b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where {N1, N2, T <: Real}

#     squared model-data misfit for pointwise data
#     controls are a vector input for Optim.jl
#     Issue: couldn't figure out how to nest with costfunction_obs!
#     Issue: why are wis and locs both needed? `gobserve` function
# """
# function costfunction_obs(uvec,Alu,b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where {N1, N2, T <: Real}

#     # data misfit and gradient

#     unvec!(u,uvec) # need to know which controls are present

#     adjustboundarycondition!(b,u) # map u -> b
#     c = steadyinversion(Alu,b,γ)  # gives the misfit

#     # observe at right spots
#     ỹ = observe(c,wis,γ)
#     n = ỹ - y

#     Jcontrol = uvec'*(Q⁻*uvec)
#     Jdata = n ⋅ (Wⁱ * n) # dot product
#     J = Jdata + Jcontrol

#     guvectmp = 2*(Q⁻*uvec)
#     gn = 2Wⁱ * n
#     gỹ = gn

#     gc = TMI.gobserve(gỹ,c,locs)
#     gb = gsteadyinversion(gc, Alu, b, γ)

#     # gadjust_boundary_condition!
#     #gu = gb

#     # for ii in 1:sum(gu.wet)
#     #     # force guvec to change?
#     #     guvec[ii] = gu.tracer[u.wet][ii]
#     #     guvec[ii] += guvectmp[ii]
#     # end
#     #end
#     guvec = []
#     return J, guvec

# end

""" 
    function costfunction(J,gJ,uvec,Alu,b,y,Wⁱ,wis,Q⁻,γ)
    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue #1: couldn't figure out how to nest with costfunction_obs!
    Issue #2: Update for BoundaryCondition types
    
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
function costfunction(uvec::Vector{T},Alu,b::BoundaryCondition{T},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where T <: Real

    # control penalty and gradient
    Jcontrol = uvec'*(Q⁻*uvec)
    guvec = 2*(Q⁻*uvec)

    # data misfit and gradient
    u = zerosurfaceboundary(γ)
    u.tracer[u.wet] = uvec
    
    b += u # easy case where u and b are on the same boundary
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
    gu = gb

    guvec .+= gu.tracer[gu.wet]
    
    return J, guvec
end

""" 
    function costfunction!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,Q⁻,γ)
    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue: couldn't figure out how to nest with costfunction_obs!
    Issue: why are wis and locs both needed? `gobserve` function

"""
function costfunction!(J,guvec,uvec::Vector{T},Alu,b::BoundaryCondition{T},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where T <: Real

    # data misfit and gradient
    u = zerosurfaceboundary(γ)
    u.tracer[u.wet] = uvec
    
    b += u # easy case where u and b are on the same boundary
    c = steadyinversion(Alu,b,γ)  # gives the misfit

    # observe at right spots
    ỹ = observe(c,wis,γ)
    n = ỹ - y

    if guvec != nothing    
        guvectmp = 2*(Q⁻*uvec)
        gn = 2Wⁱ * n

        gỹ = gn
        
        gc = gobserve(gỹ,c,locs)

        gb = gsteadyinversion(gc, Alu, b, γ)
        gu = gb 

        for ii in 1:sum(gu.wet)
            # force guvec to change?
            guvec[ii] = gu.tracer[u.wet][ii]
            guvec[ii] += guvectmp[ii]
        end
    end

    if J != nothing
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
        loc = (rand(0.0:0.1:360.0),
               rand(-90.0:0.1:90.0),
               rand(0.0:1.0:5750.0))

        iswet(loc,γ) && return loc
        println("dry point, try again")
    end # if not, then start over.
end

"""
    function varying!(du, u, p, t)
    ODE function for varying boundary cond
    Sets up dc/dt = L*C + B*f to be solved 
# Arguments
- `du`: dc/dt (must have this name for DifferentialEquations.jl to work
- `u`: C, what we are solving for 
- `p`: parameters for diffeq - must hold specified vars  
- `t`: time we are solving for (automatically determined by DE.jl)
# Output
- `du`: numerical value of LC+Bf, vector of size 74064 for 4°
"""
function varying!(du, u, p, t)
    
    #load parameters 
    Csfc,surface_ind,τ,L,B,li,LC,BF,Cb = p
    println("time = ", t)
    
    #generate Cb - interpolated surface boundary condition  
    li_t = convert(Float64, li[t])
    Cb .= (ceil(li_t)-li_t).*Csfc[Int(floor(li_t)), :] .+ (li_t-floor(li_t)).*Csfc[Int(ceil(li_t)), :]
    
    #use PreallocationTools.jl to handle Dual type in u 
    LC = get_tmp(LC, first(u)*t) 
    BF = get_tmp(BF, first(u)*t) 

    #Figure out what u is at surface 
    u_sfc = @view u[surface_ind]

    #Inplace math to make faster 
    mul!(LC, L, u) 
    mul!(BF, B, -(u_sfc.-Cb)./τ) 
    @. du = LC + BF 
    nothing
end

"""
    function ces_ncwrite(γ,time,sol_array)
    Write .nc file output for commonerasim.jl 
# Arguments
- `γ`: 
- `time`: vector of time values 
- `sol_array`: solution array in form time x lat x lon x depth - must match γ + time 
# Output
- saves .nc file titled "ces_output.nc" in data array 
"""
function ces_ncwrite(γ,time,sol_array)
    file = pkgdatadir() * "/ces_output.nc"
    ds = NCDataset(file,"c")

    #define dimensions 
    defDim(ds,"lon", size(γ.lon)[1])
    defDim(ds,"lat",size(γ.lat)[1])
    defDim(ds,"depth",size(γ.depth)[1])
    defDim(ds,"time",size(time)[1])

    #write theta output variable 
    v = defVar(ds,"theta",Float64, ("time","lon","lat","depth"))
    v[:,:,:,:] = sol_array

    #write dimensions as variables (this might not be kosher...) 
    vlon = defVar(ds,"lon",Float64, ("lon",))
    vlon[:] = γ.lon
    vlat = defVar(ds,"lat",Float64,("lat",))
    vlat[:] = γ.lat
    vtime = defVar(ds,"time",Float64,("time",))
    vtime[:] = time
    vdepth = defVar(ds,"depth",Float64,("depth",))
    vdepth[:] = γ.depth

    v.attrib["title"] = "output of commonerasim.jl" 
    v.attrib["units"] = "potential temperature anomaly"
    close(ds)
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

""" 
Save TMI configuration to NetCDF format for non-proprietary access
"""
function config2nc(TMIversion,A,γ,L,B)

    # make new netcdf file.
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    isfile(filenetcdf) && rm(filenetcdf)

    grid2nc(TMIversion,γ)
    
    matfields2nc(TMIversion,γ)

    !isnothing(A) && watermassmatrix2nc(TMIversion,A)

    !isnothing(L) && circulationmatrix2nc(TMIversion,L,γ)

    !isnothing(B) && boundarymatrix2nc(TMIversion,B)

    #= is this part of the config? Or should it go to
     a separate output? It is similar to the output fields above. Probably should be considered part of the config. =#
    regions2nc(TMIversion,γ)

    optim2nc(TMIversion)

end

"""
Save grid dictionaries of attributes for writing to NetCDF file
"""
function griddicts(γ)
    # update names and types in dictionary
    
    TMIgrids = Dict("lon" => γ.lon,
                    "lat" => γ.lat,
                    "depth" => γ.depth)
    
    TMIgridsatts = Dict("lon" => Dict("longname" => "Longitude", "units" => "°E"),
                        "lat" => Dict("longname" => "Latitude", "units" => "°N"),
                        "depth" => Dict("longname" => "depth", "units" => "m"))

    return TMIgrids, TMIgridsatts

end

"""
Read 3D fields from mat file and save to NetCDF file.
"""
function matfields2nc(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    filemat = pkgdatadir("TMI_"*TMIversion*".mat")
    vars = matread(filemat)

    TMIgrids, TMIgridsatts = griddicts(γ)

    T = eltype(γ.lon) # does the eltype of longitude have to equal the tracer eltype?
    #T =  Float64

    varlist = Dict("dP"=>"qPO₄","q"=>"qPO₄",
                   "Tobs"=>"θ","Tmod"=>"θ","Tlgm"=>"θ",
                   "Terr"=>"σθ",
                   "Sobs"=>"Sp","Smod"=>"Sp","Slgm"=>"Sp",
                   "Serr"=>"σSp",
                   "O18obs"=>"δ¹⁸Ow","O18mod"=>"δ¹⁸Ow","O18lgm"=>"δ¹⁸Ow",
                   "O18err"=>"σδ¹⁸Ow",
                   "Pobs"=>"PO₄","Pmod"=>"PO₄","Plgm"=>"PO₄",
                   "Perr" => "σPO₄",
                   "Nobs"=>"NO₃","Nmod"=>"NO₃","Nlgm"=>"NO₃",
                   "Nerr" => "σNO₃",
                   "Oobs"=>"O₂","Omod"=>"O₂","Olgm"=>"O₂",
                   "Oerr"=>"σO₂",
                   "C13obs"=>"δ¹³C","C13mod"=>"δ¹³C","C13lgm"=>"δ¹³C",
                   "C13err" =>  "σδ¹³C")

    # iterate over all possible variables listed above
    Izyx = cartesianindex(filemat)
    TMIfields = Dict{String,Array{T,3}}()
    for (kk,vv) in varlist
        haskey(vars,kk) ? push!(TMIfields, vv => tracerinit(vars[kk], Izyx, γ.wet)) : nothing
    end

    # also save fields that are stored in the x struct, if they exist
    if haskey(vars,"x")
        for (kk,vv) in varlist
            haskey(vars["x"],kk) ? push!(TMIfields, vv => tracerinit(vars["x"][kk], Izyx, γ.wet)) : nothing
        end
    end
    
    TMIfieldsatts = fieldsatts()

    # iterate in TMIgrids Dictionary to write to NetCDF.
    for (varname,varvals) in TMIfields
        
        nccreate(filenetcdf,varname,"lon",γ.lon,TMIgridsatts["lon"],"lat",γ.lat,TMIgridsatts["lat"],"depth",γ.depth,TMIgridsatts["depth"],atts=TMIfieldsatts[varname])
        println("write ",varname)
        ncwrite(varvals,filenetcdf,varname)

    end
end

"""
All variable names and attributes.
Useful for writing NetCDF files.
"""
fieldsatts() = 
    Dict("θ" => Dict("longname" => "potential temperature", "units" => "°C"),
         "σθ" => Dict("longname" => "1σ standard error in potential temperature", "units" => "°C"),
         "Sp" => Dict("longname" => "practical salinity", "units" => "PSS-78"),
         "σSp" => Dict("longname" => "1σ standard error in practical salinity", "units" => "PSS-78"),
         "δ¹⁸Ow" => Dict("longname" => "oxygen-18 to oxygen-16 ratio in seawater", "units" => "‰ VSMOW"),
         "σδ¹⁸Ow" => Dict("longname" => "1σ standard error in oxygen-18 to oxygen-16 ratio in seawater", "units" => "‰ VSMOW"),
         "PO₄" => Dict("longname" => "phosphate", "units" => "μmol/kg"),
         "σPO₄" => Dict("longname" => "1σ standard error in phosphate", "units" => "μmol/kg"),
         "qPO₄" => Dict("longname" => "local source of phosphate", "units" => "μmol/kg"),
         "NO₃" => Dict("longname" => "nitrate", "units" => "μmol/kg"),
         "σNO₃" => Dict("longname" => "1σ standard error in nitrate", "units" => "μmol/kg"),
         "O₂" => Dict("longname" => "dissolved oxygen", "units" => "μmol/kg"),
         "σO₂" => Dict("longname" => "1σ standard error in dissolved oxygen", "units" => "μmol/kg"),
         "δ¹³C" => Dict("longname" => "carbon-13 to carbon-12 ratio in DIC", "units" => "‰ PDB"),
         "σδ¹³C" => Dict("longname" => "1σ standard error fin carbon-13 to carbon-12 ratio in DIC", "units" => "‰ PDB"),
         "F₀" => Dict("longname" => "normalized mass flux out of gridcell", "units" => "(kg seawater/s)/(kg gridcell)"))


"""
Read vectors from mat file, translate to 3D,
 and save surface field to NetCDF file.
"""
function regions2nc(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    filemat = pkgdatadir("TMI_"*TMIversion*".mat")

    # region names
    # didn't figure out how to use an ordered dict, instead use a tuple
    list = ("GLOBAL","ANT","SUBANT",
            "NATL","NPAC","TROP","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

    regionname = Dict("GLOBAL" => "globally uniform",
                      "ANT" => "Antarctic",
                      "SUBANT" => "Subantarctic",
                      "NATL" => "North Atlantic",
                      "NPAC" => "North Pacific",
                      "TROP" => "tropical and subtropical",
                      "ARC" => "Arctic",
                      "MED" => "Mediterranean",
                      "ROSS" => "Ross Sea sector",
                      "WED" => "Weddell Sea sector",
                      "LAB" => "Labrador and Irminger Seas",
                      "GIN" => "Greenland-Iceland-Norwegian Seas",
                      "ADEL" => "Adélie Land sector",
                      "SUBANTATL" => "Atlantic-sector Subantarctic",
                      "SUBANTPAC" => "Pacific-sector Subantarctic",
                      "SUBANTIND" => "Indian-sector Subantarctic",
                      "TROPATL" => "tropical and subtropical Atlantic",
                      "TROPPAC" => "tropical and subtropical Pacific",
                      "TROPIND" => "tropical and subtropical Indian")
    
    matobj = matopen(filemat)
    if haskey(matobj,"d_all")
        d_all = read(matobj,"d_all")
        close(matobj)
    else
        return
    end

    # a kludge for now
    T = eltype(γ.lon)
    
    # iterate over all regions in d_all
    Izyx = cartesianindex(filemat)
    regions = Dict{String,Array{T,2}}()
    regionatts = Dict{String,Dict{String,String}}()
    
    for rr = 1:size(d_all,2)
        # 3D fields in zyx vector format
        # are changed to 3D xyz format
        d = tracerinit(d_all[:,rr],Izyx,γ.wet)

        # just save the surface 2D field
        push!(regions, list[rr] => d[:,:,1])
        
        push!(regionatts, list[rr] =>
         Dict("longname" => regionname[list[rr]]*" surface region", "units" => "[]"))
    end

    TMIgrids, TMIgridsatts = griddicts(γ)

    # iterate in regions Dictionary to write to NetCDF.
    for (varname,varvals) in regions
        dvarname = "d_"*varname
        nccreate(filenetcdf,dvarname,"lon",γ.lon,TMIgridsatts["lon"],"lat",γ.lat,TMIgridsatts["lat"],atts=regionatts[varname])
        println("write ",dvarname)
        ncwrite(varvals,filenetcdf,dvarname)

    end
end

function watermassmatrix2nc(TMIversion,A)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    i, j, m = findnz(A)
    nelements = length(i)

    # add the circulation matrix: problem can't store sparse matrix.
    varname= "m"
    elementatts = Dict("longname" => "TMI sparse matrix element number")
    matts =  Dict("longname" => "TMI water-mass fraction (sparse matrix values)", "units" => "(kg source)/(kg total)")
    nccreate(filenetcdf,varname,"A_element",1:nelements,elementatts,atts=matts)
    println("write ",varname)
    ncwrite(m, filenetcdf,varname)
    
     varname= "Arow"
     destatts = Dict("longname" => "gridcell number of destination (row value)")
     nccreate(filenetcdf,varname,"A_element",1:nelements,elementatts,atts=destatts)
    println("write ",varname)
    ncwrite(i, filenetcdf,varname)

    varname= "Acol"
    sourceatts = Dict("longname" => "gridcell number of source (column value)")
    nccreate(filenetcdf,varname,"A_element",1:nelements,elementatts,atts=sourceatts)
    println("write ",varname)
    ncwrite(j, filenetcdf,varname)

end

"""
Save optimization parameters to NetCDF file)

Future considerations: split into 2 functions
1) read from mat
2) save to nc
"""
function optim2nc(TMIversion)

    filemat = pkgdatadir("TMI_"*TMIversion*".mat")
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")

    matobj = matopen(filemat)
    if haskey(matobj,"fval")
        J̃ = read(matobj,"fval")
        iteratts = Dict("longname" => "iteration number")
        Jatts =  Dict("longname" => "cost function value", "units" => "[]")

        varname = "J" # J̃ not output to screen properly
        println("write ",varname)
        nccreate(filenetcdf,varname,"iter",1:length(J̃),iteratts,atts=Jatts)
    end
    if haskey(matobj,"u")
        ũ = read(matobj,"u")
        iteratts = Dict("longname" => "control element number")
        uatts =  Dict("longname" => "control vector", "units" => "[]")

        varname = "ũ"
        println("write ",varname)
        nccreate(filenetcdf,varname,"control_element",1:length(ũ),iteratts,atts=uatts)
    end
    close(matobj)

end

"""
Save circulation matrix `L` to NetCDF file.
"""
function circulationmatrix2nc(TMIversion,L,γ)

    T = eltype(L)
    fullmatrix = false # more efficient to just save F₀, then modify A to get L 
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    if !fullmatrix
        F₀ = tracerinit(γ.wet,T)
        for nd ∈ eachindex(γ.I)
            # normalized mass flux out of gridcell is found on diagonal
            F₀[γ.I[nd]] = -L[nd,nd]
        end

        TMIgrids, TMIgridsatts = griddicts(γ)
        TMIfieldsatts = fieldsatts()
        varname = "F₀"
        nccreate(filenetcdf,varname,"lon",γ.lon,TMIgridsatts["lon"],"lat",γ.lat,TMIgridsatts["lat"],"depth",γ.depth,TMIgridsatts["depth"],atts=TMIfieldsatts[varname])
        println("write ",varname)
        ncwrite(F₀,filenetcdf,varname)

    else # fullmatrix = true, makes a bigger nc file

        i, j, F = findnz(L)
        nelements = length(i)

        # add the circulation matrix: problem can't store sparse matrix.
        varname = "F"
        elementatts = Dict("longname" => "TMI sparse matrix element number")
        Fatts =  Dict("longname" => "normalized mass flux (sparse matrix values)","units" => "(kg seawater/s)/(kg gridcell)")
        nccreate(filenetcdf,varname,"L_element",1:nelements,elementatts,atts=Fatts)
        println("write ",varname)
        ncwrite(F, filenetcdf,varname)
        
        varname= "Lrow"
        destatts = Dict("longname" => "gridcell number of destination (row value)")
        nccreate(filenetcdf,varname,"L_element",1:nelements,elementatts,atts=destatts)
        println("write ",varname)
        ncwrite(i, filenetcdf,varname)

        varname= "Lcol"
        sourceatts = Dict("longname" => "gridcell number of source (column value)")
        nccreate(filenetcdf,varname,"L_element",1:nelements,elementatts,atts=sourceatts)
        println("write ",varname)
        ncwrite(j, filenetcdf,varname)
    end
    
end

"""
Save boundary matrix for transient model to NetCDF file
"""
function boundarymatrix2nc(TMIversion,B)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    i, j, b = findnz(B)
    nelements = length(i)

    # add the circulation matrix: problem can't store sparse matrix.
    varname= "b"
    elementatts = Dict("longname" => "TMI boundary matrix element number")
    matts =  Dict("longname" => "Boundary matrix values", "units" => "[]")
    nccreate(filenetcdf,varname,"B_element",1:nelements,elementatts,atts=matts)
    println("write ",varname)
    ncwrite(b, filenetcdf,varname)

    varname= "Brow"
    destatts = Dict("longname" => "gridcell number of 3D field")
    nccreate(filenetcdf,varname,"B_element",1:nelements,elementatts,atts=destatts)
    println("write ",varname)
    ncwrite(i, filenetcdf,varname)

    varname= "Bcol"
    sourceatts = Dict("longname" => "gridcell number of surface boundary condition")
    nccreate(filenetcdf,varname,"B_element",1:nelements,elementatts,atts=sourceatts)
    println("write ",varname)
    ncwrite(j, filenetcdf,varname)

end

"""
Put grid properties (Cartesian index) into NetCDF file
"""
function grid2nc(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")

    linearindexatts = Dict("longname" => "linear index")
    nfld = length(γ.I)

    iatts =  Dict("longname" => "Cartesian index for x-direction", "units" => "[]")
    jatts =  Dict("longname" => "Cartesian index for y-direction", "units" => "[]")
    katts =  Dict("longname" => "Cartesian index for z-direction", "units" => "[]")

    varname = "i"
    nccreate(filenetcdf,varname,"linearindex",1:nfld,linearindexatts,atts=iatts)
    ncwrite(lonindex(γ.I),filenetcdf,varname)
    
    varname = "j"
    nccreate(filenetcdf,varname,"linearindex",1:nfld,linearindexatts,atts=jatts)
    ncwrite(latindex(γ.I),filenetcdf,varname)
    
    varname = "k"
    nccreate(filenetcdf,varname,"linearindex",1:nfld,linearindexatts,atts=katts)
    ncwrite(depthindex(γ.I),filenetcdf,varname)
    
end

"""
     function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    Read an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)
    Return a BoundaryCondition
"""
function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    file = pkgdatadir("TMI_"*TMIversion*".nc")
    tracername = "d_"*region

    # Didn't use readfiled because dsfc is 2d.
    dsfc = ncread(file,tracername)
    b = BoundaryCondition(dsfc,3,1,γ.wet[:,:,1])
    return b
end

#read surface layer
function readopt(filename,γ)
    nc = NCDataset(filename)
#    lat = nc["latitude"][:]
#    lon = nc["longitude"][:]
    time = nc["year"][:]
#    depth = nc["depth"][:]
    theta = nc["theta"][:, :, :, :]

    #flip to time descending order
    reverse!(time)
    time = convert(Vector{Int}, time)
    reverse!(theta, dims = 1)
    theta_permuted = zeros((size(theta)[1], size(γ.wet)[1], size(γ.wet)[2], size(γ.wet)[3]))
    permutedims!(theta_permuted, theta, [1,4,3,2])
    return time, theta_permuted 
end

end
