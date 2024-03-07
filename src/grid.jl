using Base: index_ndims, methods_including_ambiguous
"""
    struct Grid

    TMI grid with accounting for wet/dry points
"""
struct Grid{T,N}
    # lon::Vector{T}
    # lat::Vector{T}
    # depth::Vector{T}
    axes::NTuple{N,Vector{T}}
    wet::BitArray{N}
    interior::BitArray{N}
end

"""
function Grid(TMIfile)

    Construct the Grid given a file name

# Arguments
- `TMIfile::String`: NetCDF file name for TMI version

# Output
- `γ::Grid`: TMI grid struct
"""
function Grid(TMIfile::String; A = watermassmatrix(TMIfile))

    # get properties of grid
    labels = axislabels(TMIfile)

    ndims = length(labels)
    #dimsize = Tuple{Int64}(undef,ndims)
    dimsize = Tuple([length(labels[d]) for d in 1:ndims])
    #end
    
    # make ocean mask
    #wet = wetmask(TMIfile,length(lon),length(lat),length(depth))
    wet = wetmask(TMIfile,dimsize)

    # make interior mask
    #interior = interiormask(A,wet,length(lon),length(lat),length(depth))
    interior = interiormask(A,wet,dimsize)
    return Grid(labels,wet,interior)
end

"""
function Grid(foreign_file, maskname, lonname, latname, depthname)

    Construct the Grid from a non-TMI file given the names of relevant fields.

    Assumes that an ocean mask is available.
    Assumes an input NetCDF file.
    Assumes everything below the top layer is part of the interior. 
    Tested for Float32 fields (should work for other types).

# Arguments
- `foreign_file::String`
- `maskname::String`
- `lonname::String`
- `latname::String`
- `depthname::String`
# Output
- `γ::Grid`: TMI grid struct
"""
function Grid(foreign_file::S, maskname::S, lonname::S, latname::S, depthname::S) where S <: String
    
    # make ocean mask
    ds = Dataset(foreign_file)
    wet = Bool.(ds[maskname])::BitArray # very slow! (couple of secs), use `convert` instead?

    T = eltype(ds[lonname][1]) # not sure that this will always work
    
    lon = convert(Vector{T},ds[lonname]) 
    lat = convert(Vector{T},ds[latname])
    depth = - convert(Vector{T},ds[depthname]) # flip sign for actual "depth"

    # make interior mask: Assume no lateral or bottom boundaries (CAUTION)
    interior = deepcopy(wet)
    interior[:,:,1] .= false

    return Grid((lon,lat,depth),wet,interior)
end

"""
    Base.propertynames(γ::Grid) = (I,R,fieldnames(typeof(x))...)

    Do not store Cartesian and linear indices.
    Compute them on demand.
""" 
Base.propertynames(γ::Grid) = (:I,:R,:lon,:lat,:depth,fieldnames(typeof(γ))...)
function Base.getproperty(γ::Grid, d::Symbol)
    if d === :I
        return cartesianindex(γ.wet)
    elseif d === :R
        return linearindex(γ.wet)
    elseif d === :lon
        return γ.axes[1]
    elseif d === :lat
        return γ.axes[2]
    elseif d === :depth
        return γ.axes[3]
    else
        return getfield(γ, d)
    end
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
    R[wet]=1:sum(wet)
    return R
end

"""
    function checkgrid!(c,wet)

    perform a check of file compatibility
     with grid
"""
function checkgrid!(tracer,wet)
    if sum(isnan.(tracer[wet])) > 0
        println(sum(isnan.(tracer[wet]))," points")
        error("readfield warning: NaN on grid")
    end

    # check for non NaN or nonzero off grid
    # Need to rethink how to do this.
    if sum(isnan.(tracer[.!(wet)])) < length(isnan.(tracer[.!(wet)]))
        println("readfield warning: non-NaN value off grid")
        println("resetting to NaN")
        tracer[.!(wet)] .= NaN
    end
end

function wetmask(TMIfile::String,dimsize::NTuple)
    # read Cartesian Index from file.
    I = cartesianindex(TMIfile)
    wet = falses(dimsize)
    wet[I] .= 1
    return wet
end
function wetmask(TMIfile::String,nx,ny,nz)
    # older version that assumes N = 3
    # read Cartesian Index from file.
    I = cartesianindex(TMIfile)
    wet = falses(nx,ny,nz)
    wet[I] .= 1
    return wet
end

"""
function interiormask(A,wet,nx,ny,nz)
"""
function interiormask(A,wet,nx,ny,nz)
    interior = falses(nx,ny,nz)
    I = cartesianindex(wet)
    list = findall(.!(isone.(sum(abs.(A),dims=2))))
    interior[I[list]] .= true 
    return interior
end
function interiormask(A,wet,dimsize::NTuple)
    interior = falses(dimsize)
    I = cartesianindex(wet)
    list = findall(.!(isone.(sum(abs.(A),dims=2))))
    interior[I[list]] .= true 
    return interior
end
