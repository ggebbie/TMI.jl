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
struct BoundaryCondition{T <: Real,R <: Real, N <: Integer, B <: AbstractMatrix{T}}
    tracer::B
    i::Vector{R}
    j::Vector{R}
    k::R
    dim::N
    dimval::N
    wet::BitMatrix
    name::Symbol
    longname::String
    units::String
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, x::BoundaryCondition)
    summary(io, x); println(io)
    print(io, "Field size ")
    println(io, size(x.tracer))
    show(io,mime,heatmap(transpose(x.tracer),zlabel=x.units,title=x.longname))
end

"""
function BoundaryCondition(tracer::AbstractMatrix{T},i::Vector{R},j::Vector{R},k::R,dim::N,dimval::N,wet::BitMatrix) where T <: Real where R <: Real where N <: Integer

    Outer constructor for BoundaryCondition if there's no worry about
    tracer type, long name, or units.
# Arguments
- `tracer::AbstractMatrix{T}`
- `i::Vector{Real}`
- `j::Vector{Real}`
- `k::Real`
- `dim::Integer`
- `dimval:Integer`
- `γ::Grid`
# Output
- `b::BoundaryCondition`
"""
# an outer constructor that ignores units
function BoundaryCondition(tracer::AbstractMatrix{T},i::Vector{R},j::Vector{R},k::R,dim::N,dimval::N,γ::Grid) where T <: Real where R <: Real where N <: Integer

    return BoundaryCondition(tracer,i,j,k,dim,dimval,γ,:bc,"boundary condition","unknown") 
end

"""
    function writeboundarycondition(file,b,gamma::Grid)

    Write a BoundaryCondition to NetCDF.
 
    Use NCDatasets so that Unicode is correct

# Arguments
- `file`: TMI NetCDF file name
- `b::BoundaryCondition`: a TMI.BoundaryCondition struct
# Output
- none
# Side-effect
- write to `file`
"""
function write(file,b::BoundaryCondition{T,R,N,B}) where T <: Real where R <: Real where N <: Integer where B <: AbstractMatrix

    if T == Bool
        Tcheck = Int8
    else
        Tcheck = T
    end

    # only 2 dimensions are needed
    if b.dim == 1
        lat = b.i
        depth = b.j
        lon = b.k
        tuple2d = ("lat","depth")
    elseif b.dim == 2
        lon = b.i
        depth = b.j
        lat = b.k
        tuple2d = ("lon","depth")
    elseif b.dim == 3
        lon = b.i
        lat = b.j
        depth = b.k
        tuple2d = ("lon","lat")
    end
    Nx = length(lon)
    Ny = length(lat)
    Nz = length(depth)

    if !isfile(file)
        # create new NetCDF file
        ds = Dataset(file,"c")

        atts = TMI.gridatts()
        ds.attrib["title"] = "boundary condition"
        
        defDim(ds,"lon",Nx)
        defDim(ds,"lat",Ny)
        defDim(ds,"depth",Nz) 
        
        vlon = defVar(ds,"lon",R,["lon"],
            attrib = OrderedDict(atts["lon"]))
        vlon[:] = lon

        vlat = defVar(ds,"lat",R,["lat"],
            attrib = OrderedDict(atts["lat"]))
        vlat[:] = lat

        vdepth = defVar(ds,"depth",R,["depth"],
            attrib = OrderedDict(atts["depth"]))
        vdepth[:] = depth
        
        v = defVar(ds,String(b.name),Tcheck,tuple2d,
            attrib = OrderedDict("longname" => b.longname,
                "units" => b.units))
        v[:,:] = b.tracer

        close(ds)

    else
        # assumption: on the same grid
        ds = Dataset(file,"a")

        println(b.name)
        v = defVar(ds,String(b.name),Tcheck,tuple2d,
            attrib = OrderedDict("longname" => b.longname,
                      "units" => b.units))
        v[:,:] = b.tracer
        close(ds)
    end
        
    return nothing
end

"""
    function boundaryconditionatts(dim::Int64,dimval::Int64,γ::Grid)

       Help initialize boundary condition by getting some attributes
"""
function boundaryconditionatts(dim::Integer,dimval::Integer,γ::Grid)

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
    function zeros(dim::Int64,dimval::Int64,γ::Grid,name::Symbol,longname::String,units::String)::BoundaryCondition

       Initialize boundary condition with zeroes
# Arguments
- `dim`:
- `dimval`
- `γ::Grid`
- `name::Symbol`
- `longname::String`
- `units::String`

# Output
- `b::BoundaryCondition`
"""
function zeros(dim::I,dimval::I,γ::Grid,name::Symbol,longname::String,units::String)::BoundaryCondition where I <: Integer

    i,j,k,wet = boundaryconditionatts(dim,dimval,γ)

    tracer = Array{Float64}(undef,size(wet))
    tracer[wet] .= zero(Float64)
    tracer[.!wet] .= zero(Float64)/zero(Float64)
    b = BoundaryCondition(tracer,i,j,k,dim,dimval,wet,name,longname,units)

end

"""
    function ones(dim::Int64,dimval::Int64,γ::Grid)::BoundaryCondition

       Initialize boundary condition with ones
"""
function ones(dim::I,dimval::I,γ::Grid,name::Symbol,longname::String,units::String)::BoundaryCondition where I <: Integer

    i,j,k,wet = boundaryconditionatts(dim,dimval,γ)

    tracer = Array{Float64}(undef,size(wet))
    tracer[wet] .= ones(Float64)
    tracer[.!wet] .= zero(Float64)/zero(Float64)
    b = BoundaryCondition(tracer,i,j,k,dim,dimval,wet)

end

"""
   Get boundary condition by extracting from 3D tracer
"""
function getboundarycondition(tracer3d::Field,dim::Integer,dimval::Integer,γ::Grid)::BoundaryCondition

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


"""
    function getboundarycondition(field::Field,dim,dimval)::BoundaryCondition
   Get boundary condition by extracting from Field (i.e., 3D tracer)
# Arguments
- `field::Field`: 3D tracer field with metadata and grid
- `dim`: dimension number (1,2,3) that the boundary plane has constant value
- `dimval`: index number in dimension `dim` that defines boundary plane
# Output
- `b::BoundaryCondition`: boundary condition on a plane with metadata and grid
"""
function getboundarycondition(field::Field,dim,dimval)::BoundaryCondition

    dimsize = size(field.γ.wet)
    # dumb way to do it
    if dim == 1
        wet2d = field.γ.wet[dimval,:,:]
        tracer2d = field.tracer[dimval,:,:]
        i = field.γ.lat
        j = field.γ.depth
        k = field.γ.lon[dimval]
    elseif dim == 2
        wet2d = field.γ.wet[:,dimval,:]
        tracer2d = field.tracer[:,dimval,:]
        i = field.γ.lon
        j = field.γ.depth
        k = field.γ.lat[dimval]
    elseif dim == 3
        wet2d = field.γ.wet[:,:,dimval]
        tracer2d = field.tracer[:,:,dimval]
        i = field.γ.lon
        j = field.γ.lat
        k = field.γ.depth[dimval]
    else
        error("boundary condition not implemented in 4+ dimensions")
    end
    
    b = BoundaryCondition(tracer2d,i,j,k,dim,dimval,wet2d,
                          field.name,field.longname,field.units)

end

vec(u::BoundaryCondition) = u.tracer[u.wet]
