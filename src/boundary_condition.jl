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
struct BoundaryCondition{T <: Real,
    R <: Real,
    N,
    G <: Integer,
    B <: AbstractArray{T,N}}
    tracer::B
    axes::NTuple{N,Vector{R}}
    k::R
    dim::G
    dimval::G
    wet::BitArray{N}
    name::Symbol
    longname::String
    units::String
end

# special show function if N = 2
function Base.show(io::IO,
    mime::MIME{Symbol("text/plain")},
    x::BoundaryCondition{T,R,2}) where {T,R}
    summary(io, x); println(io)
    print(io, "Field size ")
    println(io, size(x.tracer))
    show(io,mime,
        heatmap(transpose(x.tracer),
            zlabel=x.units,
            title=x.longname)
    )
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
function BoundaryCondition(tracer::B,
    axes::NTuple{N,Vector{R}},
    k::R,
    dim::G,
    dimval::G,
    wet::BitArray{N}) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}

    return BoundaryCondition(tracer,axes,k,dim,dimval,wet,:bc,"boundary condition","unknown") 
end

"""
    function write(file,b)

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
function boundaryconditionatts(dim::Integer,dimval::Integer,γ::Grid{R,N}) where {R,N}

    for n in 1:N 
        if n == dim
            k = γ.axes[n][dimval]
            ind = deleteat!(collect(1:N),n)
            # boundary axes
            baxes = γ.axes[ind]
            wet2d = copy(selectdim(γ.wet,dim,dimval))
            return baxes,k,wet2d
        end
    end
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
    baxes,k,wet = boundaryconditionatts(dim,dimval,γ)
    tracer = Array{Float64}(undef,size(wet))
    tracer[wet] .= zero(Float64)
    tracer[.!wet] .= zero(Float64)/zero(Float64)
    return BoundaryCondition(tracer,baxes,k,dim,dimval,wet,name,longname,units)
end

"""
    function ones(dim::Int64,dimval::Int64,γ::Grid)::BoundaryCondition

       Initialize boundary condition with ones
"""
function ones(dim::I,dimval::I,γ::Grid,name::Symbol,longname::String,units::String)::BoundaryCondition where I <: Integer
    baxes,k,wet = boundaryconditionatts(dim,dimval,γ)
    tracer = Array{Float64}(undef,size(wet))
    tracer[wet] .= ones(Float64)
    tracer[.!wet] .= zero(Float64)/zero(Float64)
    return BoundaryCondition(tracer,baxes,k,dim,dimval,wet,name,longname,units)
end

"""
   Get boundary condition by extracting from N-dimensional tracer and returning (N-1)-dimensional array
"""
function getboundarycondition(field::Field{T,R,N},dim::Integer,dimval::Integer,γ::Grid)::BoundaryCondition where {T<:Real,R<:Real,N}

    for n in 1:N 
        if n == dim
            k = γ.axes[n][dimval]
            ind = deleteat!(collect(1:N),n)
            # boundary axes
            baxes = γ.axes[ind]
            wet_nminus1 = copy(selectdim(γ.wet,dim,dimval))
            tracer_nminus1 = copy(selectdim(field.tracer,
                dim,
                dimval))

            return BoundaryCondition(tracer_nminus1,
                baxes,
                k,
                dim,
                dimval,
                wet_nminus1)
        end
    end
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
getboundarycondition(field::Field,dim::Integer,dimval::Integer) =
    getboundarycondition(field,dim,dimval,field.γ)
    
vec(u::BoundaryCondition) = u.tracer[u.wet]

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
function surfacepatch(lonbox,latbox,γ::Grid{R,3})::BoundaryCondition where R

    # input must be self consistent
    ((lonbox[1] > lonbox[2]) || (latbox[1] > latbox[2])) && error("surfacepatch: upper bound less than lower bound")

    # reference longitude to closest central longitude
    lon_central = (lonbox[1] + lonbox[2])/2

    # shift longitudes so that wraparound values are far away
    lon_shifted = replace(x -> (x > lon_central + 180) ? x - 360 : x, γ.lon)
    replace!(x -> (x < lon_central - 180) ? x + 360 : x, lon_shifted)
    
    # preallocate
    patch = zeros(γ.wet[:,:,1])

    # can you add a logical to a Float64? yes, it's 1.0
    [patch[i,j] +=  latbox[1] ≤ γ.lat[j] ≤ latbox[2] && lonbox[1] ≤ lon_shifted[i] ≤ lonbox[2] for i in eachindex(γ.lon) for j in eachindex(γ.lat)]

    # 3,1 to identify surface
    b = BoundaryCondition(patch,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:patch,"rectangular patch","none")
    
    return b
end

# define the correct dimension and index for each control plane
# maybe someday find a way to hide γ
zerosurfaceboundary(γ::Grid,name=:none,longname="unknown",units="unknown") = zeros(3,1,γ,name,longname,units)::BoundaryCondition

zeronorthboundary(γ,name=:none,longname="unknown",units="unknown") = zeros(2,maximum(latindex(γ.I)),γ,name,longname,units)::BoundaryCondition

zeroeastboundary(γ,name=:none,longname="unknown",units="unknown") = zeros(1,maximum(lonindex(γ.I)),γ,name,longname,units)::BoundaryCondition

zerosouthboundary(γ,name=:none,longname="unknown",units="unknown") = zeros(2,1,γ,name,longname,units)::BoundaryCondition

zerowestboundary(γ,name=:none,longname="unknown",units="unknown") = zeros(1,1,γ,name,longname,units)::BoundaryCondition

onesurfaceboundary(γ,name=:none,longname="unknown",units="unknown") = ones(3,1,γ,name,longname,units)::BoundaryCondition

onenorthboundary(γ,name=:none,longname="unknown",units="unknown") = ones(2,maximum(latindex(γ.I)),γ,name,longname,units)::BoundaryCondition

oneeastboundary(γ,name=:none,longname="unknown",units="unknown") = ones(1,maximum(lonindex(γ.I)),γ,name,longname,units)::BoundaryCondition

onesouthboundary(γ,name=:none,longname="unknown",units="unknown") = ones(2,1,γ,name,longname,units)::BoundaryCondition

onewestboundary(γ,name=:none,longname="unknown",units="unknown") = ones(1,1,γ,name,longname,units)::BoundaryCondition

getsurfaceboundary(c::Field) = getboundarycondition(c,3,1)::BoundaryCondition
getnorthboundary(c::Field) = getboundarycondition(c,2,maximum(latindex(c.γ.I)))::BoundaryCondition
geteastboundary(c::Field) = getboundarycondition(c,1,maximum(lonindex(c.γ.I)))::BoundaryCondition
getsouthboundary(c::Field) = getboundarycondition(c,2,1)::BoundaryCondition
getwestboundary(c::Field) = getboundarycondition(c,1,1)::BoundaryCondition

""" 
    function setboundarycondition!(d::Field,b::BoundaryCondition)
    apply boundary condition to the equation constraints
# Arguments
- `d`::Field, equation constraints (i.e., right hand side)
- `b`::BoundaryCondition
"""
function setboundarycondition!(d::Field,b::BoundaryCondition)
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
        gb = BoundaryCondition(gd.tracer[b.dimval,:,:],b.axes,b.k,b.dim,b.dimval,b.wet)
    elseif b.dim == 2
        #gb.tracer = gd.tracer[:,b.dimval,:]
        gb = BoundaryCondition(gd.tracer[:,b.dimval,:],b.axes,b.k,b.dim,b.dimval,b.wet)
    elseif b.dim == 3
        #gb.tracer .+= gd.tracer[:,:,b.dimval] 
        gb = BoundaryCondition(gd.tracer[:,:,b.dimval],b.axes,b.k,b.dim,b.dimval,b.wet)
    else
        error("controls not implemented for 4+ dimensions")
    end
    return gb
end

"""
    function setboundarycondition!(d::Field{T},b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}}) where {N, T <: Real}

    set all boundary conditions in a Named Tuple
"""
function setboundarycondition!(d::Field,b::NamedTuple)
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
function gsetboundarycondition(gd::Field{T},b::NamedTuple) where T <: Real

    gb1 = Vector{BoundaryCondition{T}}(undef,length(keys(b)))
    for (ii,vv) in enumerate(b)
        gb1[ii] = gsetboundarycondition(gd,vv)
    end

    # https://discourse.julialang.org/t/construct-namedtuple-dynamically/15394/7
    gb = (;zip(keys(b), gb1)...)
                                       
    return gb
end

"""
    function adjustboundarycondition!(b::Union{BoundaryCondition,NamedTuple},u::Union{BoundaryCondition,NamedTuple})

    adjust all boundary conditions b that are described in u
"""
adjustboundarycondition(b::BoundaryCondition,u::BoundaryCondition) = b + u
function adjustboundarycondition(b::Union{BoundaryCondition,NamedTuple},u::Union{BoundaryCondition,NamedTuple}) 
    bnew = deepcopy(b)
    adjustboundarycondition!(bnew,u)
    return bnew
end

"""
    function adjustboundarycondition!(b::Union{BoundaryCondition,NamedTuple},u::Union{BoundaryCondition,NamedTuple})

    adjust all boundary conditions b that are described in u

    warning: if u doesn't contain any boundary condition adjustments,
    nothing will change.
"""
function adjustboundarycondition!(b::BoundaryCondition,u::BoundaryCondition)
    # write it out so b changes when returned
    b.tracer[b.wet] += u.tracer[u.wet] 
end
function adjustboundarycondition!(b::NamedTuple,u::NamedTuple)
    # only the bkeys are certain to be type BoundaryCondition
    for bkey in keys(b)
        haskey(u,bkey) && adjustboundarycondition!(b[bkey],u[bkey]) 
    end
end

# Seems not to be general because gu overwritten.
function gadjustboundarycondition(gb::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real
    gu  = gb
    return gu
end

"""
    function gadjustboundarycondition!(b::BoundaryCondition{T},u::BoundaryCondition{T}) where T <: Real

    adjust the (one) boundary condition 
    Just copy the variable.
    Keep this function so that calling functions can look alike.
    Could probably combine with lower function, use Union type
"""
function gadjustboundarycondition!(gu::BoundaryCondition,gb::BoundaryCondition) 
    gu.tracer[gu.wet] += gb.tracer[gb.wet]
end
function gadjustboundarycondition!(gu::NamedTuple,gb::NamedTuple)
    for bkey in keys(gb)
        #gadjustboundarycondition!(gu,gb[bkey])
        haskey(gu,bkey) && gadjustboundarycondition!(gu[bkey],gb[bkey]) 
    end
end

#function gadjustboundarycondition(gb::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where {N1, N2, T <: Real}
function gadjustboundarycondition(gb::Union{NamedTuple,BoundaryCondition},u::NamedTuple) #where {N1, N2, T <: Real}
    gu = gb[keys(u)] # grab the parts of the named tuple corresponding to u
    return gu
end

function watermassmatrix(A, b::Union{BoundaryCondition,NamedTuple}, γ::Grid)
    bmask = boundarymask(b, γ)
    Abc = deepcopy(A)
    bmaskwet = bmask[γ.wet]
    for i in eachindex(bmaskwet)
        if bmaskwet[i] 
            Abc[i,:] .= 0.0
            Abc[i,i] = 1.0
        end
    end
    return Abc
end

"""
update water-mass matrix to agree with new boundaries in γ
"""
function watermassmatrix(A::SparseMatrixCSC, γ::Grid)
    bmask = (γ.wet .&&  .!γ.interior)
    Abc = deepcopy(A)
    bmaskwet = bmask[γ.wet]
    for i in eachindex(bmaskwet)
        if bmaskwet[i] 
            Abc[i,:] .= 0.0
            Abc[i,i] = 1.0
        end
    end
    return Abc
    
end

"""
     function boundarymask(b, γ::Grid)
"""
function boundarymask(b::BoundaryCondition, γ::Grid)
    boundary = falses(length.(γ.axes))
    if b.dim == 1
        boundary[b.dimval,:,:] += b.wet
    elseif b.dim == 2
        boundary[:,b.dimval,:] += b.wet 
    elseif b.dim == 3
        boundary[:,:,b.dimval] += b.wet
    else
        error("controls not implemented for 4+ dimensions")
    end
    return boundary 
end

function boundarymask(b::NamedTuple, γ::Grid)
    boundary = falses(length.(γ.axes))
    for b1 in b
        boundary = (boundary .|| boundarymask(b1, γ))
    end
    return boundary 
end

# update Grid with newly-defined boundaries
function Grid(b::Union{BoundaryCondition,NamedTuple}, γ::Grid)
    boundary = boundarymask(b, γ)

    return Grid(
        γ.axes,
        γ.wet,
        (γ.wet .&&  .!boundary),
        γ.wrap,
        γ.Δ)
end
