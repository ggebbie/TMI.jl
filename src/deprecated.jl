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
    function Γsfc 
    Γsfc anonymously calls surfacecontrol2field
"""
#Γsfc = surfacecontrol2field

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
    [sumwis[i] = Interpolations.InterpGetindex(wetwrap)[wis[i]...] for i in eachindex(wis)]

    # reconstruct the observations
    ỹ = Vector{Float64}(undef,N)
    c̃ = tracerinit(γ.wet)
    c̃[γ.wet] = cvec
    replace!(c̃,NaN=>0.0)
    cwrap = view(c̃,list,:,:)

    # divide by sum of all ocean weights so that this is still a true average
    [ỹ[i] = Interpolations.InterpGetindex(cwrap)[wis[i]...]/sumwis[i] for i in eachindex(wis)]
    return ỹ
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
Read 3D fields from mat file and save to NetCDF file.
"""
function matfields2nc_orig(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    filemat = pkgdatadir("TMI_"*TMIversion*".mat")
    vars = matread(filemat)

    TMIgrids, TMIgridsatts = griddicts(γ)

    T = eltype(γ.lon) # does the eltype of longitude have to equal the tracer eltype?
    #T =  Float64

    # iterate over all possible variables listed above
    Izyx = cartesianindex(filemat)
    TMIfields = Dict{String,Array{T,3}}()
    for (kk,vv) in mat2ncfield()
        haskey(vars,kk) ? push!(TMIfields, vv => tracerinit(vars[kk], Izyx, γ.wet)) : nothing
    end

    # also save fields that are stored in the x struct, if they exist
    if haskey(vars,"x")
        for (kk,vv) in mat2ncfield()
            println(kk)
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
     function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    Read an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)
Return a BoundaryCondition

    Version 1: operates on a 2D Float field 
"""
function surfaceregion(TMIversion::String,region::Union{String,Symbol},γ::Grid)::BoundaryCondition

    file = pkgdatadir("regions_"*TMIversion*".nc") 
    
    # version 1: region masks were stored with "d" prefix
    # new version: regions defined by region name alone
    tracername = "d_"*String(region) 

    #dsfc = ncread(file,tracername) # using NetCDF.jl
    ds = Dataset(file,"r") # using NCDatasets.jl
    v = ds[tracername]
    units = v.attrib["units"]
    longname = v.attrib["longname"]

    lon = γ.lon
    lat = γ.lat
    depth = γ.depth[1]
    mask = v[:,:] # Float

    #name = v.attrib["name"] # error
    close(ds)
    
    b = BoundaryCondition(mask,lon,lat,depth,3,1,γ.wet[:,:,1],Symbol(region),longname,units)
    return b
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

function boundaryconditionatts_old(dim::Integer,dimval::Integer,γ::Grid)

    # dimsize = size(γ.wet)
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

function getboundarycondition_old(tracer3d::Field,dim::Integer,dimval::Integer,γ::Grid)::BoundaryCondition

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

function getboundarycondition_old(field::Field,dim,dimval)::BoundaryCondition

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
