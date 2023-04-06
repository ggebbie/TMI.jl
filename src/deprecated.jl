
# """
#     function adjustboundarycondition(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where N1, N2, T <: Real

#     adjust all boundary conditions b that are described in u
# """
# function adjustboundarycondition(b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}},u::NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}) where {N1, N2, T <: Real}

#     bnew = deepcopy(b)
#     ukeys = keys(u)
#     for ukey in keys(u)
#         bnew[ukey].tracer[bnew[ukey].wet] += u[ukey].tracer[bnew[ukey].wet] 
#     end
#     return bnew
# end

# function unvec!(u::BoundaryCondition{T},uvec::Vector{T}) where T <: Real
#     I = findall(u.wet)
#     counter = 0
#     for i in I
#         counter +=1
#     # doesn't work to pass u back
#         #u.tracer[u.wet] = uvec
#         u.tracer[i] = uvec[counter]
#     end
# end
# function unvec!(u::NamedTuple,uvec::Vector) #where {N, T <: Real}

#     #counter = 0
#     for v in u
#         # if v isa BoundaryCondition
#         #     wet = v.wet
#         # elseif v isa Field
#         #     wet = v.γ.wet
#         # end
#         n = sum(wet(v))
#         v.tracer[wet(v)] = uvec[counter+1:counter+n]
#         counter += n
#     end
# end


# function zerotemplate(utemplate,uvec)
#     tracer = zeros(wet(utemplate))
#     tracer[wet(utemplate)] = uvec
#     return tracer
# end

# """
#     function unvec(u,uvec)

#     Replace u with new u
#     Undo the operations by vec(u)
#     Needs to update u because attributes of 
#     u need to be known at runtime.
# """
# function unvec(utemplate::Union{Field{T},BoundaryCondition{T}},uvec::Vector{T}) where T <: Real
#     tracer = zerotemplate(utemplate,uvec)
#     u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,wet(utemplate))
#     return u
# end
# function unvec(utemplate::BoundaryCondition{T},uvec::Vector{T}) where T <: Real
#     tracer = zeros(utemplate.wet)
#     tracer[utemplate.wet] = uvec
#     u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,utemplate.wet)
#     return u
# end
# function unvec(utemplate::Field{T},uvec::Vector{T}) where T <: Real
#     tracer = zeros(utemplate.γ)
#     tracer.tracer[utemplate.γ.wet] .= uvec
#     #u = BoundaryCondition(tracer,utemplate.i,utemplate.j,utemplate.k,utemplate.dim,utemplate.dimval,utemplate.wet)
#     return tracer
# end
# function unvec(utemplate::NamedTuple{<:Any,NTuple{N,BoundaryCondition{T}}},uvec::Vector{T}) where {N, T <: Real}
#     counter = 0
#     vals = Vector{BoundaryCondition}(undef,N)
#     for (ii,vv) in enumerate(utemplate)
#         n = sum(vv.wet)
#         vals[ii] = unvec(vv, uvec[counter+1:counter+n])
#         counter += n
#     end
#     u = (;zip(keys(utemplate), vals)...)
#     return u
# end
# function unvec(utemplate::NamedTuple{<:Any,NTuple{N,Field{T}}},uvec::Vector{T}) where {N, T <: Real}
#     counter = 0
#     vals = Vector{Field}(undef,N)
#     for (ii,vv) in enumerate(utemplate)
#         n = sum(vv.γ.wet)
#         vals[ii] = unvec(vv, uvec[counter+1:counter+n])
#         counter += n
#     end
#     u = (;zip(keys(utemplate), vals)...)
#     return u
# end
# function unvec(utemplate,uvec::Vector{T}) where {T <: Real}
#     counter = 0
#     vals = Vector{Any}(undef,length(utemplate))
#     for (ii,vv) in enumerate(utemplate)
#         if vv isa Field
#             n = sum(vv.γ.wet)
#         elseif vv isa BoundaryCondition
#             n = sum(vv.wet)
#         end
#         vals[ii] = unvec(vv, uvec[counter+1:counter+n])
#         counter += n
#     end
#     u = (;zip(keys(utemplate), vals)...)
#     return u
# end


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

# """
#     `function +(c::BoundaryCondition,d::BoundaryCondition)::BoundaryCondition`
#     Define addition for Fields
# """
# function Base.:+(c::T,d::T)::T where T <: Union{Source,Field,BoundaryCondition}

#     if c.wet != d.wet # check conformability
#         error("BoundaryCondition's not conformable for addition")
#     end
#     array = zeros(c.wet)
#     e = BoundaryCondition(array,c.i,c.j,c.k,c.dim,c.dimval,c.wet)
    
#     # a strange formulation to get
#     # return e to be correct
#     e.tracer[e.wet] += c.tracer[c.wet]
#     e.tracer[e.wet] += d.tracer[d.wet]
#     return e
# end

# """
#     `function +(c::Field,d::Field)::Field`
#     Define addition for Fields
# """
# function +(c::Field{T},d::Field{T})::Field{T} where T <: Real
#     # initialize output
#     if c.γ.wet != d.γ.wet # check conformability
#         error("Fields not conformable for addition")
#     end

#     if !isequal(d.units,c.units)
#         error("Units not consistent:",d.units," vs ",c.units)
#     end

#     e = zeros(d.γ,d.name,d.longname,d.units)

#     # a strange formulation to get
#     # return e to be correct
#     e.tracer[e.γ.wet] += c.tracer[c.γ.wet]
#     e.tracer[e.γ.wet] += d.tracer[d.γ.wet]
#     return e
# end

# """
#     `function -(c::Field,d::Field)::Field`
#     Define addition for Fields
# """
# function -(c::Field{T},d::Field{T})::Field{T} where T <: Real
#     # initialize output
#     if c.γ.wet !== d.γ.wet # check conformability
#         error("Fields not conformable for addition")
#     end
#     e = zeros(d.γ)
#     e.tracer[e.γ.wet] += c.tracer[c.γ.wet]
#     e.tracer[e.γ.wet] -= d.tracer[d.γ.wet]
#     return e
# end
