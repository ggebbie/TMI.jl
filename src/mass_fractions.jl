"""
struct MassFraction

store mass fractions in a Field-like array

- `fraction::Array{T,3}`
- `γ::Grid`
- `name::Symbol`
- `longname::String`
- `units::String`
- `position::CartesianIndex{3}`
"""
struct MassFraction{T <: Real,R <: Real,N,F <: AbstractArray{T,N}}
    fraction::F
    γ::Grid{R,N}
    name::Symbol
    longname::String
    units::String
    position::CartesianIndex{N}
end

function MassFraction(A,
    γ::Grid,
    Δ::CartesianIndex;
    longname = "mass fraction from neighbor")

    # allocate masks
    #ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))

    axislabels = γ.axes
    ndims = length(axislabels)
    ngrid =
        Tuple([length(axislabels[d]) for d in 1:ndims])

    wet = falses(ngrid)
    m   = NaN*ones(ngrid)
    R   = copy(γ.R)
    Iint = cartesianindex(γ.interior)
    for I in Iint
        #Istep = I + step
        Istep, inbounds = step_cartesian(I::CartesianIndex, Δ::CartesianIndex, γ::Grid)

        if inbounds
            if γ.wet[Istep]
                wet[I] = true
                m[I] = A[R[I],R[Istep]]
            end
        end
    end

    # make a field
    return MassFraction(m,
        #Grid(γ.lon,γ.lat,γ.depth,
        Grid(γ.axes,
            wet,
            γ.interior,
            γ.wrap,
            γ.Δ),
        :m_ij,
        longname,
        "unitless",
        Δ)
end



Base.vec(m::MassFraction) = m.fraction[m.γ.wet]
Base.length(m::MassFraction) = sum(m.γ.wet)
Base.maximum(m::MassFraction) = maximum(m.fraction[m.γ.wet])
Base.minimum(m::MassFraction) = minimum(m.fraction[m.γ.wet])

wet(m::MassFraction) = m.γ.wet

function Base.vec(u::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    T = eltype(values(u)[1].fraction)
    #T = eltype(u)
    uvec = Vector{T}(undef,0)
    for v in u
        #append!(uvec,v.tracer[v.wet])
        append!(uvec,vec(v))
    end
    return uvec
end

function unvec!(u::MassFraction,uvec::Vector{T}) where T <: Real
    I = findall(wet(u)) # findall seems slow
    for (ii,vv) in enumerate(I)
        u.fraction[vv] = uvec[ii]
    end
end

function unvec!(u::NamedTuple{names, <:Tuple{Vararg{MassFraction}}},uvec::Vector) where names
    nlo = 1
    nhi = 0
    for v in u
        nhi += sum(wet(v))
        unvec!(v,uvec[nlo:nhi])
        nlo = nhi + 1
    end
end

function unvec(u₀::NamedTuple{names, <:Tuple{Vararg{MassFraction}}},uvec::Vector) where names
    u = deepcopy(u₀)
    unvec!(u,uvec)
    return u
end

function Base.similar(v::Vector{<:MassFraction})
    # return MassFraction[similar(mf) for mf in v]
    return map(similar, v)
end

function Base.similar(nt::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    return map(similar, nt)
end

function Base.similar(mf::MassFraction) where T
    similar_m_fraction = similar(mf.fraction)
    similar_m_fraction[.!(mf.γ.interior)] .= NaN

    return MassFraction(
        similar_m_fraction,
        mf.γ,
        mf.name,
        mf.longname,
        mf.units,
        mf.position
    )
end

"""
function massfractions(c::NamedTuple, w::NamedTuple; alg = :local)

Create NamedTuple of mass fractions from observations `c`

Doesn't produce a `MassFraction` struct and thus
is named in lower case.

# Arguments
- `c::NamedTuple`: input observations
- `w::NamedTuple`: scale size of observations (used if obs not fit exactly)
- `alg=:local`: default algorithm is `:local`
# Output
- `m::NamedTuple`: collection of mass fractions
"""
#function massfractions(c::NamedTuple, w::NamedTuple; alg = :local) # assume: `Field`s inside
function massfractions(c, w; alg = :local) # assume: `Field`s inside
    if alg == :local
        return local_solve(c::NamedTuple, w::NamedTuple)
    else
        println("TMI.massfractions: global solution for TMI water-mass matrix not implemented yet")
    end
end

local_solve() = nothing
local_quadprog() = nothing

"""
function `step_cartesian(I, Δ, γ)`

# Arguments
- `I::CartesianIndex`: starting point
- `Δ::CartesianIndex`: step
- `γ::Grid`: TMI-defined grid

# Output
- `Istep::CartesianIndex`: new location
- `inbounds::Bool`: inside the domain bounds?
"""
function step_cartesian(I::CartesianIndex{N},
    Δ::CartesianIndex{N},
    γ::Grid{R,N}) where {R,N}

    #R = copy(γ.R)
    Istep = I + Δ

    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # number of steps outside of bounds
    Ihi = max(Ilast,Istep)-Ilast
    Ilo = Ifirst - min(Ifirst,Istep)
    inbounds = iszero(Ihi) && iszero(Ilo)

    if iszero(Ihi) && iszero(Ilo)
        return Istep, true
    else
        # allocate masks
        ngrid =
            Tuple([length(γ.axes[d]) for d in 1:N])

        #ngrid = (length(γ.lon),
        #    length(γ.lat),
        #    length(γ.depth))

        wrapstep = zeros(Int,length(γ.wrap))
        for idim in eachindex(γ.wrap)
            # check upper bound
            if Ihi[idim]>0 && γ.wrap[idim]
                wrapstep[idim] = -ngrid[idim] # wr
            end

            # check lower bound
            if Ilo[idim]>0 && γ.wrap[idim]
                wrapstep[idim] = ngrid[idim] # wr
            end
        end

        # take another step to wrap around
        Istep += CartesianIndex(Tuple(wrapstep))
        
        Ihi = max(Ilast,Istep)-Ilast
        Ilo = Ifirst - min(Ifirst,Istep)
        inbounds = iszero(Ihi) && iszero(Ilo)
        return Istep, inbounds
    end
end

"""
function `neighbor_indices(n::Integer)`

Direction (step) of neighbors away from
a central point. Choose n = 6 (default) or n=26. 
# Argument
- `n=6`: max number of neighbors
# Output
- `In::Vector{CartesianIndex}`: indices of neighbors
"""
function neighbor_indices(n=6)
    if n == 6
        return [CartesianIndex(0,1,0),
            CartesianIndex(1,0,0),
            CartesianIndex(0,-1,0),
            CartesianIndex(-1,0,0),
            CartesianIndex(0,0,-1),
            CartesianIndex(0,0,1)]
    elseif n == 26
        return [CartesianIndex(i,j,k) for i=-1:1 for j = -1:1 for k=-1:1 if sum(abs(i)+abs(j)+abs(k)) > 0]
    end
end

massfractions_north(A, γ) = MassFraction(A, γ, CartesianIndex(0,1,0);
    longname = "mass fraction from northern neighbor")

massfractions_east(A, γ) = MassFraction(A, γ, CartesianIndex(1,0,0);
    longname = "mass fraction from eastern neighbor")

massfractions_south(A, γ) = MassFraction(A, γ, CartesianIndex(0,-1,0);
    longname = "mass fraction from southern neighbor")
    
massfractions_west(A, γ) = MassFraction(A, γ, CartesianIndex(-1,0,0);
    longname = "mass fraction from western neighbor")

massfractions_up(A, γ) = MassFraction(A, γ, CartesianIndex(0,0,-1);
    longname = "mass fraction from upper neighbor")
    
massfractions_down(A, γ) = MassFraction(A, γ, CartesianIndex(0,0,1);
    longname = "mass fraction from lower neighbor")

function tracer_contribution(c::Field,m::Union{MassFraction,NamedTuple})

    # allocate masks
    ngrid = (length(c.γ.lon), length(c.γ.lat), length(c.γ.depth))
    # tracer contribution: units C*kg
    mc = Field(zeros(ngrid),
        c.γ,
        :mc,
        "tracer contribution",
        c.units)

    tracer_contribution!(mc,c,m)
    return mc
end

function tracer_contribution!(mc::Field,c::Field,m::MassFraction)

    # loop over all locations with "m" values
    Im = cartesianindex(m.γ.wet)
    for I in Im
        #Istep = I + step
        Istep, _ = step_cartesian(
            I::CartesianIndex,
            m.position,
            c.γ)
        #; wrap=(true,false,false))

        mc.tracer[I] += m.fraction[I]*
            (c.tracer[Istep] - c.tracer[I])
    end
end
function tracer_contribution!(mc::Field,c::Field, m::NamedTuple)
    for m1 in m
        tracer_contribution!(mc,c,m1)
    end
end

"""
function `neighbors(m::NamedTuple,γ::Grid)`

How many neighbors does each grid cell have?
# Arguments
- `m::NamedTuple`: input mass fractions to obtain their stencil (opportunity to simplify)
- `γ::TMI.Grid`
# Output
- `n::Field`: integer number of neighbors
"""
function neighbors(m::Union{Vector,NamedTuple},
    γ::Grid{R,N}) where {R,N}

    ngrid =
        Tuple([length(γ.axes[d]) for d in 1:N])
    #ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    if R == Float32
        n = zeros(Int32,ngrid)
    elseif R == Float64
        n = zeros(Int64,ngrid)
    else
        error("TMI.neighbors, type "*T*" not implemented")
    end
    for m1 in m
        n += m1.γ.wet
    end
    
    # make a field
    return Field(n,
        γ,
        :n,
        "neighbors",
        "unitless")
end
"""
function `neighbors(γ::Grid;
    longname = "number of neighbors")`

How many neighbors does each grid cell have?
Conceptually, it only depends on the grid, but this
algorithm is slower than the one that takes mass
fractions as input.

# Arguments
- `γ::TMI.Grid`
# Output
- `n::Field`: integer number of neighbors
"""
function neighbors(γ::Grid{R,N};
    longname = "number of neighbors") where {R,N}

    # allocate masks
    #ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    ngrid =
        Tuple([length(γ.axes[d]) for d in 1:N])

    if R == Float32
        n = zeros(Int32,ngrid)
    elseif R == Float64
        n = zeros(Int64,ngrid)
    else
        error("TMI.neighbors, type "*T*" not implemented")
    end
    
    Iint = cartesianindex(γ.interior)
    for d in γ.Δ
        for I in Iint
            Istep, inbounds = step_cartesian(I,d,γ)
            if inbounds
                n[I] += γ.wet[Istep]
            end
        end
    end

    return Field(n,γ,:n,longname,"unitless")
end

function massfractions_isotropic(γ::Grid{R,3}) where R
    # only works for a 3D grid 
    # get a sample with zeros
    nfield = sum(γ.wet)
    A = spzeros(nfield,nfield)

    m = (north = massfractions_north(A,γ),
        east   = massfractions_east(A,γ),
        south  = massfractions_south(A,γ),
        west   = massfractions_west(A,γ),
        up     = massfractions_up(A,γ),
        down   = massfractions_down(A,γ))

    # get number of neighbors for each location
    n = neighbors(m,γ)

    for m1 in m
        Im = cartesianindex(m1.γ.wet)
        for I in Im
            m1.fraction[I] = 1.0/n.tracer[I]
        end
    end
    return m
end
function massfractions_isotropic(γ::Grid{R,1}) where R
    # only works for a 3D grid 
    # get a sample with zeros
    nfield = sum(γ.wet)
    A = spzeros(nfield,nfield)

    m = [TMI.MassFraction(A, γ, d,
        longname = "mass fraction")
        for d in γ.Δ]

    # get number of neighbors for each location
    n = neighbors(m,γ)

    for m1 in m
        Im = cartesianindex(m1.γ.wet)
        for I in Im
            m1.fraction[I] = 1.0/n.tracer[I]
        end
    end
    return m
end


"""
`function watermassmatrix(m::Union{NamedTuple,Vector}, γ::Grid)`

Produce water-mass matrix from mass fractions and grid.

# Arguments
- `m::NamedTuple`: collection of `MassFraction`s
- `γ::TMI.Grid`

# Output
- `A`: sparse water-mass matrix
"""
function watermassmatrix(m::Union{NamedTuple,Vector}, γ::Grid)

    # get total size of mass fractions

    #ksfc = 1 # surface could be at a different level
    #nsfc = sum(γ.wet[:,:,ksfc]) # surface points

    nfield = sum(γ.wet)
    nm = 0
    for m1 in m
        nm += length(m1)
    end

    ilist = Array{Int64}(undef,nm+nfield)
    jlist = Array{Int64}(undef,nm+nfield)
    mlist = Array{Float64}(undef,nm+nfield)

    # ones on diagonal
    for ii in 1:nfield
        ilist[ii] = ii
        jlist[ii] = ii
        mlist[ii] = 1.0
    end
    
    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    # allocate water mass matrix
    nfield = sum(γ.wet)
    # Note, flipped sign of interior points
    #A = spdiagm(nfield,nfield,ones(nfield))

    counter = nfield
    for I in Iint
        for m1 in m
        #for I in Iint
            nrow = R[I]
            if m1.γ.wet[I]
                Istep, _ = step_cartesian(I, m1.position, γ)
                counter += 1
                ilist[counter] = nrow
                jlist[counter] = R[Istep]
                mlist[counter] = -m1.fraction[I]
                #A[nrow,R[Istep]] = -m1.fraction[I]
            end
        end
    end
    return sparse(ilist,jlist,mlist)
end

function gwatermassmatrix(gA, m::Union{NamedTuple,Vector}, γ::Grid)

    nfield = sum(γ.wet)

    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    # allocate water mass matrix
    nfield = sum(γ.wet)
    # Note, flipped sign of interior points
    #A = spdiagm(nfield,nfield,ones(nfield))
    gm = similar(m)
    [gm1.fraction .= NaN for gm1 in gm]

    counter = nfield
    for I in Iint
        for gm1 in gm
        #for I in Iint
            nrow = R[I]
            if gm1.γ.wet[I]
                Istep, _ = step_cartesian(I, gm1.position, γ)
                counter += 1
                gm1.fraction[I] = -gA[nrow,R[Istep]]
            end
        end
    end

    return gm
end


function local_watermass_matrix(c::NamedTuple,
    m::NamedTuple,
    I::CartesianIndex,
    nlocal::Real)

    #γ = first(c).γ
    ncol   = nlocal # number of neighbors
    nrow   = length(c) #+ 1 # add mass conservation

    # allocate tracer matrix
    A = zeros(nrow,ncol) # then overwrite later
    
    # loop through all mass fractions
    #for i1 in eachindex(m)
    i = 0; j = 0
    for m1 in m

        if m1.γ.wet[I]
            j += 1
            i = 0
            # grab the tracer values to fill a column
            # is this the correct γ?
            Istep, _ = step_cartesian(I, m1.position, m1.γ)

            #if γ.wet[Istep]: should automagically be true
            for c1 in c
                i += 1
                A[i,j] = c1.tracer[Istep] - c1.tracer[I]
            end
        end
    end
    return A
end


"""
function `local_watermass_matrix(c::NamedTuple,
    m::NamedTuple,
    I::CartesianIndex,
    neighbors::Field)`
Find local water-mass matrix with singularity checker
(`true` if one neighbor only has a single connection
to the rest of the ocean)
# Arguments
- `c::NamedTuple`: input tracers
- `m::NamedTuple`: mass fractions for grid stencil
- `I::CartesianIndex`: local "location"
- `neighbors::Field`: integer number of neighbors
# Output
- `A::Matrix`: local water-mass matrix
- `single_connection::Bool`: true if flagged for singularity warning
"""
function local_watermass_matrix(c::NamedTuple,
    m::Union{NamedTuple,Vector},
    I::CartesianIndex,
    neighbors::Field)

    ncol   = neighbors.tracer[I] # number of neighbors
    nrow   = length(c)  + 1 # add mass conservation

    single_connection = false # warning of singularity if one of the neighbors is singly connected
    
    # allocate tracer matrix
    A = ones(nrow,ncol) # then overwrite later
    
    i = 0; j = 0
    for m1 in m

        if m1.γ.wet[I]
            j += 1
            i = 0
            # grab the tracer values to fill a column
            # is this the correct γ?
            Istep, _ = step_cartesian(I,m1.position,m1.γ)
            if neighbors.tracer[Istep] == 1
                single_connection = true
            end
            
            #if γ.wet[Istep]: should automagically be true
            for c1 in c
                i += 1
                A[i,j] = c1.tracer[Istep] - c1.tracer[I]
            end
        end
    end
    return A, single_connection
end
