"""
    struct MassFractions

If 3D, then names known ahead of time and
can use `struct` instead of NamedTuple

# Attributes
    `north::Field{T}`
    `east::Field{T}`
    `south::Field{T}`
    `west::Field{T}`
    `up::Field{T}`
    `down::Field{T}`
"""
struct MassFractions{T <: Real}
    fraction::Array{T,3}
    γ::Grid
    name::Symbol
    longname::String
    units::String
    position::CartesianIndex{3}
end

function MassFractions(A,
    γ::Grid,
    Δ::CartesianIndex;
    wrap=(true, false, false),
    longname = "mass fraction from neighbor")

    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    # allocate masks
    ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    wet = falses(ngrid)
    m   = NaN*ones(ngrid)
    
    for I in Iint
        #Istep = I + step
        Istep, inbounds = step_cartesian(I::CartesianIndex, Δ::CartesianIndex, γ::Grid; wrap=(true,false,false))

        if inbounds
            if γ.wet[Istep]
                wet[I] = true
                m[I] = A[R[I],R[Istep]]
            end
        end
    end

    # make a field
    return MassFractions(m,
        Grid(γ.lon,γ.lat,γ.depth,
            wet,
            γ.interior),
        :m_north,
        longname,
        "unitless",
        Δ)
end

"""
function step_cartesian(I, Δ; wrap=(true,false,false))

# Arguments
- `I::CartesianIndex`: starting point
- `Δ::CartesianIndex`: step
- `γ::Grid`: TMI-defined grid
- `wrap::Tuple{Bool}`: does the dimension wrap around?

# Output
- `Istep::CartesianIndex`: new location
- `inbounds::Bool`: inside the domain bounds?
"""
function step_cartesian(I::CartesianIndex,
    Δ::CartesianIndex,
    γ::Grid;
    wrap=(true,false,false))

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
        ngrid = (length(γ.lon),
            length(γ.lat),
            length(γ.depth))

        wrapstep = zeros(Int,length(wrap))
        for idim in eachindex(wrap)
            # check upper bound
            if Ihi[idim]>0 && wrap[idim]
                #Istep[idim] .-= ngrid[idim] # wr
                wrapstep[idim] = -ngrid[idim] # wr
            end

            # check lower bound
            if Ilo[idim]>0 && wrap[idim]
                #Istep[idim] .+= ngrid[idim] # wr
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

massfractions_north(A, γ) = MassFractions(A, γ, CartesianIndex(0,1,0);
    wrap=(true, false, false), longname = "mass fraction from northern neighbor")

massfractions_east(A, γ) = MassFractions(A, γ, CartesianIndex(1,0,0);
    wrap=(true, false, false), longname = "mass fraction from eastern neighbor")

massfractions_south(A, γ) = MassFractions(A, γ, CartesianIndex(0,-1,0);
    wrap=(true, false, false), longname = "mass fraction from southern neighbor")
    
massfractions_west(A, γ) = MassFractions(A, γ, CartesianIndex(-1,0,0);
    wrap=(true, false, false), longname = "mass fraction from western neighbor")

massfractions_up(A, γ) = MassFractions(A, γ, CartesianIndex(0,0,-1);
    wrap=(true, false, false), longname = "mass fraction from upper neighbor")
    
massfractions_down(A, γ) = MassFractions(A, γ, CartesianIndex(0,0,1);
    wrap=(true, false, false), longname = "mass fraction from lower neighbor")

function tracer_contribution(c::Field,m::Union{MassFractions,NamedTuple})

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

function tracer_contribution!(mc::Field,c::Field,m::MassFractions)

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
function neighbors(m::NamedTuple,γ::Grid)

    ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    n = zeros(Int64,ngrid)
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

function massfractions_isotropic(γ::Grid)

    # get a sample with zeros
    nfield = sum(γ.wet)
    A = spzeros(nfield,nfield)

    m = (north = TMI.massfractions_north(A,γ),
        east   = TMI.massfractions_east(A,γ),
        south  = TMI.massfractions_south(A,γ),
        west   = TMI.massfractions_west(A,γ),
        up     = TMI.massfractions_up(A,γ),
        down   = TMI.massfractions_down(A,γ))

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
function massfractions(c::NamedTuple; alg = :local)

Create NamedTuple of mass fractions from observations `c`

Doesn't produce a `MassFractions` struct and thus
is named in lower case.

# Arguments
- `c::NamedTuple`: input observations
- `alg=:local`: default algorithm is `:local`
# Output
- `m::NamedTuple`: collection of mass fractions
"""
function massfractions(c::NamedTuple; alg = :local) # assume: `Field`s inside
    if alg == :local
        return local_solve(c::NamedTuple)
    else
        println("TMI.massfractions: global solution for TMI water-mass matrix not implemented yet")
    end
end

"""
function local_solve(c::NamedTuple)

Given tracers, solve for mass fractions using
a repeated local algorithm.

`local_solve!` is the workhorse for this algorithm and
specifies a default inverse tapering parameter
of `α=100_000`.

# Arguments
- `c::NamedTuple`: `Field` tracers (observations or model output)
# Output
- `m̃::NamedTuple`: `MassFractions` in a NamedTuple collection
"""
function local_solve(c::NamedTuple) # assume: `Field`s inside
    γ = first(c).γ # assumption: all grids match up
    m̃ = TMI.massfractions_isotropic(γ) # good first guess
    TMI.local_solve!(m̃,c)
    return m̃
end
function local_solve!(m::NamedTuple,c::NamedTuple; alg = :quadprog)

    γ = first(c).γ
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    neighbors   = TMI.neighbors(m,γ)
    nrow   = length(c) + 1 # add mass conservation

    # allocate maximum needed
    nmax = maximum(neighbors)
    l = zeros(nmax)
    u = ones(nmax)
    b = vcat(zeros(nrow-1),1.0)
    mlocal = zeros(nmax)
    nlocal = zeros(nrow)
    ϵ = 1e-8 # for checking tolerances
    
    for I in Iint

        # well-mixed first guess
        ncol = neighbors.tracer[I]
        m0 = ones(ncol) ./ ncol

        # does it already fit the data?
        #Alocal = local_watermass_matrix(c,m,I,ncol.tracer[I])
        Alocal, single_connection = local_watermass_matrix(c,m,I,neighbors)

        n0 = b - Alocal*m0

        if sum(abs.(n0[1:ncol])) < 1.0e-8 # something small
            mlocal[1:ncol] = m0
        else
    
            # Invert! by maximizing mixing and fitting tracers/mass perfectly

            # attempts to fit tracers and mass conservation perfectly
            #m_local[1:nlocal] = x0 + Alocal'*((Alocal*Alocal')\noise)
            mlocal[1:ncol] = m0 + Alocal\n0
            nlocal[1:nrow] = b - Alocal*mlocal[1:ncol]

            # check fit and check non-negativity
            if (sum(abs.(nlocal)) > 1.0e-8) || !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ) 
                # In less perfect cases, consider forcing m_f = sum_i=1^(f-1) m_i
                # With fewer tracers, incorporate a first guess with horizontal (perhaps) motion
                # May require non-negative least-squares or quadratic programmin (JuMP)

                # quadratic programming
                model = Model(HiGHS.Optimizer);
                set_silent(model);
                #nn = size(Alocal,2)
                @variable(model, l[i] <= x[i = 1:ncol] <= u[i] ) 
                #@constraint(model, sum(x) == 1.0)
                @constraint(model, Alocal*x == b)
                @objective(model, Min, sum((x.-m0).^2))
                optimize!(model)
                if termination_status(model) != OPTIMAL
                    println("not optimal")
                end
                mlocal[1:ncol] = value.(x)
            end
        
            # if !single_connection 
            #     m_local = Alocal\b
            #     if sum(abs.(m_local)) > 1.0000001
            #         println("loc",I)
            #         println(sum(abs.(m_local)))
            #     end
            # elseif single_connection # problem: some cells have only one neighbor
            #     #elseif ncol.tracer[I] == 1
            #     println(single_connection)
            #     m_local = (Alocal'*Alocal+LinearAlgebra.I./α)\(Alocal'*b + x0./α)
            
            # else
            #     println("no neighbors at all?")
            # end

        end
        
        # Save into proper mass fractions
        i = 0
        for m1 in m
            if m1.γ.wet[I]
                i += 1
                # should not need to check bounds 
                m1.fraction[I] = mlocal[i]
            end
        end
    end
end

function local_watermass_matrix(c::NamedTuple,
    m::NamedTuple,
    I::CartesianIndex,
    nlocal::Real)

    #γ = first(c).γ
    ncol   = nlocal # number of neighbors
    nrow   = length(c) + 1 # add mass conservation

    # allocate tracer matrix
    A = ones(nrow,ncol) # then overwrite later
    
    # loop through all mass fractions
    #for i1 in eachindex(m)
    i = 0; j = 0
    for m1 in m

        if m1.γ.wet[I]
            j += 1
            i = 0
            # grab the tracer values to fill a column
            # is this the correct γ?
            Istep, _ = step_cartesian(I,
                m1.position,
                m1.γ,
                wrap=(true,false,false))

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
    m::NamedTuple,
    I::CartesianIndex,
    neighbors::Field)

    #γ = first(c).γ
    ncol   = neighbors.tracer[I] # number of neighbors
    nrow   = length(c) + 1 # add mass conservation

    single_connection = false # warning of singularity if one of the neighbors is singly connected
    
    # allocate tracer matrix
    A = ones(nrow,ncol) # then overwrite later
    
    # loop through all mass fractions
    #for i1 in eachindex(m)
    i = 0; j = 0
    for m1 in m

        if m1.γ.wet[I]
            j += 1
            i = 0
            # grab the tracer values to fill a column
            # is this the correct γ?
            Istep, _ = step_cartesian(I,
                m1.position,
                m1.γ,
                wrap=(true,false,false))

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

"""
function watermassmatrix(m::NamedTuple, γ::Grid;
    wrap=(true, false, false))

Produce water-mass matrix from mass fractions and grid.

# Arguments
- `m::NamedTuple`: collection of `MassFraction`s
- `γ::TMI.Grid`
- `wrap::(true, false, false)`: does the particular dimension have a wraparound (periodic) domain?

# Output
- `A`: sparse water-mass matrix
"""
function watermassmatrix(m::NamedTuple, γ::Grid;
    wrap=(true, false, false))

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
                Istep, _ = step_cartesian(I,
                    m1.position,
                    γ;
                    wrap=(true,false,false))

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

Base.vec(m::MassFractions) = m.fraction[m.γ.wet]
Base.length(m::MassFractions) = sum(m.γ.wet)
Base.maximum(m::MassFractions) = maximum(m.fraction[m.γ.wet])
Base.minimum(m::MassFractions) = minimum(m.fraction[m.γ.wet])
