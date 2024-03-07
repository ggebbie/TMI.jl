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
struct MassFraction{T <: Real}
    fraction::Array{T,3}
    γ::Grid
    name::Symbol
    longname::String
    units::String
    position::CartesianIndex{3}
end


function MassFraction(A,
    γ::Grid,
    Δ::CartesianIndex;
    wrap=(true, false, false),
    longname = "mass fraction from neighbor")

    # allocate masks
    ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    wet = falses(ngrid)
    m   = NaN*ones(ngrid)
    R   = copy(γ.R)
    Iint = cartesianindex(γ.interior)
    for I in Iint
        #Istep = I + step
        Istep, inbounds = step_cartesian(I::CartesianIndex, Δ::CartesianIndex, γ::Grid; wrap=wrap)

        if inbounds
            if γ.wet[Istep]
                wet[I] = true
                m[I] = A[R[I],R[Istep]]
            end
        end
    end

    # make a field
    return MassFraction(m,
        Grid(γ.lon,γ.lat,γ.depth,
            wet,
            γ.interior),
        :m_north,
        longname,
        "unitless",
        Δ)
end

"""
function `step_cartesian(I, Δ; wrap=(true,false,false))`

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
                wrapstep[idim] = -ngrid[idim] # wr
            end

            # check lower bound
            if Ilo[idim]>0 && wrap[idim]
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
    wrap=(true, false, false), longname = "mass fraction from northern neighbor")

massfractions_east(A, γ) = MassFraction(A, γ, CartesianIndex(1,0,0);
    wrap=(true, false, false), longname = "mass fraction from eastern neighbor")

massfractions_south(A, γ) = MassFraction(A, γ, CartesianIndex(0,-1,0);
    wrap=(true, false, false), longname = "mass fraction from southern neighbor")
    
massfractions_west(A, γ) = MassFraction(A, γ, CartesianIndex(-1,0,0);
    wrap=(true, false, false), longname = "mass fraction from western neighbor")

massfractions_up(A, γ) = MassFraction(A, γ, CartesianIndex(0,0,-1);
    wrap=(true, false, false), longname = "mass fraction from upper neighbor")
    
massfractions_down(A, γ) = MassFraction(A, γ, CartesianIndex(0,0,1);
    wrap=(true, false, false), longname = "mass fraction from lower neighbor")

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
function neighbors(m::NamedTuple,γ::Grid{T}) where T

    ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    if T == Float32
        n = zeros(Int32,ngrid)
    elseif T == Float64
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
    Δ=neighbor_indices(6),
    wrap=(true, false, false),
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
function neighbors(γ::Grid{T};
    Δ=neighbor_indices(6),
    wrap=(true, false, false),
    longname = "number of neighbors") where T

    # allocate masks
    ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    if T == Float32
        n = zeros(Int32,ngrid)
    elseif T == Float64
        n = zeros(Int64,ngrid)
    else
        error("TMI.neighbors, type "*T*" not implemented")
    end
    
    Iint = cartesianindex(γ.interior)
    for d in Δ
        for I in Iint
            Istep, inbounds = step_cartesian(I::CartesianIndex,
                d::CartesianIndex,
                γ::Grid;
                wrap=wrap)

            if inbounds
                n[I] += γ.wet[Istep]
            end
        end
    end

    return Field(n,γ,:n,longname,"unitless")
end


function massfractions_isotropic(γ::Grid)

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
function massfractions(c::NamedTuple, w::NamedTuple; alg = :local) # assume: `Field`s inside
    if alg == :local
        return local_solve(c::NamedTuple, w::NamedTuple)
    else
        println("TMI.massfractions: global solution for TMI water-mass matrix not implemented yet")
    end
end

"""
function local_solve(c::NamedTuple, w::NamedTuple)

Given tracers, solve for mass fractions using
a repeated local algorithm.

`local_solve!` is the workhorse for this algorithm and
specifies a default inverse tapering parameter
of `α=100_000`.

# Arguments
- `c::NamedTuple`: `Field` tracers (observations or model output)
- `w::NamedTuple`: scale size of tracers (for weighting observational fits)
# Output
- `m̃::NamedTuple`: `MassFraction` in a NamedTuple collection
"""
function local_solve(c::NamedTuple, w::NamedTuple) # assume: `Field`s inside
    γ = first(c).γ # assumption: all grids match up
    m̃ = TMI.massfractions_isotropic(γ) # good first guess
    #TMI.local_solve!(m̃,c,w)
    TMI.local_solve_old!(m̃, c, w)
    return m̃
end
function local_solve!(m::NamedTuple,c::NamedTuple, w::NamedTuple; alg = :quadprog)

    γ = first(c).γ
    wlocal = [w1 for w1 in w] # unpack Named Tuple
    Iint = cartesianindex(γ.interior)
    n   = TMI.neighbors(m,γ)
    nrow   = length(c) + 1 # add mass conservation

    # allocate maximum needed
    nmax = maximum(n)
    mlocal = zeros(nmax)
    nlocal = zeros(nrow)
    ϵ = 1e-8 # for checking tolerances
    b = vcat(zeros(nrow-1),1.0)
    
    for I in Iint

        # well-mixed first guess
        ncol = n.tracer[I]
        m0local = ones(ncol) ./ ncol
        
        # does it already fit the data?
        #Alocal = local_watermass_matrix(c,m,I,ncol.tracer[I])
        Alocal, single_connection = local_watermass_matrix(c,m,I,n)

        n0 = Alocal*m0local # data misfit (mass excluded)

        if sum(abs.(n0)) < ϵ # something small
            mlocal[1:ncol] = m0local
        else
    
            # Invert! by maximizing mixing and fitting tracers/mass perfectly
            # attempts to fit tracers and mass conservation perfectly
            Alocal2 = vcat(Alocal,ones(1,size(Alocal,2)))
            mlocal[1:ncol] = m0local + Alocal2\vcat(n0,0.0)
            nlocal[1:nrow] = vcat(Alocal*mlocal[1:ncol],
                1-sum(mlocal[1:ncol]))

            # check fit and check non-negativity
            if single_connection ||
                ((sum(abs.(nlocal)) > ϵ) ||
                    !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ))

                model, x = local_quadprog(Alocal2,m0local)

                if termination_status(model) != OPTIMAL
                    println("2nd try not feasible ",I)
                    # relax tracer assumption
                    model2, x2 = local_quadprog(Alocal,
                        m0local[1:ncol],wlocal)
                    if termination_status(model2) != OPTIMAL
                        println("WARNING: NEVER FOUND A BETTER SOLUTION!")
                        mlocal[1:ncol] = m0local
                    else
                        mlocal[1:ncol] = value.(x2)
                    end
                else
                    mlocal[1:ncol] = value.(x)
                end
            end
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

function local_quadprog(m0local::Vector, Alocal::Matrix, b::Vector)
    #model =
    #     Model(optimizer_with_attributes(COSMO.Optimizer, "verbose" => false,
    #         "eps_abs" => 1e-4, "eps_rel" = 1e-4))

    model = Model(COSMO.Optimizer);
    #model = Model(HiGHS.Optimizer);
    set_time_limit_sec(model, 0.1)
    set_silent(model);
    ncol = length(m0local)
    @variable(model, 0.0 <= x[i = 1:ncol] <= 1.0 ) 
    @constraint(model, Alocal*x == b)
    #@constraint(model, sum(x) == 1.0)
    @objective(model, Min,
        local_objective_mixing(x,m0local))
    optimize!(model)
    return model, x 
end

local_objective_mixing(x,m0local) =
    sum((x.-m0local).^2)
    
function local_quadprog(Alocal::Matrix,m0local::Vector,w::Vector)
    #settings = COSMO.Settings(verbose = false, eps_abs = 1e-2, eps_rel = 1e-2)
    model =
         Model(optimizer_with_attributes(COSMO.Optimizer, "verbose" => false,
             "eps_abs" => 1e-4, "eps_rel" => 1e-4))
    set_time_limit_sec(model, 0.1)

    #model = Model(COSMO.Optimizer);
    #set_silent(model);
    ncol = length(m0local) 
    @variable(model, 0.0 <= x[i = 1:ncol] <= 1.0 ) 
    @constraint(model, sum(x) == 1.0)
    #@objective(model, Min,
    #    sum(((Alocal[1:end-1,:]*x)./w).^2))
    #    println(local_objective_obs(Alocal,x,w))

    # cut off mass conservation from this objective function
    @objective(model, Min,
        local_objective_obs(Alocal[1:end-1,:],x,w[1:end-1]))
    optimize!(model)
    return model, x
end

local_objective_obs(Alocal,x,w) =
    sum(((Alocal*x)./w).^2)

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

function local_solve_old!(m::NamedTuple, c::NamedTuple, w::NamedTuple ; alg = :quadprog)

    γ = first(c).γ
    #Rfull = CartesianIndices(γ.wet)
    Iint = cartesianindex(γ.interior)

    n   = TMI.neighbors(m,γ)
    nrow   = length(c) + 1 # add mass conservation

    # allocate maximum needed
    nmax = maximum(n)
    b = vcat(zeros(nrow-1),1.0)
    mlocal = zeros(nmax)
    nlocal = zeros(nrow)
    ϵ = 1e-6 # tolerance of mass conservation
    Jlimit = 1e-4 # scaled cost function target
    wlocal = vcat([w1 for w1 in w],ϵ)

    for I in Iint

        # well-mixed first guess
        ncol = n.tracer[I]
        m0local = ones(ncol) ./ ncol
        Alocal, single_connection = local_watermass_matrix_old(c,m,I,n)
        n0 = b - Alocal*m0local
        J0 = sum((n0./wlocal).^2)/nrow

        if J0 < Jlimit # hit target already?
            mlocal[1:ncol] = m0local
            #println("first guess good enough @ ",I)
        else
            mlocal[1:ncol] = m0local + Alocal\n0
            nlocal[1:nrow] = b - Alocal*mlocal[1:ncol]
            Jlocal = sum((nlocal./wlocal).^2)/nrow

            # check connection, fit, and non-negativity
            if single_connection ||
                ((Jlocal > Jlimit) ||
                    !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ))

                # quadratic programming
                # println("run local quadprog @ ",I)
                model, x = local_quadprog(m0local, Alocal, b)

                if termination_status(model) != OPTIMAL
                    model2, x2 = local_quadprog(Alocal,
                        m0local, wlocal)
                    if termination_status(model2) != OPTIMAL
                        println("TMI.local solve: WARNING: NEVER FOUND A BETTER SOLUTION! @ ",I)
                        mlocal[1:ncol] = m0local
                    else
                        # println("satisfied only mass cons @ ",I)
                        mlocal[1:ncol] = value.(x2)
                    end
                else
                    #println("solved for both perfect data and mass cons @ ",I)
                    mlocal[1:ncol] = value.(x)
                end
            end
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
    nrow   = length(c) # + 1 # add mass conservation

    single_connection = false # warning of singularity if one of the neighbors is singly connected
    
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
function local_watermass_matrix_old(c::NamedTuple,
    m::NamedTuple,
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

Base.vec(m::MassFraction) = m.fraction[m.γ.wet]
Base.length(m::MassFraction) = sum(m.γ.wet)
Base.maximum(m::MassFraction) = maximum(m.fraction[m.γ.wet])
Base.minimum(m::MassFraction) = minimum(m.fraction[m.γ.wet])
