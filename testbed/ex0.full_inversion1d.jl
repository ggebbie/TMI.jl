import Pkg; Pkg.activate(".")

using Revise
using TMI
using SparseArrays

import TMI:MassFraction
import TMI:step_cartesian

# Add this after your MassFraction struct definition

# Copy constructor that preserves tracked types
function Base.similar(mf::MassFraction, ::Type{T}) where T
    return MassFraction(
        similar(mf.fraction, T, size(mf.fraction)),
        mf.γ,
        mf.name,
        mf.longname,
        mf.units,
        mf.position
    )
end

function Base.similar(mf::MassFraction) where T
    return MassFraction(
        similar(mf.fraction),
        mf.γ,
        mf.name,
        mf.longname,
        mf.units,
        mf.position
    )
end

function Base.similar(v::Vector{<:MassFraction})
    # return MassFraction[similar(mf) for mf in v]
    return map(similar, v)
end

# No type argument
function Base.similar(nt::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    return map(similar, nt)
end

ngrid = (500) # number of grid cells
xmax = 1000.0 # domain size 
lon = collect(range(0.0,1000.0,length=ngrid[1]))
tracer = collect(1.0.-lon./xmax)

axes = (lon,)
wet = trues(ngrid)
interior = copy(wet)
interior[begin] = false
interior[end] = false

wrap = (false,)
Δ = [CartesianIndex(1,),CartesianIndex(-1,)]
γ = Grid(axes,wet,interior,wrap,Δ)

m = massfractions_isotropic(γ)
[m1.fraction .*= 5*i for (i, m1) in enumerate(m)]
A = watermassmatrix(m, γ)

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





function Jsumquare(m, γ, m_template)
    #assume at most six neighbors 
    nwet = sum(γ.wet)
    wet = γ.wet
    counter = 0
    times_counter = 1

    m_new = [similar(m1, eltype(m)) for m1 in m_template]
    # m_new = deepcopy(m_template)
    # Copy template values into the new arrays
    for (i, mf_new) in enumerate(m_new)
        mf_new.fraction .= m_template[i].fraction
    end
    m_length = [1]
    for m1 in m_new
        push!(m_length, length(m1))
    end

    m_length = cumsum(m_length)
    for m1 in m_new
        m1.fraction[m1.γ.wet] .= m[m_length[times_counter]:m_length[times_counter+1]-1]
        times_counter += 1
        counter +=nwet
    end
    A = watermassmatrix(m_new, γ)
    return sum(A.^2)
end

function Jsumquare(m)
    A = watermassmatrix(m, γ)
    sum(A.^2)
end

function gJsumsquare(gJ, m)
    nfield = sum(γ.wet)
    nm = 0
    for m1 in m
        nm += length(m1)
    end

    ilist = Array{Int64}(undef,nm+nfield) #n_neighbors + diagnonal terms
    jlist = Array{Int64}(undef,nm+nfield)
    mlist = Array{Float64}(undef,nm+nfield)

        # ones on diagonal
    for ii in 1:nfield
        ilist[ii] = ii
        jlist[ii] = ii
        mlist[ii] = 0.0
    end
    

    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    # allocate water mass matrix
    nfield = sum(γ.wet)

    A = watermassmatrix(m, γ)
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
                mlist[counter] = 2 * A[nrow,R[Istep]] * gJ #no minus sign because its just entries of A (don't take m into account)
                #A[nrow,R[Istep]] = -m1.fraction[I]
            end
        end
    end

    gA = sparse(ilist,jlist,mlist)
    gm = gwatermassmatrix(gA, m, γ)
    return gm

end


m_vec = vcat([vec(m1) for m1 in m][:]...)
i, j, mA = findnz(A)

Jsumquare(m_vec, γ, m)
Jsumquare(m)

using ReverseDiff
using BenchmarkTools
dJ = ReverseDiff.gradient(mvec -> Jsumquare(mvec, γ, m), m_vec)
dJ_adj = gJsumsquare(1., m)
dJ_adj_vec =  vcat([vec(m1) for m1 in dJ_adj][:]...)

@btime ReverseDiff.gradient(mvec -> Jsumquare(mvec, γ, m), m_vec);
@btime gJsumsquare(1., m);

@btime watermassmatrix(m, γ);

@btime gwatermassmatrix(A, m, γ);

i, j, z = findnz(A)