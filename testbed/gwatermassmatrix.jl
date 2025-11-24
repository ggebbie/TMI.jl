import Pkg; Pkg.activate(".")

using Revise
using TMI
using SparseArrays

import TMI:MassFraction
import TMI:step_cartesian
import TMI:wet
import TMI:unvec
import TMI:unvec!


# Copy similar that preserves tracked types
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

function Base.similar(nt::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    return map(similar, nt)
end


vec(m::MassFraction) = m.tracer[m.γ.wet]
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



TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion);
m = massfractions_isotropic(γ)
mvec = vec(m)

tmp = unvec(m, mvec)[:up].fraction .- m[:up].fraction


tmp = 0.0 .* m[:down].γ.wet
[tmp .+= m[k].γ.wet for k in keys(m)]
sum(tmp .> 0.0)
sum(γ.wet)
sum(γ.interior)

sum(filter(!isnan, tmp))

cartesianindex(m[:down].γ.wet)


Base.vec(m::MassFraction) = m.fraction[m.γ.wet]

mvec = vec(m[:down])



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

function gsteadyinversion(gc::Field,c::Field,A, Alu, b::Union{BoundaryCondition,NamedTuple},γ::Grid;q=nothing,r=1.0) #where T <: Real
    #println("running adjoint steady inversion")
    gd = Alu' \ gc

    #ga = gd * c', but sum A entries are fixed (1, or 0)
    #just need to consider the rows/columns that propogate into the adjoint
    rows, cols, _ = findnz(A) 

    # only consider the off-diagonal entries
    rows_offdiag = rows[offdiag_mask]
    offdiag_mask = rows .!= cols 
    cols_offdiag = cols[offdiag_mask]
    vals_offdiag = -vec(dbar)[rows_offdiag] .* vec(c)[cols_offdiag]
    gA = sparse(rows_offdiag, cols_offdiag, vals_offdiag, size(A,1), size(A,2))

    gb = gsetboundarycondition(gd,b)

    if !isnothing(q)
        gq = gsetsource(gd,q,r)
        return gb,gq, gA
    else
        return gb, gA
    end

end




function prior_boundary_cost(uvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return uvec' * Qⁱᵤ * uvec #a shortcut since (u₀ + u) - u₀ = u
end

function gprior_boundary_cost(uvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    gJ = 1 #should always be 1
    guvec = 2 * Qⁱᵤ * uvec * gJ 
    return guvec
end

function prior_source_cost(qvec::Vector,q₀::Union{BoundaryCondition,NamedTuple},Qⁱₛ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return qvec' * Qⁱₛ * qvec #a shortcut since (u₀ + u) - u₀ = u
end

function gprior_boundary_cost(qvec::Vector,q₀::Union{BoundaryCondition,NamedTuple},Qⁱₛ::Symmetric)
    gJ = 1 #should always be 1
    gqvec = 2 * Qⁱₛ * qvec * gJ 
    return gqvec
end


function steadyinversion_residual(c::Field, c_obs::Union{Vector, Field}, 
                                  locs::{Nothing, Vector{Tuple{T,T,T}}}, 
                                  γ::Grid)::Vector{T} where T <: Real
    if c_obs isa Field
        r = c - c_obs
        return r
    elseif c_obs isa Vector
        y = observe(c, locs, γ)
        r = y .- c_obs
        return r
    end
end


function gsteadyinversion_residual(dr::Union{Vector, Field},
                                    c::Field,
                                  locs::Vector{Tuple{T,T,T}}, 
                                  γ::Grid)::Vector{T} where T <: Real
    if dr isa Field
        gc = dr
        return gc
    elseif dr isa Vector
        gy = dr 
        gc = gobserve(gy,c,locs) #this recalculates weights everytime its called. could be improved
        return gc
    end
end

function model_observation_cost(r::Union{Vector, Field},Wⁱ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return r ⋅ (Wⁱ * r) # dot product, gives total misfit

end

function gmodel_observation_cost(r::Vector,Wⁱ::Symmetric)
    gJ = 1 #should always be 1
    grvec = 2 * Wⁱ * r * gJ 
    return grvec
end




TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion);
m = massfractions_isotropic(γ)

y = (θ =  readfield(TMIfile, "θ", γ), 
    θ2 =  readfield(TMIfile, "θ", γ))

A * y

bθ = getsurfaceboundary(y.θ)


nb = length(vec(bθ))
X = diagm(ones(nb))
(X * bθ) - bθ

using BenchmarkTools

f1(bθ, X) = vec(bθ)' * (X *  vec(bθ))
f2(bθ, X) = vec(bθ)' * X *  vec(bθ)


@btime dot(vec(bθ), X, vec(bθ))
@btime f1(bθ, X) 
@btime f2(bθ, X) 