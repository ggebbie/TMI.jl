#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example : Find the distribution of a tracer given:              %
%       (a) the pathways described by A,                          %
%       (b) interior sources and sinks given by dC,               % 
%           that best fits observations, Cobs,                    %
%   and (c) inequality constraints on the tracer concentration.   %
%                                                                  %
% Mathematically, minimize J = (C-Cobs)^T W (C-Cobs) subject to    %
%                         AC = d + Gamma u                         %
%  where u is the estimated change in surface concentration.    % 
%
% See Supplementary Section 2, Gebbie & Huybers 2011.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

using Revise
using TMI
using SparseArrays
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

m = massfractions_isotropic(γ)

#
function Jsumquare(m::Vector, γ, m_template)
    #assume at most six neighbors 
    nwet = sum(γ.wet)
    wet = γ.wet
    counter = 0
    times_counter = 1

    for m1 in m_new
        m1.fraction[m1.γ.wet] .= m[m_length[times_counter]:m_length[times_counter+1]]
        times_counter += 1
        counter +=nwet
    end
    A = watermassmatrix(m, γ)
    sum(A.^2)
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
    gm = deepcopy(m)
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
                #A[nrow,R[Istep]] = -m1.fraction[I]
            end
        end
    end

    return gm
end



function Jsumquare(m::Vector, γ, m_template)
    #assume at most six neighbors 
    nwet = sum(γ.wet)
    wet = γ.wet
    counter = 0
    times_counter = 1
    m_new = deepcopy(m_template)

    m_length = [1]
    for m1 in m_new
        push!(m_length, length(m1))
    end
    m_length = cumsum(m_length)
    for m1 in m_new
        # print(times_counter)
        m1.fraction[m1.γ.wet] .= m[m_length[times_counter]:m_length[times_counter+1]-1]
        times_counter += 1
        counter +=nwet
    end
    A = watermassmatrix(m_new, γ)
    return sum(A.^2)
end

function Jsumquare(m::NamedTuple)

    A = watermassmatrix(m, γ)
    sum(A.^2)
end

function gJsumsquare(gJ, m::NamedTuple)
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
                mlist[counter] = -2 * gJ
                #A[nrow,R[Istep]] = -m1.fraction[I]
            end
        end
    end

    gA = sparse(ilist,jlist,mlist)
    gm = gwatermassmatrix(gA, m, γ)
    return gm

end


Jsumquare(m)
m_vec = vcat([vec(m1) for m1 in m][:]...)
Jsumquare(m_vec, γ, m)


# using ForwardDiff
using Zygote

dJ = Zygote.gradient(mvec -> Jsumquare(mvec, γ, m), m_vec)

dJ_adj = gJsumsquare(1., m)

