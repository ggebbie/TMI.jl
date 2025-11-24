import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI

TMIversion = versionlist()[6] # G14 has no remote mass fractions
A, Alu, γ, TMIfile, L, B = config(TMIversion);

γ.I

neighbors
function watermassmatrix2massfraction(A, γ::Grid)
        #adhoc, could also use neighbor_indices(n::Integer)
        m = (north = massfractions_north(A,γ),
        east   = massfractions_east(A,γ),
        south  = massfractions_south(A,γ),
        west   = massfractions_west(A,γ),
        up     = massfractions_up(A,γ),
        down   = massfractions_down(A,γ))

        #add a check for bottom spreading param (i.e., ) A ≠ watermassmatrix(m, γ) for 
        # the newly formed m. 
        # could also use a check based on neighbors(m::Union{Vector,NamedTuple}, γ::Grid{R,N})

end


function update_inversion(cobs::NamedTuple, W⁻¹::NamedTuple{Matrix}, A_prior, γ, ; 
    update_mass_fraction = true, #variant of massfractions
    update_boundary_conditions = true, #variant of sparse data map, 
    src_index::Union{Nothing, Vector{CartesianIndex{3}}} = nothing,
    dst_index::Union{Nothing, Vector{CartesianIndex{3}}} = nothing
    )
    ntracers = length(cobs)
    m = watermassmatrix2massfraction(A, γ) #::Union{Vector,NamedTuple}
    m = ()

end


function get_interpolation_weights(locs::Vector{Tuple{Float64,Float64,Float64}}, γ::Grid)
        N = length(locs)
        wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
        [wis[i] = interpindex(locs[i],γ) for i in 1:N]
        return wis
end

function massfractions_cost()

end
function nonconservativetracer_cost()

end
function model_data_misfit_cost(cobs::Union{Field, Vector}, Alu, W⁻¹::NamedTuple{Matrix}, 
                           b₀::NamedTuple{BoundaryCondition}; 
                           u::Vector = nothing, 
                           cobs_locs::Vector{Tuple{Float64,Float64,Float64}} = nothing)
    #W⁻¹ is the tracer weighting 
    #COST wrt input boundary conditions and mass matrix  
    # # get weighted interpolation indices
    J = 0.0
    # cobs_adjusted = adjustsource(q,cobs) need to adjust for sources/sinks? 

    for varname in keys(cobs)
        #b += u # easy case where u and b are on the same boundary
        b = isnothing(u) ? b₀ : adjustboundarycondition(b₀,u) 
        ϕmodel = steadyinversion(Alu,b,γ)
        cobsvar = cobs[varname]

        if typeof(cobsvar) isa Field 
            n =  vec(ϕtrue - cobsvar)  # gives the misfit. could also use observe, but might be slower
        elseif typeof(cobsvar) isa Vector  
            wis = get_interpolation_weights(cobs_locs, γ)
            ϕmodel_sample = observe(ϕmodel,wis,γ)
            n = ϕmodel_sample .- cobsvar

        end
        J += (n' * W⁻¹ * n)

        #g = square_matrix(length of cobs * (length of cobsvar + length of u ))
        #gradient calculated for each tracer. 
        #    dc = 2*W{nc}*misfit; % forcing by tracer misfit.
        # gc[tracer_silver] = R' \ (P' * (L' \ (U' \ (Q' * dc)))); 
    end

    return J
end