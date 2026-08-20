"""
    Observations(y, σ, γ; locations=nothing, L=nothing, Q=nothing)

Store tracer observations and their statistical description on a grid.

# Arguments
- `y::NamedTuple`: tracer `Field`s or vectors of observed values
- `σ`: scalar, `Field`, vector, or tracer `NamedTuple` of uncertainties
- `γ::Grid`: model grid
- `locations`: optional coordinate vector or tracer `NamedTuple` of coordinates
- `L`: optional scalar, `Field`, vector, or tracer `NamedTuple` of decorrelation scales
- `Q`: optional `AbstractMatrix` or tracer `NamedTuple` of precision matrices

# Output
- `observations::Observations`: normalized observations and block precision matrix
"""
struct Observations{ObservationValues,Uncertainties,LengthScales,Locations,
    InterpolationIndices,GridType,PrecisionMatrix}
    y::ObservationValues
    σ::Uncertainties
    L::LengthScales
    locations::Locations
    interpolation_indices::InterpolationIndices
    γ::GridType
    Q::PrecisionMatrix
    function Observations(y::NamedTuple{names,Y}, σ, γ::Grid;
        locations=nothing, L=nothing, Q=nothing,
    ) where {names,Y<:Tuple{Vararg{Union{Field,AbstractVector}}}}
        if σ isa NamedTuple
            σ = NamedTuple{keys(y)}(map(name -> σ[name], keys(y)))
        else
            σ = map(_ -> σ, y)
        end
        L, locations, supplied_Q = map(specification -> begin
            if specification isa NamedTuple
                NamedTuple{keys(y)}(map(name -> get(specification, name, nothing),
                    keys(y)))
            else
                map(_ -> specification, y)
            end
        end, (L, locations, Q))
        indices = map(support -> begin
            if isnothing(support)
                nothing
            else
                interpindex(support, γ)
            end
        end, locations)
        precisions = map((observed, uncertainty, scale, indices, precision) -> begin
            if !isnothing(precision)
                sparse(precision)
            elseif isnothing(scale)
                if uncertainty isa Real
                    number_of_values = observed isa Field ? (isnothing(indices) ?
                        sum(γ.interior) : length(indices)) : length(observed)
                    uncertainty_values = fill(uncertainty, number_of_values)
                elseif uncertainty isa Field && isnothing(indices)
                    uncertainty_values = uncertainty.tracer[γ.interior]
                elseif uncertainty isa Field
                    uncertainty_values = observe(uncertainty, indices, γ)
                else
                    uncertainty_values = uncertainty
                end
                sparse(Diagonal(inv.(uncertainty_values .^ 2)))
            else
                if uncertainty isa Real
                    uncertainty_values = fill(uncertainty, sum(γ.interior))
                elseif uncertainty isa Field
                    uncertainty_values = uncertainty.tracer[γ.interior]
                else
                    uncertainty_values = uncertainty
                end
                if scale isa Real
                    scale_values = fill(scale, sum(γ.interior))
                elseif scale isa Field
                    scale_values = scale.tracer[γ.interior]
                else
                    scale_values = scale
                end
                gaussianprecision(uncertainty_values, scale_values, γ)
            end
        end, y, σ, L, indices, supplied_Q)
        precision = blockdiag(Base.values(precisions)...)
        return new{typeof(y),typeof(σ),typeof(L),typeof(locations),
            typeof(indices),typeof(γ),typeof(precision)}(y, σ, L, locations,
            indices, γ, precision)
    end
end

function vec(observations::Observations)
    γ = observations.γ
    values = map((value, indices) -> begin
        if !(value isa Field)
            value
        elseif isnothing(indices)
            value.tracer[γ.interior]
        else
            observe(value, indices, γ)
        end
    end, observations.y, observations.interpolation_indices)
    return vec(values)
end

"""
    steadyinversion(Alu, b, q, stoichiometry, γ)
    steadyinversion(Alu, b, q, observations; stoichiometry)

Solve multiple steady tracers and optionally sample them as observations.

# Arguments
- `Alu`, `b`, `q`: factored transport, boundaries, and sources
- `stoichiometry::NamedTuple`: stoichiometric dependencies in the coupled tracer equations
- `γ`, `observations`: model grid and observation sampling

# Output
- `modeled`: modeled tracer fields or observed values
"""
function steadyinversion(Alu, b::NamedTuple, q::NamedTuple,
    stoichiometry::NamedTuple, γ::Grid)
    tracer_names = keys(stoichiometry)
    boundaries = NamedTuple{tracer_names}(map(name -> b[name], tracer_names))
    sources = map(coefficients -> begin
        if isnothing(coefficients)
            nothing
        else
            reduce(+, map(name -> coefficients[name] * q[name], keys(coefficients)))
        end
    end, stoichiometry)
    return steadyinversion(Alu, boundaries, sources, γ)
end
function steadyinversion(Alu, b::NamedTuple, q::NamedTuple,
    observations::Observations; stoichiometry::NamedTuple)
    modeled = steadyinversion(Alu, b, q, stoichiometry, observations.γ)
    γ = observations.γ
    return map((field, indices) -> begin
        isnothing(indices) ? field.tracer[γ.interior] : observe(field, indices, γ)
    end, modeled, observations.interpolation_indices)
end

function Base.show(io::IO, observations::Observations)
    gridded = count(isnothing, observations.interpolation_indices)
    print(io, "Observations($(length(observations.y)) tracers; ",
        "$gridded gridded, $(length(observations.y) - gridded) sparse; ",
        "$(length(vec(observations))) values)")
end
