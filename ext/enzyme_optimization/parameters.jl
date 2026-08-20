"""
    Inversion(observations; boundary_conditions, sources, stoichiometry,
        mass_fractions)

Construct an inversion from observations and fixed or controlled parameters.

# Arguments
- `observations::Observations`: observed values and precision
- `boundary_conditions::NamedTuple`, `sources::NamedTuple`: fixed controls or `(prior, σ, L, lb, ub)`
- `stoichiometry::NamedTuple`: stoichiometric dependencies in the coupled tracer equations
- `mass_fractions::NamedTuple`: fixed fractions or `(prior=m₀, σ=σ_m, lb=..., ub=...)`

# Output
- `inversion::Inversion`: physical controls, bounds, priors, and precision matrices
"""
struct Inversion{ObservationData,Stoichiometry,ControlVector,ControlBounds,
    ControlRanges,ControlTemplates,BoundaryIndices,BoundaryValues,ControlPriors,
    PriorVector,BoundaryPrecision,SmoothnessPrecision,SourcePrecision,
    MassFractionPrecision}
    observations::ObservationData
    stoichiometry::Stoichiometry
    x₀::ControlVector
    bounds::ControlBounds
    control_ranges::ControlRanges
    control_templates::ControlTemplates
    boundary_control_indices::BoundaryIndices
    boundary_template_values::BoundaryValues
    control_priors::ControlPriors
    x_prior::PriorVector
    Q_b::BoundaryPrecision
    Q_s::SmoothnessPrecision
    Q_q::SourcePrecision
    Q_m::MassFractionPrecision
end
_controlvalues(template, value::Real) = vec(value * one(template))
_controlvalues(::Any, values) = vec(values)
_specificationvalues(::NamedTuple{(),Tuple{}}, ::Symbol, ::Real) = Float64[]
_specificationvalues(settings::NamedTuple, name::Symbol, default::Real) = vcat(
    map(setting -> _controlvalues(setting.prior, get(setting, name, default)), settings)...)

function _splitcontrols(entries::NamedTuple)
    names = Tuple(keys(entries))
    controlled_names = Tuple(filter(name -> entries[name] isa NamedTuple, names))
    controlled = NamedTuple{controlled_names}(map(name -> entries[name], controlled_names))
    templates = map(entry -> entry isa NamedTuple ? entry.prior : entry, entries)
    return controlled, templates
end
function _independentprecision(template, uncertainty)
    σ = _controlvalues(template, uncertainty)
    return sparse(Diagonal(inv.(σ .^ 2)))
end
_blockdiagonal(::NamedTuple{(),Tuple{}}) = spzeros(Float64, 0, 0)
_blockdiagonal(blocks::NamedTuple) = blockdiag(Base.values(blocks)...)

function Inversion(observations::Observations; boundary_conditions::NamedTuple,
    sources::NamedTuple, stoichiometry::NamedTuple, mass_fractions::NamedTuple)
    boundary_settings, boundary_templates = _splitcontrols(boundary_conditions)
    source_settings, source_templates = _splitcontrols(sources)
    b₀ = map(setting -> setting.prior, boundary_settings)
    q₀ = map(setting -> setting.prior, source_settings)
    source_templates = map(source -> adjustsource(source, zero(source)), source_templates)
    if haskey(mass_fractions, :prior)
        prior = mass_fractions.prior
        mass = (prior=prior, template=prior, initial=vec(prior),
            lower=_controlvalues(prior, get(mass_fractions, :lb, -Inf)),
            upper=_controlvalues(prior, get(mass_fractions, :ub, Inf)),
            precision=_independentprecision(prior, mass_fractions.σ))
    else
        mass = (prior=nothing, template=mass_fractions, initial=Float64[],
            lower=Float64[], upper=Float64[], precision=spzeros(0, 0))
    end

    b₀_values, q₀_values = vec(b₀), vec(q₀)
    x₀ = vcat(b₀_values, zero(q₀_values), mass.initial)
    number_of_b, number_of_q = length(b₀_values), length(q₀_values)
    control_ranges = (boundary=1:number_of_b,
        source=(number_of_b + 1):(number_of_b + number_of_q),
        mass_fraction=(number_of_b + number_of_q + 1):length(x₀))

    b_lower = _specificationvalues(boundary_settings, :lb, -Inf)
    b_upper = _specificationvalues(boundary_settings, :ub, Inf)
    q_lower = _specificationvalues(source_settings, :lb, -Inf) .- q₀_values
    q_upper = _specificationvalues(source_settings, :ub, Inf) .- q₀_values
    lower_bounds = vcat(b_lower, q_lower, mass.lower)
    upper_bounds = vcat(b_upper, q_upper, mass.upper)

    σ_b = map(setting -> _controlvalues(setting.prior, setting.σ), boundary_settings)
    σ_b = map((setting, values) -> unvec(setting.prior, values),
        boundary_settings, σ_b)
    Q_b = map(_independentprecision, b₀, map(setting -> setting.σ, boundary_settings))
    laplacian = surfacelaplacianmatrix(observations.γ)
    L_b = map(setting -> begin
        scale = get(setting, :L, nothing)
        if scale isa Real
            scale * one(setting.prior)
        elseif scale isa Field
            getsurfaceboundary(scale)
        else
            scale
        end
    end, boundary_settings)
    Q_s = map((σ, scale) -> begin
        if isnothing(scale)
            spzeros(eltype(vec(σ)), length(σ), length(σ))
        else
            boundary_smoothness_precision_matrix(σ, scale, laplacian)
        end
    end, σ_b, L_b)
    Q_q = map(_independentprecision, q₀, map(setting -> setting.σ, source_settings))

    tracer_names = keys(observations.y)
    stoichiometry = NamedTuple{tracer_names}(map(
        name -> get(stoichiometry, name, nothing), tracer_names))
    boundary_control_indices, offset = Int[], 0
    for name in keys(boundary_templates)
        count = length(boundary_templates[name])
        if haskey(boundary_settings, name)
            append!(boundary_control_indices, offset .+ (1:count))
        end
        offset += count
    end
    control_templates = (boundary_conditions=boundary_templates, sources=source_templates,
        mass_fractions=mass.template)
    control_priors = (boundary=b₀, source=q₀, mass_fraction=mass.prior)
    return Inversion(observations, stoichiometry, x₀,
        (lower=lower_bounds, upper=upper_bounds), control_ranges,
        control_templates, boundary_control_indices, vec(boundary_templates),
        control_priors, copy(x₀), _blockdiagonal(Q_b),
        _blockdiagonal(Q_s), _blockdiagonal(Q_q), mass.precision)
end
function Base.show(io::IO, inversion::Inversion)
    print(io, "Inversion($(length(inversion.x₀)) controls, ",
        "$(length(inversion.observations.y)) observed tracers)")
end
unvec(::Nothing, template::NamedTuple, ::AbstractVector) = template
unvec(::NamedTuple, template::NamedTuple, x::AbstractVector) = unvec(template, x)
unvec!(::NamedTuple, ::AbstractVector, ::Nothing) = nothing
unvec!(mass_fractions::NamedTuple, x::AbstractVector, ::NamedTuple) =
    unvec!(mass_fractions, x)

"""
    unvec(inversion::Inversion, x)

Reconstruct boundary conditions, sources, and mass fractions from controls.

# Arguments
- `inversion`, `x`: inversion description and optimizer control vector

# Output
- `controls`: reconstructed boundary conditions, sources, and mass fractions
"""
function unvec(inversion::Inversion, x::AbstractVector)
    ranges, templates = inversion.control_ranges, inversion.control_templates
    boundary_values = copy(inversion.boundary_template_values)
    boundary_values[inversion.boundary_control_indices] .= x[ranges.boundary]
    boundary_conditions = unvec(templates.boundary_conditions, boundary_values)
    sources = adjustsource(inversion.control_priors.source,
        x[ranges.source], templates.sources)
    mass_fractions = unvec(inversion.control_priors.mass_fraction,
        templates.mass_fractions, x[ranges.mass_fraction])
    return (; boundary_conditions, sources, mass_fractions)
end
