

"""
    check_shared_references(nt::NamedTuple, name::String)

Check if any fields in a `NamedTuple` point to the same object in memory. This
is a utility to prevent bugs where in-place updates to one control parameter
accidentally modify another. Throws an error if shared references are found.
"""
function check_shared_references(nt::NamedTuple, name::String)
    vals = filter(!isnothing, collect(values(nt)))
    n = length(vals)

    for i in 1:n
        for j in (i+1):n
            if vals[i] === vals[j]
                key_names = collect(keys(nt))
                # Find which keys share the reference
                shared_keys = String[]
                for k in 1:length(nt)
                    if !isnothing(nt[k]) && (nt[k] === vals[i])
                        push!(shared_keys, string(key_names[k]))
                    end
                end
                error("Shared reference in '$name'. Fields $(join(shared_keys, ", "))\n" *
                      "point to same object. Use deepcopy for independent copies.")
            end
        end
    end
end
check_shared_references(x, name::String) = nothing  # Non-NamedTuple, no check needed


"""
    _assemble_precision_matrix(val, len; name="<unknown>") -> AbstractMatrix

Builds a precision matrix (inverse covariance) from a configuration value.

# Arguments
- `val`: Scalar variance or covariance matrix.
- `len`: Dimension of the control vector.
- `name`: (Optional) Variable name for error reporting.
"""
function _assemble_precision_matrix(val::Union{Real, AbstractVector, Symmetric, Diagonal}, len::Int; name::String="<unknown>")
    if val isa Real
        return Diagonal((1.0 / val) .* ones(len))
    elseif val isa AbstractVector
        length(val) != len && error("Vector variance '$name' length mismatch. Exp: $len, Got: $(length(val)).")
        return Diagonal(1.0 ./ val)
    elseif val isa Union{Symmetric, Diagonal}
        return inv(val)
    else
        error("Invalid type for '$name' uncertainty config: $(typeof(val)).")
    end
end

"""
    _build_tracer_precision_matrix(config::NamedTuple, controls_nt::NamedTuple) -> NamedTuple

Constructs precision matrices for tracer-like control variables.

# Arguments
- `config`: Configuration for a specific control component.
- `controls_nt`: `NamedTuple` of control variables.
"""
function _build_tracer_precision_matrix(config::NamedTuple, controls_nt::NamedTuple)
    Q_config = get(config, :covariance, get(config, :variance, nothing))
    ctrl_keys = filter(k -> !isnothing(controls_nt[k]), keys(controls_nt))

    if isnothing(Q_config)
        @warn "No uncertainty config. Defaulting components to variance 1.0."
        return NamedTuple{ctrl_keys}(map(name -> Diagonal(ones(length(vec(controls_nt[name])))), ctrl_keys))
    end

    return NamedTuple{ctrl_keys}(
        map(ctrl_keys) do name
            len = length(vec(controls_nt[name]))
            val = get(Q_config, name, 1.0)
            !haskey(Q_config, name) && @warn "No uncertainty config for '$name'. Defaulting to variance 1.0."
            return _assemble_precision_matrix(val, len; name=string(name))
        end
    )
end

"""
    _build_massfrac_precision_matrix(config::NamedTuple, controls_nt::NamedTuple) -> AbstractMatrix

Constructs the precision matrix for mass fraction control variables.

# Arguments
- `config`: `mass_fraction` configuration tuple.
- `controls_nt`: `NamedTuple` of `MassFraction` control variables.
"""
function _build_massfrac_precision_matrix(config::NamedTuple, controls_nt::NamedTuple)
    Q_config = get(config, :covariance, get(config, :variance, nothing))
    
    if isnothing(Q_config)
        @warn "No uncertainty config for mass_fraction. Defaulting variance 1.0."
        return Diagonal(ones(length(vec(controls_nt))))
    end

    if Q_config isa Real || Q_config isa Union{AbstractVector, Symmetric, Diagonal}
        return _assemble_precision_matrix(Q_config, length(vec(controls_nt)); name="mass_fraction")
    
    elseif Q_config isa NamedTuple
        vars_vecs = map(keys(controls_nt)) do name
            len = length(vec(controls_nt[name]))
            val = get(Q_config, name, 1.0)
            !haskey(Q_config, name) && @warn "No uncertainty config for mass fraction '$name'. Defaulting variance 1.0."
            !(val isa Real) && error("Mass fraction NamedTuple variance: values must be scalar.")
            (1.0 / val) .* ones(len)
        end
        return Diagonal(vcat(vars_vecs...))
    else
        error("Invalid type for mass_fraction uncertainty config: $(typeof(Q_config)).")
    end
end

# Helper to generate a bound NamedTuple for a given control set
function _generate_control_bounds(controls_nt::NamedTuple, raw_bounds::Union{NamedTuple, Nothing}, default_val::Real)
    map(keys(controls_nt)) do name
        ctrl = getproperty(controls_nt, name)
        if !isnothing(ctrl)
            val = isnothing(raw_bounds) ? default_val : get(raw_bounds, name, default_val)
            return vec(ctrl) .* 0.0 .+ val # Create a Vector of correct size filled with the bound
        else
            return nothing # If control is nothing, bound is nothing
        end
    end |> NamedTuple{keys(controls_nt)} # Ensure NamedTuple keys are preserved
end