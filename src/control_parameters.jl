"""
    check_shared_references(nt::NamedTuple, name::String)

Error if non-nothing fields in a NamedTuple share the same memory address. Use to
catch cases where two controls accidentally reference the same object, causing
in-place updates to affect multiple fields.

# Arguments
- `nt`: NamedTuple to check
- `name`: Name used in the error message
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
                error("Shared reference detected in $name: fields $(join(shared_keys, ", ")) " *
                      "point to the same object in memory. Use deepcopy to create independent copies.")
            end
        end
    end
end
check_shared_references(x, name::String) = nothing  # Non-NamedTuple, no check needed

"""
    ControlParameters

Container for optimization control variables (boundary conditions, sources, mass fractions) plus optional first guesses and weighting matrices. Supports caching a flattened control vector and lengths for slicing/unvec operations.

# Fields (all may be `nothing`)
- `γ`: grid (required)
- `ub`, `uq`, `m`: NamedTuples of boundary conditions, sources, and mass fractions
- `u₀`, `q₀`, `m₀`: first-guess values for the above
- `Qᵤ`, `Qₛ`, `Qₘ`: inverse covariance matrices for control errors
- `ub_bounds`, `uq_bounds`, `m_bounds`: optional bounds (NamedTuples keyed like `ub`/`uq`/`m`) for box constraints
- `vector`: cached flattened controls (if `cache_vector=true`)
- `ub_length`, `uq_length`, `m_length`: cached segment lengths (always stored)
- `mass_fraction_steps`: cached precomputed mass-fraction step targets (if enabled)
- `A`: cached water-mass matrix (if enabled)
- `dub`, `duq`: cached boundary/source perturbations
- `gm`: cached mass-fraction gradient buffer
- `b`, `q`: cached boundary/source fields initialized from priors
"""
struct ControlParameters{UB,UQ,M,U0,Q0,M0,QU,QS,QM,UBB,UQB,MB,V,UBL,QL,ML,S,AT,G<:Grid,DUB,DUQ,GDUB,GDUQ,GM,BU,BQ}
    ub::UB; uq::UQ; m::M
    u₀::U0; q₀::Q0; m₀::M0
    Qᵤ::QU; Qₛ::QS; Qₘ::QM
    ub_bounds::UBB; uq_bounds::UQB; m_bounds::MB
    vector::V
    ub_length::UBL; uq_length::QL; m_length::ML
    mass_fraction_steps::S
    A::AT
    γ::G
    dub::DUB
    duq::DUQ
    gdub::GDUB
    gduq::GDUQ
    gm::GM
    b::BU
    q::BQ
end

check_bounds_keys(::Nothing, ::Any, ::String) = nothing
function check_bounds_keys(bounds::NamedTuple, controls::NamedTuple, name::String)
    keys(bounds) == keys(controls) || error("$name keys must match control keys")
    return nothing
end


function ControlParameters(; γ::Grid,
                            ub = nothing, uq = nothing, m = nothing,
                            u₀ = nothing, q₀ = nothing, m₀ = nothing,
                            Qᵤ = nothing, Qₛ = nothing, Qₘ = nothing,
                            ub_bounds = nothing, uq_bounds = nothing, m_bounds = nothing,
                            cache_vector = false,
                            cache_mass_fraction_steps = true,
                            cache_A = true)
    # Ensure control entries and priors are independent objects.
    for (val, name) in ((ub, "ub"), (uq, "uq"), (m, "m"),
                        (u₀, "u₀"), (q₀, "q₀"), (m₀, "m₀"),
                        (Qᵤ, "Qᵤ"), (Qₛ, "Qₛ"))
        check_shared_references(val, name)
    end

    # Validate optional bounds.
    !isnothing(ub_bounds) && isnothing(ub) && throw(ArgumentError("ub_bounds provided but ub is nothing"))
    !isnothing(uq_bounds) && isnothing(uq) && throw(ArgumentError("uq_bounds provided but uq is nothing"))
    !isnothing(m_bounds) && isnothing(m) && throw(ArgumentError("m_bounds provided but m is nothing"))
    isnothing(ub) || check_bounds_keys(ub_bounds, ub, "ub_bounds")
    isnothing(uq) || check_bounds_keys(uq_bounds, uq, "uq_bounds")
    isnothing(m) || check_bounds_keys(m_bounds, m, "m_bounds")

    # Optionally cache precomputed mass-fraction steps and water-mass matrix.
    m_steps = nothing
    if cache_mass_fraction_steps
        isnothing(m) && throw(ArgumentError("cache_mass_fraction_steps=true requires m"))
        m_steps = precompute_mass_fraction_steps(m, γ)
    end

    A_cached = nothing
    if cache_A
        isnothing(m) && throw(ArgumentError("cache_A=true requires m"))
        if cache_mass_fraction_steps
            A_cached = watermassmatrix(m, γ, m_steps)
        else
            A_cached = watermassmatrix(m, γ)
        end
    end

    # Optionally cache flattened controls for Optim/ForwardDiff style interfaces.
    vector = cache_vector ? vectorize_controls(ub, uq, m) : nothing

    # Segment lengths let `unvec` slice an incoming control vector without recomputing sizes.
    ub_length = isnothing(ub) ? nothing : length(vec(ub))
    uq_length = isnothing(uq) ? nothing : length(vec(uq))
    m_length = isnothing(m) ? nothing : length(vec(m))

    dub_cached = isnothing(ub) ? nothing : deepcopy(ub)
    duq_cached = isnothing(uq) ? nothing : deepcopy(uq)

    gdub_cached = isnothing(ub) ? nothing : deepcopy(ub)
    gduq_cached = isnothing(uq) ? nothing : deepcopy(uq)
    gm_cached = isnothing(m) ? nothing : deepcopy(m)

    b_cached = isnothing(u₀) ? nothing : deepcopy(u₀)
    q_cached = isnothing(q₀) ? nothing : deepcopy(q₀)
    return ControlParameters{typeof(ub), typeof(uq), typeof(m),
                             typeof(u₀), typeof(q₀), typeof(m₀),
                             typeof(Qᵤ), typeof(Qₛ), typeof(Qₘ),
                             typeof(ub_bounds), typeof(uq_bounds), typeof(m_bounds),
                             typeof(vector),
                             typeof(ub_length), typeof(uq_length), typeof(m_length),
                             typeof(m_steps), typeof(A_cached),
                             typeof(γ), typeof(dub_cached), typeof(duq_cached),
                             typeof(gdub_cached), typeof(gduq_cached),
                             typeof(gm_cached),
                             typeof(b_cached), typeof(q_cached)}(
        ub, uq, m, u₀, q₀, m₀, Qᵤ, Qₛ, Qₘ, ub_bounds, uq_bounds, m_bounds,
        vector, ub_length, uq_length, m_length,
        m_steps, A_cached, γ, dub_cached, duq_cached, gdub_cached, gduq_cached, gm_cached, b_cached, q_cached)
end


"""
    vectorize_controls(ub, uq, m) -> Vector

Flatten control parameters into a single vector by concatenating vectorized non-nothing fields.
"""
function vectorize_controls(ub::NamedTuple, uq::NamedTuple, m::NamedTuple)
    return vcat(vec.([ub, uq, m])...)
end

"""
    vectorize_controls(cp::ControlParameters) -> Vector

Extract or compute the flattened vector representation from a ControlParameters object.
If the vector is already stored, return it. Otherwise, compute it from the control fields.
"""
function vec(cp::ControlParameters)
    if !isnothing(cp.vector)
        return cp.vector
    else
        return vectorize_controls(cp.ub, cp.uq, cp.m)
    end
end

"""
    unvec(cp, vector) -> (ub, uq, m)

Reconstruct boundary condition, source, and mass fraction NamedTuples from a flat control vector.
"""
function unvec(cp::ControlParameters, vector::AbstractVector)
    # Calculate boundaries
    u_len = cp.ub_length
    q_len = cp.uq_length
    m_len = cp.m_length
    
    idx = 1
    ub = nothing; uq = nothing; m = nothing
    
    if !isnothing(cp.ub)
        ub = unvec(cp.ub, vector[idx:idx+u_len-1])
        idx += u_len
    end
    
    if !isnothing(cp.uq)
        uq = unvec(cp.uq, vector[idx:idx+q_len-1])
        idx += q_len
    end
    
    if !isnothing(cp.m)
        m = unvec(cp.m, vector[idx:idx+m_len-1])
    end

    return ub, uq, m
end

"""
    unvec!(cp, vector) -> (ub, uq, m)

In-place version of `unvec` that writes values from `vector` directly into the
control containers stored on `cp` (`cp.ub`, `cp.uq`, `cp.m`) without creating slices.
Only non-`nothing` fields are updated. Returns the mutated containers for convenience.
"""
function unvec!(cp::ControlParameters, vector::AbstractVector)
    idx = 1

    if !isnothing(cp.ub)
        idx = unvec!(cp.ub, vector; idx = idx, return_idx = true)
    end

    if !isnothing(cp.uq)
        idx = unvec!(cp.uq, vector; idx = idx, return_idx = true)
    end

    if !isnothing(cp.m)
        unvec!(cp.m, vector; idx = idx)
    end

    return cp.ub, cp.uq, cp.m
end

function _vectorize_bounds(bnds::Tuple{T,T}, controls::NamedTuple) where {T<:Real}
    lb = Float64[]
    ub = Float64[]

    lower, upper = bnds

    for control in controls
        n = length(vec(control))
        append!(lb, fill(T(lower), n))
        append!(ub, fill(T(upper), n))
    end

    return lb, ub
end

function _vectorize_bounds(bnds::NamedTuple, controls::NamedTuple)
    keys(bnds) == keys(controls) || error("bounds keys must match controls")

    lb = Float64[]
    ub = Float64[]

    for key in eachindex(controls)
        control = controls[key]
        isnothing(control) && continue #if control is nothing (e.g., no source) just skip. 
 
        bnd = bnds[key]
        n = length(vec(control))

        if control isa NamedTuple
            if bnd isa Tuple
                lb_sub, ub_sub = _vectorize_bounds(bnd, control)
                append!(lb, lb_sub); append!(ub, ub_sub)
            elseif bnd isa NamedTuple
                keys(bnd) == keys(control) || error("nested bounds keys for $key must match control vector keys")
                lb_sub, ub_sub = _vectorize_bounds(bnd, control)
                append!(lb, lb_sub); append!(ub, ub_sub)
            end
        else
            lower, upper = bnd
            (lower isa Real && upper isa Real) || error("bounds for $key must be scalars")
            append!(lb, fill(Float64(lower), n)); append!(ub, fill(Float64(upper), n))
        end
    end

    return lb, ub
end

"""
    vectorize_control_parameter_bounds(cp::ControlParameters) -> (lb, ub)

Flatten bound tuples (`ub_bounds`, `uq_bounds`, `m_bounds`) into lower/upper
vectors aligned with `vec(cp)`. Bounds may be a single `(lower, upper)` tuple
applied everywhere, a NamedTuple of such tuples (matching keys), or nested
NamedTuples mirroring nested controls. Missing entries default to unbounded.
"""
function vectorize_control_parameter_bounds(cp::ControlParameters)
    total_len = 0
    total_len += isnothing(cp.ub_length) ? 0 : cp.ub_length
    total_len += isnothing(cp.uq_length) ? 0 : cp.uq_length
    total_len += isnothing(cp.m_length) ? 0 : cp.m_length

    lower_bounds = Vector{Float64}(undef, total_len)
    upper_bounds = Vector{Float64}(undef, total_len)

    idx = 1

    if !isnothing(cp.ub)
        lower_ub, upper_ub = _vectorize_bounds(cp.ub_bounds, cp.ub)
        len = length(lower_ub)
        lower_bounds[idx:idx+len-1] .= lower_ub
        upper_bounds[idx:idx+len-1] .= upper_ub
        idx += len
    end
    if !isnothing(cp.uq)
        lower_uq, upper_uq = _vectorize_bounds(cp.uq_bounds, cp.uq)
        len = length(lower_uq)
        lower_bounds[idx:idx+len-1] .= lower_uq
        upper_bounds[idx:idx+len-1] .= upper_uq
        idx += len
    end
    if !isnothing(cp.m)
        lower_m, upper_m = _vectorize_bounds(cp.m_bounds, cp.m)
        len = length(lower_m)
        lower_bounds[idx:idx+len-1] .= lower_m
        upper_bounds[idx:idx+len-1] .= upper_m
    end

    return lower_bounds, upper_bounds
end
