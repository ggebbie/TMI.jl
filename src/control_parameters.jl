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
- `vector`: cached flattened controls (if `cache_vector=true`)
- `ub_length`, `uq_length`, `m_length`: cached segment lengths (if `cache_lengths=true`)
- `mass_fraction_steps`: cached precomputed mass-fraction step targets (if enabled)
- `A`: cached water-mass matrix (if enabled)
- `dub`, `duq`: cached boundary/source perturbations
- `gm`: cached mass-fraction gradient buffer
- `b`, `q`: cached boundary/source fields initialized from priors
"""
struct ControlParameters{UB,UQ,M,U0,Q0,M0,QU,QS,QM,V,UBL,QL,ML,S,AT,G<:Grid,DUB,DUQ,GDUB,GDUQ,GM,BU,BQ}
    ub::UB; uq::UQ; m::M
    u₀::U0; q₀::Q0; m₀::M0
    Qᵤ::QU; Qₛ::QS; Qₘ::QM
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


function ControlParameters(; γ::Grid,
                            ub = nothing, uq = nothing, m = nothing,
                            u₀ = nothing, q₀ = nothing, m₀ = nothing,
                            Qᵤ = nothing, Qₛ = nothing, Qₘ = nothing,
                            cache_vector = false, cache_lengths = true,
                            cache_mass_fraction_steps = true,
                            cache_A = true)
    # Ensure control entries and priors are independent objects.
    for (val, name) in ((ub, "ub"), (uq, "uq"), (m, "m"),
                        (u₀, "u₀"), (q₀, "q₀"), (m₀, "m₀"),
                        (Qᵤ, "Qᵤ"), (Qₛ, "Qₛ"))
        check_shared_references(val, name)
    end

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
    ub_length = cache_lengths && !isnothing(ub) ? length(vec(ub)) : nothing
    uq_length = cache_lengths && !isnothing(uq) ? length(vec(uq)) : nothing
    m_length = cache_lengths && !isnothing(m) ? length(vec(m)) : nothing

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
                             typeof(vector),
                             typeof(ub_length), typeof(uq_length), typeof(m_length),
                             typeof(m_steps), typeof(A_cached),
                             typeof(γ), typeof(dub_cached), typeof(duq_cached),
                             typeof(gdub_cached), typeof(gduq_cached),
                             typeof(gm_cached),
                             typeof(b_cached), typeof(q_cached)}(
        ub, uq, m, u₀, q₀, m₀, Qᵤ, Qₛ, Qₘ, vector, ub_length, uq_length, m_length,
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
