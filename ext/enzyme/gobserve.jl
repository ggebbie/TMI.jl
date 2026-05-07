"""
    interpweights(wis, γ)

Concrete `interpweights` overload for precomputed interpolation indices
(`wis`). This mirrors the weight assembly used in `TMI.interpweights(loc, γ)`
but skips location-to-index conversion.

Why this exists:
- Enzyme was unstable when differentiating directly through `observe(..., wis, γ)`.
- The custom observe rule tapes wet-point weights and uses them in reverse.
"""
function interpweights(
    wis::Tuple{I, I, I},
    γ::Grid{R, 3},
) where {T <: Real, R <: Real, I <: Interpolations.WeightedAdjIndex{2, T}}
    list = vcat(1:length(γ.lon), 1)
    δ = zeros(γ.wet)
    δwrap = view(δ, list, :, :)
    wi1, wi2, wi3 = wis

    for ii in 1:2, jj in 1:2, kk in 1:2
        δwrap[wi1.istart + ii - 1, wi2.istart + jj - 1, wi3.istart + kk - 1] +=
            wi1.weights[ii] * wi2.weights[jj] * wi3.weights[kk]
    end

    s = sum(filter(!isnan, δ))
    if !iszero(s) && s < 1.0
        δ ./= s
    end
    return δ
end

"""
    augmented_primal(config, observe, Duplicated, c, locs, γ)
    augmented_primal(config, observe, Duplicated, c, wis, γ)

Forward pass for point observations under Enzyme.

Both methods compute the primal via `func.val(...)` and tape a shared
representation: one wet-point weight vector per observation. That makes the
reverse pass identical for `locs` and `wis` call paths.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(observe)},
    ::Type{<:Duplicated},
    c::Duplicated{<:Field{T}},
    locs::Const{<:Vector{<:Tuple}},
    γ::Const{<:Grid},
) where {T <: Real}
    y = func.val(c.val, locs.val, γ.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    wetmask = c.val.γ.wet
    wetweights = needs_shadow(config) ? [interpweights(loc, γ.val)[wetmask] for loc in locs.val] : nothing
    return AugmentedReturn(primal, gy, (gy, wetweights))
end

function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(observe)},
    ::Type{<:Duplicated},
    c::Duplicated{<:Field{T}},
    wis::Const{W},
    γ::Const{<:Grid},
) where {T <: Real, I <: Interpolations.WeightedAdjIndex{2, T}, W <: Vector{Tuple{I, I, I}}}
    y = func.val(c.val, wis.val, γ.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    wetweights = nothing
    if needs_shadow(config)
        wetmask = c.val.γ.wet
        wetweights = Vector{Vector{T}}(undef, length(wis.val))
        @inbounds for ii in eachindex(wis.val)
            wetweights[ii] = interpweights(wis.val[ii], γ.val)[wetmask]
        end
    end
    return AugmentedReturn(primal, gy, (gy, wetweights))
end

"""
    reverse(config, observe, Duplicated, tape, c, locs_or_wis, γ)

Reverse pass for point observations.

`tape` contains `(gy, wetweights)`, where `wetweights[i]` is the local Jacobian
of observation `i` with respect to `c.tracer[c.γ.wet]`. Because both forward
paths tape the same format, reverse accumulation is shared.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(observe)},
    ::Type{<:Duplicated},
    tape,
    c::Duplicated{<:Field},
    locs::Const{<:Vector{<:Tuple}},
    γ::Const{<:Grid},
)
    gy, wetweights = tape
    wetmask = c.val.γ.wet
    for ii in eachindex(gy)
        c.dval.tracer[wetmask] .+= gy[ii] .* wetweights[ii]
    end
    return (nothing, nothing, nothing)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(observe)},
    ::Type{<:Duplicated},
    tape,
    c::Duplicated{<:Field{T}},
    wis::Const{W},
    γ::Const{<:Grid},
) where {T <: Real, I <: Interpolations.WeightedAdjIndex{2, T}, W <: Vector{Tuple{I, I, I}}}
    gy, wetweights = tape
    wetmask = c.val.γ.wet
    for ii in eachindex(gy)
        c.dval.tracer[wetmask] .+= gy[ii] .* wetweights[ii]
    end
    return (nothing, nothing, nothing)
end
