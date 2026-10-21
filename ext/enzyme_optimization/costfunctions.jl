"""
    function costfunction(uvec, inversion, u, F)

Cost function `J = Jdata + Jb + Js + Jq + Jm` of the control vector `uvec`: the
data cost plus the prior costs, from `costterms`. Enzyme differentiates this function.

# Arguments
- `uvec`: control vector
- `inversion::Inversion`: observations, first guesses, and weighting matrices
- `u`: preallocated controls, a copy of `inversion.u`; Enzyme differentiates
  through it
- `F`: UMFPACK factorization refactored in place, see `forwardmodel`
# Output
- `J`: cost function
"""
costfunction(uvec::AbstractVector, inversion::Inversion, u, F) =
    sum(costterms(forwardmodel(uvec, inversion, u, F).c, uvec, inversion))

"""
    function costterms(c, uvec, inversion)

The data cost and the prior costs that add up to the cost function:
- `Jdata = nᵀ Wⁱ n`: data cost, with `n = ỹ - y` the modeled minus the observed
  value at each observation point and `Wⁱ` the inverse covariance of the
  observational errors
- `Jb`, `Jq`, `Jm`: prior costs of the adjustments to the boundary conditions,
  sources, and mass fractions, `Σ (u/σ)²` over their blocks of `uvec`
- `Js`: prior smoothness cost of the boundary adjustments, `uvecᵀ S uvec` with
  `S = ∇²ᵀ Diagonal(L⁴/σ²) ∇²` (`inversion.smoothness`) for each boundary
  condition with a smoothing length `L`

# Arguments
- `c::NamedTuple`: steady-state tracers
- `uvec`: control vector
- `inversion::Inversion`: observations, standard deviations `σ`, and `smoothness`
# Output
- `(; Jdata, Jb, Js, Jq, Jm)`
"""
# @inline: as a separate call, Enzyme's gradient through it allocates 30 MiB more on WOCE
@inline function costterms(c::NamedTuple, uvec::AbstractVector, inversion::Inversion)
    (; y, σ, ranges, smoothness) = inversion
    n = observe(c, y) - vec(y)
    # prior cost of the adjustments in a block of uvec, each with standard deviation σ
    priorcost(block) = sum(i -> (uvec[i] / σ[i])^2, block; init=0.0)
    Jdata = dot(n, y.Wⁱ, n)
    Jb = priorcost(ranges.b)
    Js = dot(uvec, smoothness, uvec)
    Jq = priorcost(ranges.q)
    Jm = priorcost(ranges.m)
    return (; Jdata, Jb, Js, Jq, Jm)
end
