"""
    massconservation(inversion)

Equality constraints `Σ m = 1` in every interior cell, written for the scaled
controls `z = uvec ./ σ` used by Ipopt. The constraints are linear in `z`,
`Σ m = P m₀ + jacobian z`, so the Jacobian is constant.

# Arguments
- `inversion::Inversion`: inverse problem

# Output
- the Optimization.jl arguments `cons`, `cons_j`, `cons_jac_prototype`, `lcons`,
  `ucons`; all `nothing`, Optimization.jl's defaults, when mass fractions are
  not adjusted
"""
function massconservation(inversion::Inversion)
    ranges = inversion.ranges
    isempty(ranges.m) && return (cons=nothing, cons_j=nothing, cons_jac_prototype=nothing,
        lcons=nothing, ucons=nothing)
    γ = inversion.y.γ
    R = linearindex(γ.interior)
    rows = reduce(vcat, [R[wet(m)] for m in inversion.m₀])
    P = sparse(rows, 1:length(rows), 1.0, sum(γ.interior), length(rows))
    Pm₀ = P * vec(inversion.m₀)
    jacobian = sparse(rows, collect(ranges.m), inversion.σ[ranges.m], sum(γ.interior),
        length(inversion.σ))
    cons = (Σm, z, _) -> (mul!(Σm, jacobian, z); Σm .+= Pm₀; nothing)
    cons_j = (out, _, _) -> (copyto!(nonzeros(out), nonzeros(jacobian)); nothing)
    return (cons=cons, cons_j=cons_j, cons_jac_prototype=jacobian,
        lcons=ones(size(jacobian, 1)), ucons=ones(size(jacobian, 1)))
end
