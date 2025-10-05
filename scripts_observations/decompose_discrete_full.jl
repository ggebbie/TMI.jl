# decompose_discrete_full.jl

"""
    decompose_discrete_full(
        f::AbstractArray, g::AbstractArray,
        D_old::AbstractArray{Bool}, D_new::AbstractArray{Bool},
        w_old::AbstractArray, w_new::AbstractArray
    ) -> NamedTuple

Perform a full six‐term discrete decomposition of
```text
    Σ_{i∈D_new} g[i] * w_new[i] − Σ_{i∈D_old} f[i] * w_old[i]
```

Terms:
1. Domain shift:
   Σ_{D_new−D_old} f⋅w_old − Σ_{D_old−D_new} f⋅w_old  
2. Flux change on overlap:
   Σ_{D_old∩D_new} Δf⋅w_old  
3. Weight change on overlap:
   Σ_{D_old∩D_new} f⋅Δw  
4. Cross‐overlap:
   Σ_{D_old∩D_new} Δf⋅Δw  
5. Flux & weight change on new area:
   Σ_{D_new−D_old}(Δf⋅w_old + f⋅Δw + Δf⋅Δw)

# Returns
A named tuple with fields
- domain_shift
- flux_change
- weight_change
- interaction_overlap
- new_domain_interaction
- total_decomposed
- total_direct
- D_old_minus_new
- D_new_minus_old
- D_intersect
"""
function decompose_discrete_full(
    f::AbstractArray, g::AbstractArray,
    D_old::AbstractArray{Bool}, D_new::AbstractArray{Bool},
    w_old::AbstractArray, w_new::AbstractArray
)
    # perturbations
    df = g .- f
    dw = w_new .- w_old

    # sub-domains
    D_old_minus_new = D_old .& .!D_new
    D_new_minus_old = D_new .& .!D_old
    D_intersect     = D_old .& D_new

    # helper to sum over a mask
    integrate(field, mask) = sum(field[mask])

    # (1) Domain shift
    domain_shift = integrate(f .* w_old, D_new_minus_old) -
            integrate(f .* w_old, D_old_minus_new)
    # (2) Flux change on overlap
    integrand_change = integrate(df .* w_old, D_intersect)
    # (3) Weight change on overlap
    weight_change = integrate(f .* dw, D_intersect)
    # (4) Cross-overlap
    overlap_interaction = integrate(df .* dw, D_intersect)
    # (5) Flux & weight change on new area
    new_domain_interaction = integrate(df .* w_old, D_new_minus_old) +
            integrate(f .* dw, D_new_minus_old) +
            integrate(df .* dw, D_new_minus_old)

    total_decomposed = domain_shift + integrand_change + weight_change + overlap_interaction + new_domain_interaction
    total_direct = integrate(g .* w_new, D_new) -
                   integrate(f .* w_old, D_old)

    return (
        domain_shift           = domain_shift,
        integrand_change       = integrand_change,
        weight_change          = weight_change,
        overlap_interaction    = overlap_interaction,
        new_domain_interaction = new_domain_interaction,
        total_decomposed       = total_decomposed,
        total_direct           = total_direct,
        D_old_minus_new        = D_old_minus_new,
        D_new_minus_old        = D_new_minus_old,
        D_intersect            = D_intersect
    )
end
