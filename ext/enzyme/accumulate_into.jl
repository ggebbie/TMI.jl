import Enzyme

"""
    Enzyme.accumulate_into(gab_ref::Ref{<:Union{Source,Field,BoundaryCondition}}, accumulated, ge)

Merge a gradient contribution into the accumulator for a TMI struct
(`Source`, `Field`, `BoundaryCondition`) during Enzyme's reverse pass.
For `e = a + b`, `ge = dJ/de` flows back as `ga += ge` and `gb += ge` on
wet tracer entries. Metadata (axes, masks, names, units) is not
differentiable.

A custom rule for `+` doesn't work. In `Const(a) + Duplicated(b)`, Enzyme
sees `a` and `b` as structurally identical and rejects the setup, since
`a.tracer` could be differentiated like `b.tracer` but `Const` says it
should not. `A \\ d` (see `gldiv_field.jl`) avoids this because the solve
goes through UMFPACK / SuiteSparse (no Julia code to trace) and has only
one `Field` operand.

Enzyme already differentiates the elementwise tracer math in `+`. The
only gap is the per-struct merge. This method adds `ge.tracer` into
`gab.tracer` on wet entries, leaves metadata alone, and zeros `ge.tracer`
so it is not counted twice. The same merge serves `-`, scalar `*`, and
broadcasting — Enzyme has already applied the sign or scale to `ge`
upstream.

Arguments:
- `gab_ref`: reference to the gradient accumulator (`ga` or `gb`).
- `accumulated`: Enzyme's record of gradients already merged this pass,
  used to avoid double-counting.
- `ge`: incoming gradient contribution from `e`.
"""
function Enzyme.accumulate_into(
    gab_ref::Base.RefValue{T},
    accumulated::IdDict,
    ge::T,
) where {T <: Union{Source, Field, BoundaryCondition}}
    if !haskey(accumulated, gab_ref)
        gab = gab_ref[]

        # Only tracer values are differentiable; axes, masks, names, and units
        # are metadata carried by the TMI object.
        gab.tracer[wet(gab)] .+= ge.tracer[wet(ge)]

        # Enzyme expects moved gradient contributions to be cleared from the
        # increment accumulator after they have been merged.
        ge.tracer[wet(ge)] .= zero(eltype(ge.tracer))

        # Preserve Enzyme's alias bookkeeping so the same gradient object is not
        # accumulated twice.
        accumulated[gab_ref] = (gab_ref, ge)
    end
    return accumulated[gab_ref]
end

"""
    Enzyme.accumulate_into(gab_ref::Ref{<:NamedTuple{...,Tuple{Vararg{<:Union{Source,Field,BoundaryCondition}}}}}, accumulated, ge)

Merge a gradient contribution into a `NamedTuple` of TMI structs
(`BoundaryCondition`, `Field`, `Source`) by accumulating tracer gradients entry
wise on wet points. This is needed when Enzyme differentiates code paths that
`deepcopy` a named bundle of boundary conditions.
"""
function Enzyme.accumulate_into(
    gab_ref::Base.RefValue{T},
    accumulated::IdDict,
    ge::T,
) where {names, T <: NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}}}
    if !haskey(accumulated, gab_ref)
        gab = gab_ref[]
        for k in keys(gab)
            gabk = gab[k]
            gek = ge[k]
            gabk.tracer[wet(gabk)] .+= gek.tracer[wet(gek)]
            gek.tracer[wet(gek)] .= zero(eltype(gek.tracer))
        end
        accumulated[gab_ref] = (gab_ref, ge)
    end
    return accumulated[gab_ref]
end
