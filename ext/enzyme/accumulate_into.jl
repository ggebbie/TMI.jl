import Enzyme

"""
    Enzyme.accumulate_into(gab_ref::Ref{<:Union{Source,Field,BoundaryCondition,MassFraction}}, accumulated, ge)

Merge one reverse-mode gradient contribution `ge` into a TMI accumulator `gab`.
This handles both array fields used in this codebase:
- `.tracer` for `Source`, `Field`, `BoundaryCondition`
- `.fraction` for `MassFraction`

Only wet-point values are differentiable. Metadata is copied through unchanged.
After merge, `ge` is zeroed so Enzyme will not count it twice.

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
) where {T <: Union{Source, Field, BoundaryCondition, MassFraction}}
    if !haskey(accumulated, gab_ref)
        gab = gab_ref[]

        # Only value arrays are differentiable; axes, masks, names, and units
        # are metadata carried by the TMI object.
        if gab isa MassFraction
            gab.fraction[wet(gab)] .+= ge.fraction[wet(ge)]
        else
            gab.tracer[wet(gab)] .+= ge.tracer[wet(ge)]
        end

        # Enzyme expects moved gradient contributions to be cleared from the
        # increment accumulator after they have been merged.
        if ge isa MassFraction
            ge.fraction[wet(ge)] .= zero(eltype(ge.fraction))
        else
            ge.tracer[wet(ge)] .= zero(eltype(ge.tracer))
        end

        # Preserve Enzyme's alias bookkeeping so the same gradient object is not
        # accumulated twice.
        accumulated[gab_ref] = (gab_ref, ge)
    end
    return accumulated[gab_ref]
end

"""
    Enzyme.accumulate_into(gab_ref::Ref{<:NamedTuple{...,Tuple{Vararg{<:Union{Source,Field,BoundaryCondition,MassFraction}}}}}, accumulated, ge)

Same merge operation as above, but for `NamedTuple` bundles of TMI controls.
Each entry is merged on wet points with the right storage array
(`.tracer` vs `.fraction`).
"""
function Enzyme.accumulate_into(
    gab_ref::Base.RefValue{T},
    accumulated::IdDict,
    ge::T,
) where {names, T <: NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}}}
    if !haskey(accumulated, gab_ref)
        gab = gab_ref[]
        for k in keys(gab)
            gabk = gab[k]
            gek = ge[k]
            if gabk isa MassFraction
                gabk.fraction[wet(gabk)] .+= gek.fraction[wet(gek)]
                gek.fraction[wet(gek)] .= zero(eltype(gek.fraction))
            else
                gabk.tracer[wet(gabk)] .+= gek.tracer[wet(gek)]
                gek.tracer[wet(gek)] .= zero(eltype(gek.tracer))
            end
        end
        accumulated[gab_ref] = (gab_ref, ge)
    end
    return accumulated[gab_ref]
end
