"""
    augmented_primal(config, adjustsource!, Const, q, u)

Custom reverse-mode primal rule for

    adjustsource!(q::Source, u::Source)

This rule currently supports only the non-logscale path (`q.logscale == false`
and `u.logscale == false`), which is the path used in current enzyme tests and
ex9 workflows.
"""
function augmented_primal(
    ::RevConfigWidth{1},
    func::Const{typeof(adjustsource!)},
    ::Type{<:Const},
    q::Annotation{<:Source},
    u::Annotation{<:Source},
)
    @assert !q.val.logscale "Enzyme adjustsource! rule supports q.logscale=false only"
    @assert !u.val.logscale "Enzyme adjustsource! rule supports u.logscale=false only"
    func.val(q.val, u.val)
    return AugmentedReturn(nothing, nothing, nothing)
end

"""
    reverse(config, adjustsource!, Const, tape, q, u)

For `q += u` on interior points, the local Jacobian wrt `u` is identity. The
incoming gradient on `q` is accumulated into `u`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(adjustsource!)},
    ::Type{<:Const},
    tape,
    q::Annotation{<:Source},
    u::Annotation{<:Source},
)
    if !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        qd = q isa MixedDuplicated ? q.dval[] : q.dval
        ud.tracer[ud.γ.interior] .+= qd.tracer[qd.γ.interior]
    end
    return (nothing, nothing)
end
