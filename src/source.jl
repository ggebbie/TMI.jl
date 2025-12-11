"""
    struct Source

    This structure describes Sources, which are
    similar to Fields, but they may be
    1) non-negative
    2) have only interior mask
"""
struct Source{T <: Real, A <: Real, N, F<:AbstractArray{T,N}}
    tracer::F
    γ::Grid{A,N}
    name::Symbol
    longname::String
    units::String
    logscale::Bool
end

function readsource(file,tracername,γ::Grid;logscale=false) 
    # The mode "r" stands for read-only. The mode "r" is the default mode and the parameter can be omitted.
    tracer, units, longname = _read3d(file,tracername)
    checkgrid!(tracer,γ.interior)
    if logscale
        ct = log.(tracer)
    else
        ct = tracer
    end

    return Source(ct,γ,Symbol(tracername),longname,units,logscale)
end
function readsource(matfile,matsourcename,γ::Grid,Izyx)
    # read MATLAB field and transfer zyx format to xyz

    matobj = matopen(matfile)
    varnames, xvarnames = matvarnames(matfile)

    if matsourcename in varnames
        tvar = read(matobj,matsourcename)
    elseif matsourcename in xvarnames
        tvar = read(matobj,"x")[matsourcename]
    end

    # put zyx vector into xyz 3D array
    source = sourceinit(tvar, Izyx, γ)

    # perform a check of file compatibility with grid
    if sum(isnan.(source[γ.interior])) > 0
        error("readsource warning: NaN on interior grid")
    end
    # check for non NaN or nonzero off grid
    if sum(isnan.(source[.!(γ.interior)])) < length(isnan.(source[.!(γ.interior)]))
        println("readsource warning: non-NaN value off grid")
        println("resetting to NaN")
        source[.!(γ.interior)] .= NaN
    end
    ncsourcename = mat2ncsource()[matsourcename]
    units = fieldsatts()[ncsourcename]["units"]
    longname = fieldsatts()[ncsourcename]["longname"]
    logscale = false # not implemented for true case
    close(matobj)
    return Source(source,γ,Symbol(ncsourcename),longname,units,logscale)
end
readmatsource(file,matsourcename,γ::Grid,Izyx = cartesianindex(file)) = readsource(file,matsourcename,γ,Izyx)

writesource(file,q::Source) = write(file,q) 


""" 
    function setsource!(d::Field,q::Field,r::Number)
    apply interior source q to the equation constraints d
# Arguments
- `d`::Field, equation constraints (i.e., right hand side)
- `q`::Field, interior source
- `r`::Number, default = 1.0, stoichiometric ratio
"""
function setsource!(d::Field{T},q::Field{T},r=1.0) where T<: Real
    d.tracer[d.γ.wet] -= r * q.tracer[q.γ.wet]
end
function setsource!(d::Field{T},q::Source{T},r=1.0) where T<: Real
    d.tracer[d.γ.interior] -= r * q.tracer[q.γ.interior]
end

"""
    setsource!(d::Source, q::Source; r=1.0)

Apply one `Source` onto another in-place, scaling by `r` and respecting the
interior mask.
"""
function setsource!(d::Source{T}, q::Source{T}; r=1.0) where T<: Real
    d.tracer[d.γ.interior] .-= r .* q.tracer[q.γ.interior]
end

"""
    setsource!(d::NamedTuple, q::NamedTuple; r=1.0)

Apply a collection of sources to a collection of fields/sources with matching
keys. Skips entries where `q[key]` is `nothing`.
"""
function setsource!(d::NamedTuple, q::NamedTuple; r=1.0)
    for key in keys(q)
        if haskey(d, key) && isnothing(q[key])
            setsource!(d[key], q[key]; r)
        end
    end
end

"""
    replacesource!(d::Source, q::Source)

Overwrite the tracer values of `d` with those of `q` (in-place).
"""
function replacesource!(d::Source{T}, q::Source{T}) where T<: Real
    d.tracer .= q.tracer
end

"""
    replacesource!(d::NamedTuple, q::NamedTuple)

Overwrite a collection of sources/fields in `d` with matching entries from `q`.
"""
function replacesource!(d::NamedTuple, q::NamedTuple)
    for key in keys(q)
        if haskey(d, key) && !isnothing(q[key])
            replacesource!(d[key], q[key])
        end
    end
end


"""
    setsource!(d::NamedTuple, q::NamedTuple, r::NamedTuple)

Apply a collection of sources with per-key ratios. Requires `keys(q) == keys(r)`;
only keys present in `d` are updated.
"""
function setsource!(d::NamedTuple, q::NamedTuple, r::NamedTuple)
    keys(q) == keys(r) || error("setsource!: key mismatch between q and r")
    for key in keys(q)
        haskey(d, key) && setsource!(d[key], q[key], r[key])
    end
    return d
end

"""
    function gsetsource!(gq::Field{T},gd::Field{T},r=1.0)

    Adjoint to `setsource!`
"""
function gsetsource!(gq::Field{T},gd::Field{T},r) where T<: Real
    gq.tracer[gq.γ.wet] -= r * gd.tracer[gd.γ.wet]
end
function gsetsource!(gq::Source{T},gd::Field{T},r) where T<: Real
    gq.tracer[gq.γ.interior] -= r * gd.tracer[gd.γ.interior]
end
"""
    gsetsource(gd, q, r)

Allocate and return the adjoint of `setsource!`, producing a gradient object of
the same type as `q` with values `-r * gd` on the appropriate mask.
"""
function gsetsource(gd::Field{T},q::Union{Field,Source},r) where T<: Real
    gq = 0.0 * q
    gsetsource!(gq,gd,r)
    return gq
end

function adjustsource(q₀::Union{Source,Field,NamedTuple},u::Union{Source,Field,NamedTuple})
    q = deepcopy(q₀)
    adjustsource!(q,u)
    return q
end
    
function adjustsource!(q::Field,u::Field; r::Real = 1.0)
    # write it out so b changes when returned
    q.tracer[q.γ.wet] += r .* u.tracer[u.γ.wet] 
end
function adjustsource!(q::Source,u::Source; r::Real = 1.0)
    if q.logscale && u.logscale
        q.tracer[q.γ.interior] += u.tracer[u.γ.interior]
        q.tracer[q.γ.interior] = exp.(q.tracer[q.γ.interior])
        q.logscale = false
    elseif ~q.logscale && u.logscale
        q.tracer[q.γ.interior] = log.(q.tracer[q.γ.interior])
        q.tracer[q.γ.interior] += u.tracer[u.γ.interior]
        q.tracer[q.γ.interior] = exp.(q.tracer[q.γ.interior])
    elseif q.logscale && ~u.logscale
        error("not implemented: logscale would not lead to non-negative source")
        # u.tracer[u.γ.interior] = log.(u.tracer[u.γ.interior])
        # q.tracer[q.γ.interior] += u.tracer[u.γ.interior]
        # q.tracer[q.γ.interior] = exp.(q.tracer[q.γ.interior])
    else 
        q.tracer[q.γ.interior] += r .* u.tracer[u.γ.interior]
    end        
end
function adjustsource!(q::Source{T}, u::AbstractVector{T};
                       idx::Int = 1,
                       return_idx::Bool = false,
                       r::Real = 1.0) where {T<:Real}
    mask = q.γ.interior
    data = q.tracer
    @inbounds for I in eachindex(mask)
        if mask[I]
            data[I] += r .* u[idx]
            idx += 1
        end
    end
    return return_idx ? idx : nothing
end
function adjustsource!(q::NamedTuple,
                       uvec::AbstractVector;
                       idx::Int = 1,
                       return_idx::Bool = false,
                       r::Real = 1.0)
    for v in q
        if !isnothing(v)
            idx = adjustsource!(v, uvec; idx = idx, return_idx = true, r = r)
        end
    end
    return return_idx ? idx : nothing
end
function adjustsource!(q::NamedTuple,u::NamedTuple; r::Real = 1.0) #where {N1, N2, T <: Real}
    for qkey in keys(q)
        if haskey(u,qkey)
            adjustsource!(q[qkey],u[qkey]; r = r)
        else
            error("adjustsource!: u doesn't have qkey ",qkey)
        end
    end
end
function adjustsource!(q::Union{Field,Source},u::NamedTuple; r::Real = 1.0) #where {N1, N2, T <: Real}
    qkey = :source
    if haskey(u,qkey)
        adjustsource!(q,u[qkey]; r = r)
    else
        error("adjustsource!: u doesn't have source info")
    end
end

"""
    gadjustsource!(gu::Union{Field,Source,NamedTuple},
                   gq::Union{Field,Source,NamedTuple},
                   q₀::Union{Source,NamedTuple};
                   r::Real = 1.0)

Accumulate source gradients in-place. Adds `r * gq` into `gu`, applying the
correct domain mask (wet or interior). For `Source` inputs with `logscale`
the prior `q₀` is used to linearize the log transform before accumulation.
`NamedTuple` methods dispatch over tracers.
"""
function gadjustsource!(gu::Field,gq::Field; r::Real = 1.0)
    # write it out so b changes when returned
    gu.tracer[gu.γ.wet] += r .* gq.tracer[gq.γ.wet] 
end

"""
    gadjustsource!(gu::Source, gq::Source, q₀::Source; r::Real = 1.0)

Accumulate gradients for a log-scaled `Source`. Uses the prior `q₀` to
linearize the exponential transform when `gu.logscale` is true; otherwise
adds `r * gq` over the interior mask.
"""
function gadjustsource!(gu::Source,gq::Source,q₀::Source; r::Real = 1.0)
    q = deepcopy(q₀)
    if gq.logscale && gu.logscale
        error("not implemented")
        # gu.tracer[gu.γ.interior] += gq.tracer[gq.γ.interior]
        # q.tracer[q.γ.interior] = exp.(q.tracer[q.γ.interior])
        # q.logscale = false
        gq.logscale = false
        
    elseif gu.logscale
        q.tracer[q.γ.interior] = log.(q.tracer[q.γ.interior])
        gq.tracer[gq.γ.interior] = gq.tracer[gq.γ.interior] .* exp.(q.tracer[q.γ.interior])
        gu.tracer[gu.γ.interior] += gq.tracer[gq.γ.interior]
        # no need to adjoint this next line?
    elseif gq.logscale
        error("not implemented: logscale would not lead to non-negative source")
        # u.tracer[u.γ.interior] = log.(u.tracer[u.γ.interior])
        # q.tracer[q.γ.interior] += u.tracer[u.γ.interior]
        # q.tracer[q.γ.interior] = exp.(q.tracer[q.γ.interior])
    else 
        gu.tracer[gu.γ.interior] += r .* gq.tracer[gq.γ.interior]
    end        
end

"""
    gadjustsource!(gu::NamedTuple, gq::NamedTuple, q::NamedTuple; r::Real = 1.0)

Loop over tracers in a NamedTuple, accumulating `r * gq[key]` into
`gu[key]` with the appropriate per-tracer handling (logscale or field).
"""
function gadjustsource!(gu::NamedTuple,gq::T,q::T; r::Real = 1.0) where T <: NamedTuple 
    for qkey in keys(gq)
        if haskey(gu,qkey)
            gadjustsource!(gu[qkey],gq[qkey],q[qkey]; r = r)
        end
    end
end

"""
    gadjustsource!(gu::NamedTuple, gq::Union{Source,Field}, q::Union{Source,Field}; r::Real = 1.0)

Pick out the `:source` entry from a gradient NamedTuple and accumulate the
provided gradient there (used when only a single tracer is present).

Used in sourcemap. 
"""
function gadjustsource!(gu::NamedTuple,gq::T,q::T; r::Real = 1.0) where T <: Union{Source,Field}
    qkey = :source
    if haskey(gu,qkey)
        gadjustsource!(gu[qkey],gq,q; r = r)
    end
end
