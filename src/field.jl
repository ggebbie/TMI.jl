"""
    struct Field

    This structure permits the grid to be 
    automatically passed to functions with
    the tracer field.

    This structure assumes the Tracer type to be 
    three-dimensional.

    tracer::Array{T,3}
    γ::Grid
    name::Symbol
    longname::String
    units::String
"""
struct Field{T}
    tracer::Array{T,3}
    γ::Grid
    name::Symbol
    longname::String
    units::String
end

"""
    function Field(tracer::Array{T,3},γ::Grid) where T <: Real

    Outer constructor for Field if there's no worry about
    tracer type, long name, or units.
# Arguments
- `tracer::Array{T,3}`
- `γ::Grid`
# Output
- `field::Field`
"""
#function Field(tracer::Array{T,3},γ::Grid) where T <: Real
#   return Field(tracer,γ,:none,"unknown","unknown")
#end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, x::Field)
    summary(io, x); println(io)
    print(io, "Field size ")
    println(io, size(x.tracer))
    println(io, "Surface view")
    show(io,mime,heatmap(transpose(x.tracer[:,:,1]),zlabel=x.units,title=x.longname))
end

""" 
    function zeros(γ::Grid,name=:none,longname="unknown",units="unknown")::Field

      initialize tracer field on TMI grid
      using a Field struct and constructor
# Arguments
- `γ`::TMI.Grid
# Output
- `d`::Field,  3d tracer field with NaN on dry points
"""
function zeros(γ::Grid{T},name=:none,longname="unknown",units="unknown")::Field where T <: Real

    # preallocate
    tracer = Array{T}(undef,size(γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer[γ.wet] .= zero(T)
    tracer[.!γ.wet] .= zero(T) / zero(T) # NaNs with right type
    return Field(tracer,γ,name,longname,units)
end
