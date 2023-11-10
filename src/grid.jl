
"""
    struct Grid

    TMI grid with accounting for wet/dry points
"""
struct Grid{T <: Real}
    lon::Vector{T}
    lat::Vector{T}
    depth::Vector{T}
    wet::BitArray{3}
    interior::BitArray{3}
end

"""
    Base.propertynames(γ::Grid) = (I,R,fieldnames(typeof(x))...)

    Do not store Cartesian and linear indices.
    Compute them on demand.
""" 
Base.propertynames(γ::Grid) = (:I,:R,fieldnames(typeof(γ))...)
function Base.getproperty(γ::Grid, d::Symbol)
    if d === :I
        return cartesianindex(γ.wet)
    elseif d === :R
        return linearindex(γ.wet)
    else
        return getfield(γ, d)
    end
end

"""
    function cartesianindex(wet)
    Read and assemble the grid coordinates
    according to a 3D tracer in x,y,z order
# Arguments
- `wet`: BitArray logical mask for wet points
# Output
- `I`: 3D Cartesian indices
"""
cartesianindex(wet::BitArray{3}) = findall(wet)

"""
    function linearindex(wet)
    Read and assemble the grid coordinates.
# Arguments
- `wet`: 3D mask for wet points
# Output
- `R`: array of linear indices, but not a LinearIndices type
"""
function linearindex(wet)
    R = Array{Int64,3}(undef,size(wet))
    fill!(R,0)
    # R = Array{Union{Int64,Nothing},3}(nothing,size(wet))
    R[wet]=1:sum(wet)
    # R = LinearIndices((1:maximum(it),1:maximum(jt),1:maximum(kt)));
    # R = LinearIndices((it,jt,kt));
    #Rwet = R[γ.wet]
    return R
end

"""
    function checkgrid!(c,wet)

    perform a check of file compatibility
     with grid
"""
function checkgrid!(tracer,wet)
    if sum(isnan.(tracer[wet])) > 0
        println(sum(isnan.(tracer[wet]))," points")
        error("readfield warning: NaN on grid")
    end

    # check for non NaN or nonzero off grid
    # Need to rethink how to do this.
    if sum(isnan.(tracer[.!(wet)])) < length(isnan.(tracer[.!(wet)]))
        println("readfield warning: non-NaN value off grid")
        println("resetting to NaN")
        tracer[.!(wet)] .= NaN
    end
end
