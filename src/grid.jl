
"""
    struct Grid

    TMI grid with accounting for wet/dry points
"""
struct Grid
    lon::Vector{Float64}
    lat::Vector{Float64}
    depth::Vector{Float64}
    wet::BitArray{3}
    interior::BitArray{3}
    #    I::Vector{CartesianIndex{3}} # index
    #    R::Array{Int,3}
    #    R::LinearIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}} 
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
