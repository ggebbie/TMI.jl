"""
    struct MassFractions

If 3D, then names known ahead of time and
can use `struct` instead of NamedTuple

# Attributes
    `north::Field{T}`
    `east::Field{T}`
    `south::Field{T}`
    `west::Field{T}`
    `up::Field{T}`
    `down::Field{T}`
"""
struct MassFractions{T <: Real}
    fraction::Array{T,3}
    γ::Grid
    name::Symbol
    longname::String
    units::String
    position::CartesianIndex{3}
    # north::Field{T}
    # east::Field{T}
    # south::Field{T}
    # west::Field{T}
    # up::Field{T}
    # down::Field{T}
end


function MassFractions(A,
    γ::Grid,
    Δ::CartesianIndex;
    wrap=(true, false, false),
    longname = "mass fraction from neighbor")

    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    # allocate masks
    ngrid = (length(γ.lon), length(γ.lat), length(γ.depth))
    wet = falses(ngrid)
    m   = NaN*ones(ngrid)
    
    for I in Iint
        #Istep = I + step
        Istep, inbounds = step_cartesian(I::CartesianIndex, Δ::CartesianIndex, γ::Grid; wrap=(true,false,false))

        if inbounds
            if γ.wet[Istep]
                wet[I] = true
                m[I] = A[R[I],R[Istep]]
            end
        end
    end

    # make a field
    return MassFractions(m,
        Grid(γ.lon,γ.lat,γ.depth,
            wet,
            BitArray{3}(undef,0,0,0)),
        :m_north,
        longname,
        "unitless",
        Δ)
end

"""
function step_cartesian(I, Δ; wrap=(true,false,false))

# Arguments
- `I::CartesianIndex`: starting point
- `Δ::CartesianIndex`: step
- `γ::Grid`: TMI-defined grid
- `wrap::Tuple{Bool}`: does the dimension wrap around?

# Output
- `Istep::CartesianIndex`: new location
- `inbounds::Bool`: inside the domain bounds?
"""
function step_cartesian(I::CartesianIndex,
    Δ::CartesianIndex,
    γ::Grid;
    wrap=(true,false,false))

    #R = copy(γ.R)
    Istep = I + Δ

    # for bounds checking
    Rfull = CartesianIndices(γ.wet)
    Ifirst, Ilast = first(Rfull), last(Rfull)

    # number of steps outside of bounds
    Ihi = max(Ilast,Istep)-Ilast
    Ilo = Ifirst - min(Ifirst,Istep)
        inbounds = iszero(Ihi) && iszero(Ilo)

    if iszero(Ihi) && iszero(Ilo)
        return Istep, true
    else
        # allocate masks
        ngrid = (length(γ.lon),
            length(γ.lat),
            length(γ.depth))

        wrapstep = zeros(Int,length(wrap))
        for idim in eachindex(wrap)
            # check upper bound
            if Ihi[idim]>0 && wrap[idim]
                #Istep[idim] .-= ngrid[idim] # wr
                wrapstep[idim] = -ngrid[idim] # wr
            end

            # check lower bound
            if Ilo[idim]>0 && wrap[idim]
                #Istep[idim] .+= ngrid[idim] # wr
                wrapstep[idim] = ngrid[idim] # wr
            end
        end

        # take another step to wrap around
        Istep += CartesianIndex(Tuple(wrapstep))
        
        Ihi = max(Ilast,Istep)-Ilast
        Ilo = Ifirst - min(Ifirst,Istep)
        inbounds = iszero(Ihi) && iszero(Ilo)
        return Istep, inbounds
    end
end

massfractions_north(A, γ) = MassFractions(A, γ, CartesianIndex(0,1,0);
    wrap=(true, false, false), longname = "mass fraction from northern neighbor")

massfractions_east(A, γ) = MassFractions(A, γ, CartesianIndex(1,0,0);
    wrap=(true, false, false), longname = "mass fraction from eastern neighbor")

massfractions_south(A, γ) = MassFractions(A, γ, CartesianIndex(0,-1,0);
    wrap=(true, false, false), longname = "mass fraction from southern neighbor")
    
massfractions_west(A, γ) = MassFractions(A, γ, CartesianIndex(-1,0,0);
    wrap=(true, false, false), longname = "mass fraction from western neighbor")

massfractions_up(A, γ) = MassFractions(A, γ, CartesianIndex(0,0,1);
    wrap=(true, false, false), longname = "mass fraction from upper neighbor")
    
massfractions_down(A, γ) = MassFractions(A, γ, CartesianIndex(0,0,-1);
    wrap=(true, false, false), longname = "mass fraction from lower neighbor")

function tracer_contribution(c::Field,m::MassFractions)

    # loop over all locations with "m" values
    Im = cartesianindex(m.γ.wet)

    # allocate masks
    ngrid = (length(c.γ.lon), length(c.γ.lat), length(c.γ.depth))
    # tracer contribution: units C*kg
    mc = zeros(ngrid)
    tracer_contribution!(mc,c,m)

    # make a field
    return Field(mc, c.γ,:mc,"tracer contribution",c.units)
end
function tracer_contribution!(mc::Field,c::Field,m::MassFractions)

    # loop over all locations with "m" values
    Im = cartesianindex(m.γ.wet)
    for I in Im
        #Istep = I + step
        Istep, inbounds = step_cartesian(
            I::CartesianIndex,
            m.position,
            c.γ)
        #; wrap=(true,false,false))

        mc[I] = m.fraction[I]*
            (c.tracer[Istep] - c.tracer[I])
    end

    # make a field
    return Field(mc, c.γ,:mc,"tracer contribution",c.units)
end
function tracer_contribution!(mc::Field,c::Field, m::NamedTuple)
    for m1 in m
        tracer_contribution!(mc,c,m1)
    end
end

# """
# function wetmask_massfractions(γ)

# get masks for "m" variables

# # Input
# - `γ::Grid`
# # Output
# - `wet_north`
# - `wet_east`
# - `wet_south`
# - `wet_west`
# - `wet_up`
# - `wet_down`
# """
# function wetmask_massfractions(γ)
#     # for bounds checking
#     Rfull = CartesianIndices(γ.wet)
#     Ifirst, Ilast = first(Rfull), last(Rfull)

#     # loop over interior points (dirichlet bc for surface)
#     R = cartesianindex(γ.interior)

#     # assume 3D
#     step_north = CartesianIndex(0,1,0)
#     step_east  = CartesianIndex(1,0,0)
#     step_up    = CartesianIndex(0,0,-1)

#     # allocate masks
#     nx = length(γ.lon)
#     ny = length(γ.lat)
#     nz = length(γ.depth)
    
#     wet_north = falses(nx, ny, nz)
#     wet_east = falses(nx, ny, nz)
#     wet_south = falses(nx, ny, nz)
#     wet_west = falses(nx, ny, nz)
#     wet_up = falses(nx, ny, nz)
#     wet_down = falses(nx, ny, nz)
    
#     for I in R
#         Inorth = I + step_north
#         if Inorth[2] ≤ Ilast[2] # bounds check
#             wet_north[I] = γ.wet[Inorth]
#         end

#         Ieast = I + step_east
#         if Ieast[1] ≤ Ilast[1] # bounds check
#             wet_east[I] = γ.wet[Ieast]
#         elseif Ieast[1] > Ilast[1]
#             wet_east[I] = γ.wet[Ieast - nx*step_east] # wraparound
#         end

#         Isouth = I - step_north
#         if Isouth[2] ≥ Ifirst[2] # bounds check
#             wet_south[I] = γ.wet[Isouth]
#         end

#         Iwest = I - step_east
#         if Iwest[1] ≥ Ifirst[1] # bounds check
#             wet_west[I] = γ.wet[Iwest]
#         elseif Iwest[1] < Ifirst[1]
#             wet_west[I] = γ.wet[Iwest + nx*step_east] # wraparound
#         end

#         Iup = I + step_up
#         if Iup[3] ≥ Ifirst[3] # bounds check
#             wet_up[I] = γ.wet[Iup]
#         end

#         Idown = I - step_up
#         if Idown[3] ≤ Ilast[3] # bounds check
#             wet_down[I] = γ.wet[Idown]
#         end
#     end
#     return wet_north, wet_east, wet_south, wet_west, wet_up, wet_down
# end

# """
# function massfractions_north(A,γ)

# Arguments
# - `A`::water-mass matrix
# - `γ::Grid`

# Output
# - `m::Field`: mass fraction from northern neighbor
# """
# function massfractions_north(A,γ)

#     # for bounds checking
#     Rfull = CartesianIndices(γ.wet)
#     Ifirst, Ilast = first(Rfull), last(Rfull)

#     # loop over interior points (dirichlet bc for surface)
#     R = copy(γ.R)
#     Iint = cartesianindex(γ.interior)

#     # assume 3D
#     step_north = CartesianIndex(0,1,0)

#     # allocate masks
#     nx = length(γ.lon)
#     ny = length(γ.lat)
#     nz = length(γ.depth)
    
#     wet = falses(nx, ny, nz)
#     m   = NaN*ones(nx, ny, nz)
    
#     for I in Iint
#         Inorth = I + step_north
#         if Inorth[2] ≤ Ilast[2] # bounds check
#             if γ.wet[Inorth]
#             println(I)
#             println(R[I])
#             println(Inorth)
#             println(R[Inorth])
#                 wet[I] = true
#                 m[I] = A[R[I],R[Inorth]]
#             end
#         end
#     end

#     # make a field
#     return m_north = Field(m,
#         Grid(γ.lon,γ.lat,γ.depth,
#             wet,
#             BitArray{3}(undef,0,0,0)),
#         :m_north,
#         "mass fraction from northern neighbor",
#         "unitless")
# end
