
function get_above_cube_boundary_condition(xlim, ylim, zlim, γ::Grid)
    #this is binary boundary condition
    cube_up    = TMI.above_cube(xlim, ylim, zlim, γ)
    A          = Float64.(cube_up)
    cube_flat  = sum(A, dims=(1,2))[:]
    dimval       = findlast(cube_flat .> 0)
    A[γ.wet .== 0] .= NaN
    fld        = Field(A, γ, :nothing, "upper boundary condition", "nothing")
    return TMI.getboundarycondition(fld, 3, dimval, γ)
end

function get_below_cube_boundary_condition(xlim, ylim, zlim, γ::Grid)
    cube_down  = TMI.below_cube(xlim, ylim, zlim, γ)
    A          = Float64.(cube_down)
    cube_flat  = sum(A, dims=(1,2))[:]
    dimval       = findfirst(cube_flat .> 0)
    A[γ.wet .== 0] .= NaN
    fld        = Field(A, γ, :nothing, "lower boundary condition", "nothing")
    return TMI.getboundarycondition(fld, 3, dimval, γ)
end

function get_east_of_cube_boundary_condition(xlim, ylim, zlim, γ::Grid)
    cube_east  = TMI.east_of_cube(xlim, ylim, zlim, γ)
    A          = Float64.(cube_east)
    cube_flat  = sum(A, dims=(2,3))[:]
    dimval     = findfirst(cube_flat .> 0)
    A[γ.wet .== 0] .= NaN
    fld        = Field(A, γ, :nothing, "eastern boundary condition", "nothing")
    return TMI.getboundarycondition(fld, 1, dimval, γ)
end

function get_west_of_cube_boundary_condition(xlim, ylim, zlim, γ::Grid)
    cube_west  = TMI.west_of_cube(xlim, ylim, zlim, γ)
    A          = Float64.(cube_west)
    cube_flat  = sum(A, dims=(2,3))[:]
    dimval     = findlast(cube_flat .> 0)
    A[γ.wet .== 0] .= NaN
    fld        = Field(A, γ, :nothing, "western boundary condition", "nothing")
    return TMI.getboundarycondition(fld, 1, dimval, γ)
end

function get_south_of_cube_boundary_condition(xlim, ylim, zlim, γ::Grid)
    cube_south = TMI.south_of_cube(xlim, ylim, zlim, γ)
    A          = Float64.(cube_south)
    cube_flat  = sum(A, dims=(1,3))[:]
    dimval     = findlast(cube_flat .> 0)
    A[γ.wet .== 0] .= NaN
    fld        = Field(A, γ, :nothing, "southern boundary condition", "nothing")
    return TMI.getboundarycondition(fld, 2, dimval, γ)
end


function get_north_of_cube_boundary_condition(xlim, ylim, zlim, γ::Grid)
    cube_north = TMI.north_of_cube(xlim, ylim, zlim, γ)
    A          = Float64.(cube_north)  
    cube_flat  = sum(A, dims=(1,3))[:]
    dimval     = findfirst(cube_flat .> 0)
    A[γ.wet .== 0] .= NaN
    fld        = Field(A, γ, :nothing, "northern boundary condition", "nothing")
    return TMI.getboundarycondition(fld, 2, dimval, γ)
end