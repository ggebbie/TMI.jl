#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using GeoPythonPlot # will load optional extension
using SparseArrays
GeoPythonPlot.pygui(true)
TMIversion = versionlist()[1] # G14 has no remote mass fractions
A, Alu, γ, TMIfile, L, B = config(TMIversion);

vol= cellvolume( γ)

## big control volume to stabilize results
xlim = (125,240); ylim = (30,65); zlim = (2000,3000)
#xlim = (209,211); ylim = (34,36); zlim = (2400,2600)
#xlim = (125,240); ylim = (20,40); zlim = (1500,2500)

# control volume
cube = TMI.incube(xlim,ylim,zlim, γ)
cube_up_bc = TMI.get_above_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_down_bc = TMI.get_below_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_east_bc = TMI.get_east_of_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_west_bc = TMI.get_west_of_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_north_bc = TMI.get_north_of_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_south_bc = TMI.get_south_of_cube_boundary_condition(xlim,ylim,zlim, γ)

# bcube = (;
#         up = cube_up_bc,
#         down = cube_down_bc,
#         north = cube_north_bc,
#         east = cube_east_bc,
#         south = cube_south_bc,
#         west = cube_west_bc
#         #  
#         )

function filter_A_cube_bc(cube_bc, A, γ)
    rows = 1 * selectdim(γ.R, cube_bc.dim, Int(cube_bc.k))
    rows = rows[cube_bc.tracer .> 0.0]
    Abc = deepcopy(A)
    for row in rows
        Abc[row, :] .= 0.0
        Abc[row, row] = 1.0
    end
    return dropzeros(Abc)
end

bcube = (; south = cube_south_bc, surface = zerosurfaceboundary(γ),)
Abc = filter_A_cube_bc(cube_south_bc, A, γ)
Abclu = lu(Abc)

cube_bc_invert = steadyinversion(Abclu,bcube,γ)
#inversion of a 1 boundary condition in ocean interior has 
#entries that are all less than zero, seems problematic
maximum(cube_bc_invert) 
sum(cube_bc_invert)

sum(A .- Abc)

