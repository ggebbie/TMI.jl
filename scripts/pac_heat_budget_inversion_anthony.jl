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
xlim = (100,255); ylim = (30,65); zlim = (2000,3000)

# control volume
cube = TMI.incube(xlim,ylim,zlim, γ)
cube_up_bc = TMI.get_above_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_down_bc = TMI.get_below_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_east_bc = TMI.get_east_of_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_west_bc = TMI.get_west_of_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_north_bc = TMI.get_north_of_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_south_bc = TMI.get_south_of_cube_boundary_condition(xlim,ylim,zlim, γ)

# single‐BC method
function apply_boundary_conditions_to_A_matrix(
    A::SparseMatrixCSC{T,Int}, bc::BoundaryCondition, γ::Grid) where T
    # pick out the boundary‐rows for this one BC
    rows = selectdim(γ.R, bc.dim, bc.dimval)
    rows = rows[rows .!= 0]
    # copy & modify
    A_bc = deepcopy(A)
    for row in rows
        A_bc[row, :] .= 0.0
        A_bc[row, row] = one(T)
    end

    return dropzeros(A_bc)
end

# named‐tuple‐of‐BC method
function apply_boundary_conditions_to_A_matrix(
    A::SparseMatrixCSC{T,Int}, bc_tuple::NamedTuple, γ::Grid) where T
    A_bc = deepcopy(A)
    # apply each BC in turn on the *same* A_bc
    for bc in values(bc_tuple)
        #update the A_bc matrixs
        A_bc .= apply_boundary_conditions_to_A_matrix(A_bc, bc, γ) 
    end
    return A_bc
end


function invert_for_interior_bc(bcube)
    Abc = apply_boundary_conditions_to_A_matrix(A, bcube, γ)
    Abclu = lu(Abc)
    cube_bc_invert = steadyinversion(Abclu,bcube,γ)
    return cube_bc_invert
end

bcube = (; bc = cube_up_bc)
cube_bc_invert = invert_for_interior_bc(bcube)

bcube = (; bc = cube_south_bc)
cube_bc_invert2 = invert_for_interior_bc(bcube)

println("Maximum Value of Inversion: ", maximum(cube_bc_invert))
println("Minimum Value of Inversion: ", minimum(cube_bc_invert))
println("Sum of Inversion: ", sum(cube_bc_invert))

using NaNStatistics

sum(cube_bc_invert.tracer[γ.wet])

using Plots
using NaNStatistics

γ.wet[:, :, 27]

contourf(γ.lon, γ.lat, nansum(cube_bc_invert2.tracer, dims = 3)[:, :, 1]')

contourf(γ.lon, γ.lat, γ.wet[:, :, 27]')
contourf( cube_bc_invert.tracer[:, :, 27]')
