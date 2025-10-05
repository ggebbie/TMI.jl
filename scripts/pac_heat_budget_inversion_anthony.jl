#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using GeoPythonPlot # will load optional extension
using SparseArrays
using Plots
GeoPythonPlot.pygui(true)
TMIversion = versionlist()[1] # G14 has no remote mass fractions
A, Alu, γ, TMIfile, L, B = config(TMIversion);

vol= cellvolume( γ)

# get observations at surface
# set them as surface boundary condition
y = (θ =  readfield(TMIfile, "θ", γ),
#    S = readfield(TMIfile, "Sp", γ),
    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
)

## get steady version of tracers₄
θb = getsurfaceboundary(y.θ)

## preformed phosphate
θ = steadyinversion(Alu,θb,γ)


## big control volume to stabilize results
xlim = (100,255); ylim = (30,65); zlim = (2000,3000)

# control volume
cube = TMI.incube(xlim,ylim,zlim, γ)
cubefield = Field(Float64.(cube), γ, :nothing, "nothing", "nothing")
cubefield.tracer .= cubefield.tracer .* γ.wet
cube_up_bc = TMI.get_above_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_down_bc = TMI.get_below_cube_boundary_condition(xlim,ylim,zlim, γ)
cube_south_bc = TMI.get_south_of_cube_boundary_condition(xlim,ylim,zlim, γ)

#not needed for this problem. 
#cube_north_bc = TMI.get_north_of_cube_boundary_condition(xlim,ylim,zlim, γ)
# cube_east_bc = TMI.get_east_of_cube_boundary_condition(xlim,ylim,zlim, γ)
# cube_west_bc = TMI.get_west_of_cube_boundary_condition(xlim,ylim,zlim, γ)

# single‐BC method
function apply_boundary_conditions_to_A_matrix(
    A::SparseMatrixCSC{T,Int}, bmask::Field) where T
    # pick out the boundary‐rows for this one BC
    Abc = deepcopy(A)
    vbmask = vec(bmask)
    for i in eachindex(vbmask)
        if vbmask[i] > 0.0
            Abc[i,:] .= 0.0
            Abc[i,i] = 1.0
        end
    end
    return dropzeros(Abc)
end


# single‐BC method
function apply_boundary_conditions_to_A_matrix(
    A::SparseMatrixCSC{T,Int}, bc::BoundaryCondition, γ::Grid) where T

    bmask = zeros(γ)
    setboundarycondition!(bmask, bc)

    return apply_boundary_conditions_to_A_matrix(A, bmask)
end

# named‐tuple‐of‐BC method
function apply_boundary_conditions_to_A_matrix(
    A::SparseMatrixCSC{T,Int}, bc_tuple::NamedTuple, γ::Grid) where T
    bmask = zeros(γ)
    setboundarycondition!(bmask, bc_tuple)
    return apply_boundary_conditions_to_A_matrix(A, bmask)
end


function invert_for_interior_bc(A, bcubeA::NamedTuple, bcubeeval::BoundaryCondition)
    Abc = apply_boundary_conditions_to_A_matrix(A, bcubeA, γ)
    Abclu = lu(Abc)
    cube_bc_invert = steadyinversion(Abclu,bcubeeval,γ)
    return cube_bc_invert
end

function invert_for_interior_bc(A, bcmask::Field, bccube::Field)
    Abc = apply_boundary_conditions_to_A_matrix(A, bcmask)
    Abclu = lu(Abc)
    cube_bc_invert = Abclu \ bccube
    return cube_bc_invert
end

function volume_average(field::Field, mask::Field, γ::Grid)
    wet = γ.wet
    num = [0.0]
    denom = [0.0]
    weights = cellvolume(γ).tracer .* mask.tracer
    for ijk in eachindex(wet)
        if (wet[ijk]) & (weights[ijk] > 0)
            num .+= (field.tracer[ijk] * weights[ijk])
            denom .+= weights[ijk]
        end
    end

    return num[1] / denom[1]
end

b_surface = ones(3,1,γ,:bc_surface,"Surface","nondim")


bcA  = (; up = cube_up_bc, down = cube_down_bc, south = cube_south_bc)

bc_inv = map(bc -> invert_for_interior_bc(A, bcA, bc), bcA)

Abc =apply_boundary_conditions_to_A_matrix(A, bcA, γ)
Abc_lu = lu(Abc)

key_list = collect(keys(bc_inv))
labels = uppercasefirst.(String.(key_list))
avg_mass_fraction = [ volume_average(bc_inv[k], cubefield, γ) for k in key_list ]

sum(avg_mass_fraction)

p =  bar(labels, avg_mass_fraction;
    xlabel = "Boundary",
    ylabel = "Mass Fraction in Mid-Depth N. Pac",
    color = 1:length(key_list), 
    legend = false,
    bar_width = 0.5
)

for (lbl, y) in zip(labels, avg_mass_fraction)
    annotate!(
      p,
      lbl,          
      y - 0.04, 
      text(string(round(y, digits=3)), :center, 8, valign=:bottom, 
)
    )
end


setboundarycondition(bc::NamedTuple, γ::Grid) = (bmask = zeros(γ); setboundarycondition!(bmask, bc); return bmask)
setboundarycondition(bc::BoundaryCondition, γ::Grid) = (bmask = zeros(γ); setboundarycondition!(bmask, bc); return bmask)

F₀ = readfield(TMIfile,"F₀",γ)
qPO₄ = readsource(TMIfile,"qPO₄",γ) # use this to define mixed-layer
bmask = setboundarycondition(bcA, γ)

zero_nt(nt::NamedTuple) = (; (k => v*0 for (k,v) in pairs(nt))...)
bc₀ = zero_nt(bcA)
# better to define a Source type
Iq = findall(x -> x > 0,qPO₄.tracer)


qa_bc = zeros(γ); qa_bc.tracer[Iq] = 1 ./ F₀.tracer[Iq]
for i in eachindex(bmask.tracer)
    if bmask.tracer[i] > 0.0
        vec(qa_bc.tracer)[i] = 0.0
    end
end

# zero boundary condition
meanage_NPAC = steadyinversion(Abc_lu,bc₀,γ,q=qa_bc)
volume_average(meanage_NPAC, cubefield, γ)
volume_average(meanage(TMIversion,Alu,γ), cubefield, γ)

avg_mass_fraction = map(bc_i ->  volume_average(bc_i, cubefield, γ), bc_inv)

NPAC_vol = sum(cubefield.tracer .* cellvolume(γ).tracer)
NPAC_age = volume_average(meanage_NPAC, cubefield, γ)

volume_flux_Sv = map(amf_i ->  (1e-6 .* (NPAC_vol .* amf_i)./ NPAC_age) ./ 3.154e+7, avg_mass_fraction)
volume_flux_m3 = map(amf_i ->  ( (NPAC_vol .* amf_i)./ NPAC_age) ./ 3.154e+7, avg_mass_fraction)

volume_average(field::Field, mask::BoundaryCondition, γ::Grid) = volume_average(field, setboundarycondition(mask, γ), γ)
volume_average(field::Field, bcs::NamedTuple, γ::Grid) = map(bc -> volume_average(field, bc, γ), bcs)

bc_temp = volume_average(θ, bcA, γ)
NPAC_temp = volume_average(θ, cubefield, γ)
tmp_anom = (; (k => (bc_temp[k] - NPAC_temp) for (k,v) in pairs(bc_temp))...)

volume_flux_Sv

sum(avg_mass_fraction)

tmp_flux = (; (k => (bc_temp[k] - NPAC_temp) * (volume_flux_m3[k] / NPAC_vol) for (k,v) in pairs(bc_temp))...)
tmp_flux_cK_cent = (; (k => v*3.154e+7* 1e4 for (k,v) in pairs(tmp_flux))...)

volume_average(NPAC_massfrac, bcA, γ)

bcA  = (; up = cube_up_bc, down = cube_down_bc, south = cube_south_bc)

bcA_field = setboundarycondition(bcA, γ)
bcA_field.tracer[bcA_field.tracer .> 0] .= NaN
bcA_field.tracer[bcA_field.tracer .== 0] .= 1
bcA_field.tracer[isnan.(bcA_field.tracer)] .= 0

NPAC_massfrac = invert_for_interior_bc(A, bcA_field, cubefield)

bcA  = (; up = cube_up_bc, down = cube_down_bc, south = cube_south_bc)

bcvol = setboundarycondition(bcA, γ)
bcvol.tracer[bcvol.tracer .> 0.0] .= 1.
volume_average(NPAC_massfrac, bcvol, γ) * 

bcvol