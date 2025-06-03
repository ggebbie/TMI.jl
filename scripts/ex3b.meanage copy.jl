#=
% Diagnose the mean or ideal age.
=#

import Pkg; Pkg.activate(".")

using Revise
using TMI
using Test
# using GeoPythonPlot
using Plots
#TMIversion = "modern_90x45x33_GH10_GH12"
using SparseArrays

function replace_row_with_nan!(A::SparseMatrixCSC{Float64}, row::Int)
    # Find all non-zero elements in the row
    nz_cols = findall(!iszero, A[row, :])
    
    # Replace each non-zero element with NaN
    for col in nz_cols
        A[row, col] = NaN
    end
end


function get_vertical_fluxes_from_L(L_ref, γ)
    np = size(L_ref, 1)
    L_top = fill(0.0, np)
    L_bot = fill(0.0, np)
    
    for r in 1:np
        inds = L_ref[r, :].nzind
        vals = L_ref[r, :].nzval

        (length(inds) <= 1) && continue
        
        # Get all cartesian indices
        carts = [γ.I[i] for i in inds]

        # Find reference point (corresponding to current row)

        ref =γ.I[r]

        # Calculate differences and check for vertical neighbors
        for (i, cart) in enumerate(carts)
            diff = Tuple(cart - ref)
            # Check specifically for vertical differences
            if diff == (0, 0, 1)  # Bottom neighbor
                L_bot[r] = vals[i]
            elseif diff == (0, 0, -1)  # Top neighbor
                L_top[r] = vals[i]
            end
        end
    end
    
    return L_top, L_bot
end


TMIversion = "modern_90x45x33_GH10_GH12"
#modern_90x45x33_G14_v2
A, Alu, γ, TMIfile, L, B = config(TMIversion)

bθ = getsurfaceboundary(readfield(TMIfile, "θ", γ))
θ̃ = steadyinversion(Alu,bθ,γ)

L_sym = (L .+ L') ./ 2
L_asym = L - L_sym

r = 41705
γ.I[L[r, :].nzind]

γ.I[L[r, :].nzind] .- γ.I[r]

L[r, :]
L_sym[r, :]
L_asym[r, :]

γ.I[L[r, :].nzind] .- γ.I[r]
γ.I[L_sym[r, :].nzind] .- γ.I[r]
γ.I[L_asym[r, :].nzind] .- γ.I[r]

θ_L =  L * θ̃ 
θ_sym = L_sym * θ̃ 
θ_asym = L_asym * θ̃ 


LON, LAT = (γ.lon' .* ones(length(γ.lat))), (ones(length(γ.lon))' .* γ.lat)
LAT_mask = (LAT' .> 30) .* γ.wet[:, :, 1]
PAC_mask = TMI.surfaceregion(TMIversion, "NPAC") + TMI.surfaceregion(TMIversion, "TROPPAC")
mask = PAC_mask.tracer .* LAT_mask

contourf(mask')
using Statistics

# contourf(γ.lon, γ.lat, LAT_mask')
function average_value(tracer, mask, γ)
    volume = cellvolume(γ).tracer
    weights = volume .* mask
    num = sum(replace(tracer, NaN =>0).* weights, dims = [1, 2])[:]
    denom = sum( weights, dims = [1, 2])[:]
    return num ./ denom
end

function average_value(tracer, mask)
    weights = mask
    num = sum(replace(tracer, NaN =>0).* weights, dims = [1, 2])[:]
    denom = sum( weights, dims = [1, 2])[:]
    return num ./ denom
end

xidx = 50; yidx = 32
max_z_idx = 20

plot(θ_L.tracer[xidx, yidx, max_z_idx:end], -γ.depth[max_z_idx:end], label = "Full", 
    xlabel = "deg C / year")
plot!(θ_sym.tracer[xidx, yidx, max_z_idx:end], -γ.depth[max_z_idx:end], label = "Symmetric")
plot!(θ_asym.tracer[xidx, yidx, max_z_idx:end], -γ.depth[max_z_idx:end], label = "Anti-Symmetric")

plot(average_value(θ_L.tracer, mask, γ)[max_z_idx:end], -γ.depth[max_z_idx:end], label = "Full", 
    xlabel = "deg C / year")
plot!(average_value(θ_sym.tracer, mask, γ)[max_z_idx:end], -γ.depth[max_z_idx:end], label = "Symmetric")
plot!(average_value(θ_asym.tracer, mask, γ)[max_z_idx:end], -γ.depth[max_z_idx:end], label = "Anti-Symmetric")


θL_directional = L .* vec(θ̃)'

θL_directional_sym = L_sym .* vec(θ̃)' #these "fluxes" need to referenced by the average value of the NPAC
θL_directional_asym = L_asym .* vec(θ̃)'

θL_asym_top, θL_asym_bot = get_vertical_fluxes_from_L(θL_directional_asym, γ)
θL_sym_top, θL_sym_bot = get_vertical_fluxes_from_L(θL_directional_sym, γ)

# θL_top, θL_bot = get_vertical_fluxes_from_L(θL_directional, γ)
# L_top, L_bot = get_vertical_fluxes_from_L(L, γ)

LON, LAT = (γ.lon' .* ones(length(γ.lat))), (ones(length(γ.lon))' .* γ.lat)
LAT_mask = (LAT' .> 30) .* γ.wet[:, :, 27]
PAC_mask = TMI.surfaceregion(TMIversion, "NPAC") + TMI.surfaceregion(TMIversion, "TROPPAC")
mask = PAC_mask.tracer .* LAT_mask .* cellarea(γ).tracer



θL_asym_bot_NPAC = unvec(θ̃, θL_asym_bot).tracer[:, :, 27:28]
θL_sym_bot_NPAC = unvec(θ̃, θL_sym_bot).tracer[:, :, 27:28]

LON, LAT = (γ.lon' .* ones(length(γ.lat))), (ones(length(γ.lon))' .* γ.lat)
LAT_mask = (LAT' .> 30) .* γ.wet[:, :, 28]
PAC_mask = TMI.surfaceregion(TMIversion, "NPAC") + TMI.surfaceregion(TMIversion, "TROPPAC")
mask = PAC_mask.tracer .* LAT_mask .* cellarea(γ).tracer
θL_asym_top_NPAC = unvec(θ̃, θL_asym_top).tracer[:, :, 27:28]
θL_sym_top_NPAC = unvec(θ̃, θL_asym_top).tracer[:, :, 27:28]



unvec(θ̃, θL_asym_top).tracer[50, 32, :][end-5:end]
unvec(θ̃, θL_asym_bot).tracer[50, 32, :][end-6:end-1]

unvec(θ̃, θL_sym_top).tracer[50, 32, :][end-5:end]
unvec(θ̃, θL_sym_bot).tracer[50, 32, :][end-6:end-1]

unvec(θ̃, θL_sym_top).tracer[50, 32, :][end-5:end]

average_value(θL_asym_top_NPAC, mask ) + average_value(θL_asym_bot_NPAC, mask)
 #C/year -> cK/cent
average_value(θL_sym_top_NPAC, mask ) + average_value(θL_sym_bot_NPAC, mask)


θ_L.tracer[50, 32, :][end-5:end]
 
unvec(θ̃,L_bot).tracer[50, 32, :][end-6:end-1]
unvec(θ̃,L_top).tracer[50, 32, :][end-5:end]

unvec(θ̃,θL_bot).tracer[50, 32, :][end-6:end-1] .- unvec(θ̃,θL_top).tracer[50, 32, :][end-5:end]

plot(unvec(θ̃,θL_asym_top).tracer[50, 32, 2:end] .+ unvec(θ̃,θL_asym_top).tracer[50, 32, 1:end-1])