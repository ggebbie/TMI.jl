# Compress LGM solution by rescaling upper 1000m for sea level drop.

import Pkg; Pkg.activate(".")

using TMI, Revise
using LinearAlgebra, Statistics, SparseArrays
using Plots
using DataFrames
using JLD2
using Test

include("compress_utils.jl")

COMPRESSION_DEPTH = 1000.0
COMPRESSION_DEPTH_INT = Int(round(COMPRESSION_DEPTH))
drop_height = 125.0 #select the height that you would like to compress the ocean by

glacial_versions = filter(v->occursin("LGM",v), versionlist())

# Configuration
# TMI_VERSION = glacial_versions[2]
for TMI_VERSION in glacial_versions[1:end]
    println("Reprocessing $(TMI_VERSION)"); sleep(10.0)
    # Load grids and fields
    A_no_drop, _, γ_no_drop, TMIfile_no_drop, L, B = config(TMI_VERSION * "_nosealeveldrop")

    Δz_no_drop = layerthickness(γ_no_drop)
    zf_no_drop = get_z_faces(Δz_no_drop)
    z_centers_no_drop = get_z_centers(zf_no_drop)
    nz_no_drop = length(z_centers_no_drop)

    # Setup uncertainty weights, try to match surface and lower 1000.0 meters 
    kglacial = findfirst(==(drop_height),γ_no_drop.depth)

    compression_boundary_index = Base.findlast(z_centers_no_drop .<= COMPRESSION_DEPTH)

    H_upper = sum(Δz_no_drop[1:compression_boundary_index])
    compression_scaling = (H_upper - drop_height) / H_upper

    Δz_sealevel_dropped = copy(Δz_no_drop)
    Δz_sealevel_dropped[1:compression_boundary_index] .*= compression_scaling

    zf_sealevel_dropped = get_z_faces(Δz_sealevel_dropped)
    zc_sealevel_dropped = get_z_centers(zf_sealevel_dropped)

    ## generate a new wet mask by using modern ocean depth to set the criteria ##
    deptho_nodrop = sum((γ_no_drop.wet .* reshape(Δz_no_drop, 1, 1, nz_no_drop)), dims = 3)[:, :, 1] #calculate the ocean depth on modern grid
    drop_height_mask = (deptho_nodrop .>= drop_height) #find layers that cannot by compressed 
    wet_sealevel_dropped = (drop_height_mask .* γ_no_drop.wet) #generate new wet mask 
        
    ## generate a compressed grid ##
    b_surface_sealevel_dropped = deepcopy(ones(3, 1, γ_no_drop, :bc_surface, "Surface", "nondim"))
    b_surface_sealevel_dropped.tracer .= replace(x -> x == 1 ? 1.0 : NaN, 1.0 .* wet_sealevel_dropped[:, :, 1]) #update the boundary conditions 
    b_surface_sealevel_dropped.wet .= copy(wet_sealevel_dropped[:, :, 1]) #update the wet mask for boundary condition 

    γ_sealevel_dropped = deepcopy(Grid(b_surface_sealevel_dropped, γ_no_drop)) # deepcopy necessary but don't know why
    γ_sealevel_dropped.depth .= copy(zc_sealevel_dropped) #update z centers
    γ_sealevel_dropped.wet .= copy(wet_sealevel_dropped) #update wet mask

    ## Modify A matrix to redistribute mass originating from shelf points that are now masked  ##
    rows_modern = γ_no_drop.R
    rows_compressed_in_modern = rows_modern[γ_sealevel_dropped.wet]
    rows_modern = rows_modern[γ_no_drop.wet]

    A_sealevel_dropped = copy(A_no_drop)[rows_compressed_in_modern, rows_compressed_in_modern] #subset A matrix with ocean points in compressed grid
    rows = γ_sealevel_dropped.R[γ_sealevel_dropped.wet]
    new_land_rows = rows[sum(A_sealevel_dropped[rows, :], dims = 2) .< -1e-14] #determine where mass fractions do not sum to 1

    for landrow in new_land_rows #update mass fractions in A matrix if previously contained a shallow point
        tmprow = copy(A_sealevel_dropped[landrow, :])
        tmprow[landrow] = 0.0 #zero out the diagonal entry
        tmprow ./= sum(tmprow) #renormalize mass fractions
        tmprow[landrow] = -1.0
        A_sealevel_dropped[landrow, :] = tmprow
    end

    @test minimum(γ_sealevel_dropped.R[γ_sealevel_dropped.wet]) == 1
    @test maximum(γ_sealevel_dropped.R) == sum(γ_sealevel_dropped.wet)
    @test isapprox(maximum(sum(A_sealevel_dropped,dims=2)),1.0)
    @test minimum(sum(A_sealevel_dropped,dims=2))> -1e-14

    σepth = zero(γ_sealevel_dropped.depth)
    σepth[γ_sealevel_dropped.depth .> COMPRESSION_DEPTH] .= 0.01 
    σepth[γ_sealevel_dropped.depth .<= COMPRESSION_DEPTH] .= 5
    σepth[1] = 0.01

    θ_no_drop = readfield(TMIfile_no_drop, "θ", γ_no_drop)
    θ_no_drop_remapped = readfield(TMIfile_no_drop,"θ", γ_sealevel_dropped) #read the old temperature field on the new grid

    # Validate surface BCs
    bc_no_drop = getsurfaceboundary(θ_no_drop).tracer
    bc_sealevel_dropped = getsurfaceboundary(θ_no_drop_remapped).tracer
    max_diff = maximum(abs.(filter(!isnan, bc_no_drop .- bc_sealevel_dropped)))
    max_diff < 1e-10 || error("Surface BC difference not zero after remapping: $max_diff")

    LGM_theta_σ = nan_field(θ_no_drop_remapped, γ_sealevel_dropped)
    for ii in γ_sealevel_dropped.I
        LGM_theta_σ.tracer[ii] = σepth[ii[3]]
    end

    W⁻ = (1 / sum(γ_sealevel_dropped.wet)) .* Diagonal(1 ./ LGM_theta_σ.tracer[γ_sealevel_dropped.wet].^2)

    # Optimize surface boundary conditions
    println("Optimizing surface BCs...")
    b = zerosurfaceboundary(γ_sealevel_dropped)
    u = getsurfaceboundary(θ_no_drop_remapped)
    Alu_comp = lu(A_sealevel_dropped)

    out, f, fg, fg! = steadyclimatology(Alu_comp,b,u,
                                        θ_no_drop_remapped,
                                        W⁻,γ_sealevel_dropped; 
                                        iterations = 250, method = :bounded, u_min = -2.4, u_max = 32.0)

    # Generate final field
    θ̃_approx = rename_field(steadyinversion(Alu_comp, unvec(u, out.minimizer), γ_sealevel_dropped), γ_sealevel_dropped, :θ)

    # sum(cellvolume(γ_sealevel_dropped))

    firstwet = findfirst(sum(γ_sealevel_dropped.wet, dims = (1,2))[:, :, :][:] .> 0.)    
    error_df = DataFrame(
        Region = ["Surface", "Surface-$(COMPRESSION_DEPTH_INT)m", "Deep (>$(COMPRESSION_DEPTH_INT)m)"],
        RMSE = [
            compute_rmse(θ_no_drop_remapped, θ̃_approx, γ_sealevel_dropped.depth .== γ_sealevel_dropped.depth[firstwet], γ_sealevel_dropped),
            compute_rmse(θ_no_drop_remapped, θ̃_approx, γ_sealevel_dropped.depth[firstwet] .<= γ_sealevel_dropped.depth .<= 1000, γ_sealevel_dropped),
            compute_rmse(θ_no_drop_remapped, θ̃_approx, γ_sealevel_dropped.depth .>= 1000, γ_sealevel_dropped)
        ],
        MAE = [
            compute_mae(θ_no_drop_remapped, θ̃_approx, γ_sealevel_dropped.depth .== γ_sealevel_dropped.depth[firstwet], γ_sealevel_dropped),
            compute_mae(θ_no_drop_remapped, θ̃_approx, γ_sealevel_dropped.depth[firstwet] .<= γ_sealevel_dropped.depth .<= 1000, γ_sealevel_dropped),
            compute_mae(θ_no_drop_remapped, θ̃_approx, γ_sealevel_dropped.depth .>= 1000, γ_sealevel_dropped)
        ]
    )

    # Mean temperature statistics
    temp_df = DataFrame(
        Case = ["LGM No Drop", "LGM Optimized"],
        Mean_Ocean_Temp = [compute_mean_temp(θ_no_drop),compute_mean_temp(θ̃_approx)],
        Mean_Surface_Temp = [compute_surface_temp(θ_no_drop, γ_no_drop), compute_surface_temp(θ̃_approx, γ_sealevel_dropped) ]
    )

    relabeled_TMI_Version = "$(TMI_VERSION)_upper$(COMPRESSION_DEPTH_INT)mrelabeled"
    filename = "TMI_$(relabeled_TMI_Version).nc"

    # Save DataFrames for later analysis
    jld2_file = TMI.pkgdatadir("postprocessing_analysis_$(relabeled_TMI_Version).jld2")
    jldsave(jld2_file; error_df, temp_df)

    # Write output
    TMIfile_relabeled = TMI.pkgdatadir(filename)
    isfile(TMIfile_relabeled) && rm(TMIfile_relabeled)
    writefield(TMIfile_relabeled, θ̃_approx)
    TMI.grid2nc(relabeled_TMI_Version, γ_sealevel_dropped)
    TMI.watermassmatrix2nc(relabeled_TMI_Version, A_sealevel_dropped)

    println(error_df)
    println(temp_df)
end