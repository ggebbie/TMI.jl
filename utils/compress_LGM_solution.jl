# Compress LGM solution by rescaling upper 1000m for sea level drop.

import Pkg; Pkg.activate(".")

using TMI, Revise
using LinearAlgebra, Statistics, SparseArrays
using Plots
using DataFrames
using JLD2

include("compress_utils.jl")

COMPRESSION_DEPTH = 1000.0
COMPRESSION_DEPTH_INT = Int(round(COMPRESSION_DEPTH))

glacial_versions = filter(v->occursin("LGM",v), versionlist())


# Configuration
for TMI_VERSION in glacial_versions[2:end]
    println("Reprocessing $(TMI_VERSION)"); sleep(10.0)
    # Load grids and fields
    A_no_drop, _, γ_no_drop, TMIfile_no_drop, L, B = config(TMI_VERSION * "_nosealeveldrop")
    A_sealevel_dropped, _, γ_sealevel_dropped, TMIfile_sealevel_dropped, L, B = config(TMI_VERSION) #may need to download manually

    # Setup uncertainty weights, try to match surface and lower 1000.0 meters 
    kglacial = findfirst(==(125.0),γ_no_drop.depth)
    σepth = zero(γ_sealevel_dropped.depth)
    σepth[γ_sealevel_dropped.depth .> COMPRESSION_DEPTH] .= 0.01 
    σepth[γ_sealevel_dropped.depth .<= COMPRESSION_DEPTH] .= 5
    σepth[kglacial] = 0.01
    σepth[1:(kglacial-1)] .= NaN

    #using this to set surface to 2.5 meters instead of 0 meters, doesn't really change anything though
    z_centers_no_drop = get_z_centers(get_z_faces(layerthickness(γ_no_drop)))
    θ_sealevel_dropped = readfield(TMIfile_sealevel_dropped, "θ", γ_sealevel_dropped)
    θ_no_drop = readfield(TMIfile_no_drop, "θ", γ_no_drop)
    θ_no_drop_remapped = nan_field(θ_no_drop, γ_sealevel_dropped)

    # Interpolate modern field onto compressed grid
    imax, jmax, kmax = size(θ_no_drop.tracer)

    for idx = 1:imax, jdx = 1:jmax
        ref_col = θ_no_drop.tracer[idx, jdx, :]
        ref_wet = .!isnan.(ref_col)
        any(ref_wet) || continue

        comp_col = θ_sealevel_dropped.tracer[idx, jdx, :]
        comp_wet = .!isnan.(comp_col)
        sum(comp_wet) > 0 || continue

        z_comp = γ_sealevel_dropped.depth[findall(comp_wet)]
        upper_comp = z_comp .<= COMPRESSION_DEPTH
        lower_comp = z_comp .> COMPRESSION_DEPTH
        upper_mod = z_centers_no_drop .<= COMPRESSION_DEPTH
        lower_mod = z_centers_no_drop .> COMPRESSION_DEPTH

        new_values = similar(z_comp)

        # Upper ocean: rescale depths
        if sum(upper_comp) > 0 && sum(upper_mod .& ref_wet) > 0
            z_rescaled = rescale_vertical(z_centers_no_drop[upper_mod], extrema(z_comp[upper_comp]))
            new_values[upper_comp] = interpolate_with_nans(
                z_rescaled[ref_wet[upper_mod]], ref_col[upper_mod][ref_wet[upper_mod]],
                z_comp[upper_comp], extrapolate_val=Flat()
            )
        end

        # Lower ocean: no rescaling
        if sum(lower_comp) > 0 && sum(lower_mod .& ref_wet) > 0
            new_values[lower_comp] = interpolate_with_nans(
                z_centers_no_drop[lower_mod][ref_wet[lower_mod]], ref_col[lower_mod][ref_wet[lower_mod]],
                z_comp[lower_comp], extrapolate_val=Flat()
            )
        end

        new_col = fill(NaN, kmax)
        new_col[comp_wet] = new_values
        sum(isnan.(new_col[comp_wet])) == 0 || error("NaNs at ($idx, $jdx)")
        θ_no_drop_remapped.tracer[idx, jdx, :] = new_col
    end

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

    out, f, fg, fg! = steadyclimatology(Alu_comp,b,u,θ_no_drop_remapped,
                                        W⁻,γ_sealevel_dropped; iterations = 1000, method = :bounded, u_min = -2.4, u_max = 30.0)


    # Generate final field
    θ̃_approx = rename_field(steadyinversion(Alu_comp, unvec(u, out.minimizer), γ_sealevel_dropped), γ_sealevel_dropped, :θ)

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

    compressed_TMI_Version = "$(TMI_VERSION)_upper$(COMPRESSION_DEPTH_INT)mcompressed"
    filename = "TMI_$(compressed_TMI_Version).nc"

    # Save DataFrames for later analysis
    jld2_file = TMI.pkgdatadir("postprocessing_analysis_$(compressed_TMI_Version).jld2")
    jldsave(jld2_file; error_df, temp_df)

    # Write output
    TMIfile_compressed = TMI.pkgdatadir(filename)
    isfile(TMIfile_compressed) && rm(TMIfile_compressed)
    writefield(TMIfile_compressed, θ̃_approx)
    TMI.grid2nc(compressed_TMI_Version, γ_sealevel_dropped)
    TMI.watermassmatrix2nc(compressed_TMI_Version, A_sealevel_dropped)

    println(error_df)
    println(temp_df)
end