# Setup environment
import Pkg; Pkg.activate(".")

using Revise, TMI, Test, Plots
using Interpolations, Statistics
using LinearAlgebra, NCDatasets
# using LaTeXStrings

# Helper functions
last_non_nan(v) = findlast(!isnan, v) !== nothing ? v[findlast(!isnan, v)] : NaN
last_non_nan(m::Matrix) = last_non_nan.(eachcol(m))
global_ocean_average(x::AbstractArray, γ::Grid) = sum(replace(x, NaN=>0.0) .* cellvolume(γ).tracer) / sum(cellvolume(γ).tracer)
global_surface_average(x::AbstractArray, γ::Grid) = sum(replace(x[:, :, 1], NaN=>0.0) .* cellvolume(γ).tracer[:, :, 1]) / sum(cellvolume(γ).tracer[:, :, 1])

function extract_GH19_year(year::Int)
    # fname = "../GH19.jl/data/Theta_OPT-0015.nc"
    fname = "../GH19.jl/data/Theta_EQ-1750.nc"
    dataset = NCDataset(fname, "r")
    PI_index = findall(dataset["year"][:] .== year)[1]
    PI_theta = permutedims(dataset["theta"][PI_index, :, :, :], [3, 2, 1])
    return PI_theta
end

function bootstrap_PI_lgm_differences(N_sample, Nboot; sampling_method=:uniform, PI_year=1850)
    # Configure LGM dataset
    TMIversion_lgm = "LGM_90x45x33_G14"
    A_lgm, Alu_lgm, γ_lgm, TMIfile_lgm, L_lgm, B_lgm = config(TMIversion_lgm)
    LGM_theta = readfield(TMIfile_lgm, "θ", γ_lgm).tracer
    
    # Configure PI dataset
    TMIversion_PI = "modern_180x90x33_GH11_GH12"
    A_PI, Alu_PI, γ_PI, TMIfile_PI, L_PI, B_PI = config(TMIversion_PI)
    PI_theta_1850 = extract_GH19_year(PI_year)

    # Initialize bootstrap profiles dictionary
    boot_profiles = Dict(
        "PI_surface" => zeros(N_sample, Nboot), "PI_bottom" => zeros(N_sample, Nboot), 
        "LGM_surface" => zeros(N_sample, Nboot), "LGM_bottom" => zeros(N_sample, Nboot)
    )

    for ib in 1:Nboot
        # Generate sampling locations
        locs = Vector{Tuple{Float64, Float64}}(undef, N_sample)
        [locs[i] = wetsurfacelocation(γ_lgm, γ_PI; sampling_method=sampling_method) for i in eachindex(locs)]

        # Sample LGM and PI temperature profiles
        y_lgm, _, _, _ = random_profiles(TMIversion_lgm, "θ", γ_lgm, N_sample; locs=locs)
        boot_profiles["LGM_surface"][:, ib] .= y_lgm[1, :] 
        boot_profiles["LGM_bottom"][:, ib] .= last_non_nan(y_lgm)

        y_PI, _, _, _ = random_profiles(PI_theta_1850, γ_PI, N_sample; locs=locs)
        boot_profiles["PI_surface"][:, ib] .= y_PI[1, :] 
        boot_profiles["PI_bottom"][:, ib] .= last_non_nan(y_PI)
    end

    return boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm
end

function generate_temperature_difference_plot(boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm, 
                                             sampling_method, N_sample, Nboot, output_filename)
    # Calculate average differences
    delta_sst = mean(boot_profiles["PI_surface"] .- boot_profiles["LGM_surface"], dims=1)[:]
    delta_mot = mean(boot_profiles["PI_bottom"] .- boot_profiles["LGM_bottom"], dims=1)[:]

    # Define consistent plot style for all plots
    plot_style = Dict(
        :xlabel => r"\$\\overline{\\Delta{SST}}(^\\circ C)\$", 
        :ylabel => "\$\\overline{\\Delta{MOT}}(^\\circ C)\$",
        :title => "\$\\Delta\$Ocean Temp. (LGM vs. PI)", 
        :markerstrokewidth => 0.1, :markersize => 3,
        :legend => :bottomright, 
        :right_margin => 2Plots.mm, :left_margin => 2Plots.mm,
        :bottom_margin => 2Plots.mm,  # Add extra margin for the legend
        :size => (700, 700), :dpi => 1000,
        :titlefontsize=>18, :guidefontsize=>15,
        :tickfontsize=>13, :legendfontsize=>10,        
        :xlims => (-2, 5), :ylims => (-2, 5),
    )
    
    sampling_description = sampling_method == :uniform ? "uniform grid sampling" : "area-weighted spherical sampling"
    
    # Create the plot
    p = scatter(delta_sst, delta_mot, alpha = 0.6,
        label="Bootstrap avgs, $Nboot dots total, each dot is avg of $N_sample samples\n($sampling_description)"; plot_style...)

    # Add mean and global average points
    highlight = Dict(:markerstrokewidth => 2, :markersize => 10, :order => 2)
    
    scatter!([mean(delta_sst)], [mean(delta_mot)], 
        label="Mean diff averaged over all samples", color=:orangered4; highlight...)
    
    global_sst_diff = global_surface_average(PI_theta_1850, γ_PI) - global_surface_average(LGM_theta, γ_lgm)
    global_mot_diff = global_ocean_average(PI_theta_1850, γ_PI) - global_ocean_average(LGM_theta, γ_lgm)
    
    scatter!([global_sst_diff], [global_mot_diff], 
        label="TMI Volume Weighted Differences\n(EQ-1750 minus G14)", color=:green; highlight...)
    
    savefig(output_filename)
    return p
end

# Run analyses
N_sample = 30
Nboot = 2000
PI_year =1850

#Running bootstrap with uniform sampling
boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm = bootstrap_PI_lgm_differences(
    N_sample, Nboot, sampling_method=:uniform, PI_year=PI_year
)
generate_temperature_difference_plot(
    boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm,
    :uniform, N_sample, Nboot, "LGM_PI_Temp_diff_weighted_uniformly_sampled.png"
)

#"Running bootstrap with spherical sampling
boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm = bootstrap_PI_lgm_differences(
    N_sample, Nboot, sampling_method=:spherical, PI_year=PI_year
)
generate_temperature_difference_plot(
    boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm,
    :spherical, N_sample, Nboot, "LGM_PI_Temp_diff_weighted_spherical_sampled.png"
)