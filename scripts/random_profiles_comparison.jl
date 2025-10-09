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
    fname = "../GH19.jl/data/Theta_OPT-0015.nc"
    dataset = NCDataset(fname, "r")
    println(dataset["year"][:])
    PI_index = findall(dataset["year"][:] .== year)[1]
    PI_theta = permutedims(dataset["theta"][PI_index, :, :, :], [3, 2, 1])
    return PI_theta
end


function generate_temperature_difference_plot(boot_profiles, LGM_theta, PI_theta_1850, γ_PI, γ_lgm, 
                                             sampling_method, N_sample, Nboot, output_filename)
    # Calculate average differences
    delta_sst = mean(boot_profiles["PI_surface"] .- boot_profiles["LGM_surface"], dims=1)[:]
    delta_mot = mean(boot_profiles["PI_bottom"] .- boot_profiles["LGM_bottom"], dims=1)[:]

    # Define consistent plot style for all plots
    plot_style = Dict(
        :xlabel => "\$\\overline{\\Delta{SST}}(^\\circ C)\$", 
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
    plot!(collect(plot_style[:xlims]), collect(plot_style[:xlims]), label="1:1",
          linestyle=:dash, linewidth=2, color=:black)

    # Add mean and global average points
    highlight = Dict(:markerstrokewidth => 2, :markersize => 10, :order => 2)
    
    scatter!([mean(delta_sst)], [mean(delta_mot)], 
        label="Mean diff averaged over all samples", color=:orangered4; highlight...)
    
    global_sst_diff = global_surface_average(PI_theta_1850, γ_PI) - global_surface_average(LGM_theta, γ_lgm)
    global_mot_diff = global_ocean_average(PI_theta_1850, γ_PI) - global_ocean_average(LGM_theta, γ_lgm)
    
    scatter!([global_sst_diff], [global_mot_diff], 
        label="TMI Volume Weighted Differences\n(OPT-0015 minus G14)", color=:green; highlight...)

    savefig(output_filename)
    return p
end

# Run analyses
N_sample = 30
Nboot = 100
PI_year =20 #15 CE Equillibriation

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