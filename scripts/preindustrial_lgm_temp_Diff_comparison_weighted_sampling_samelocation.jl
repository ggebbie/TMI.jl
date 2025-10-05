import Pkg; Pkg.activate(".")
using Revise
using TMI
using Test
using Interpolations
using Statistics
using LinearAlgebra
Pkg.activate("./scripts")
Pkg.add("Plots")
Pkg.add("PyPlot")
using Plots
Pkg.instantiate()
Pkg.activate("/Users/paban23/Library/CloudStorage/OneDrive-UniversityofConnecticut/PDrive/PhD/Classes/Self_Study/Software_Development_Building_Powerful_Applications/Git_Practice/Software_Development/GH19.jl")
#Pkg.precompile() 
using GH19



Bootst = 2000 

function bootstrap_stats_same_locs(TMIversion_lgm, TMIversion_mod,TMIversion_PI, γ_lgm, γ_mod,γ_PI, θtrue, N, Bootst)
        mean_ys_lgm = Float64[]
        std_ys_lgm = Float64[]
        mean_y_os_lgm = Float64[]
        std_y_os_lgm = Float64[]
        
        mean_ys_mod = Float64[]
        std_ys_mod = Float64[]
        mean_y_os_mod = Float64[]
        std_y_os_mod = Float64[]

        mean_ys_gh19 = Float64[]
        mean_y_os_gh19 = Float64[]
    
        for _ in 1:Bootst
            # Generate random locations and weights for this iteration
            latitudes = collect(minimum(γ_lgm.lat):0.1:maximum(γ_lgm.lat))
            lat_weights = TMI.calculate_lat_weights(latitudes)
            locs = [TMI.surface_wetlocation(γ_lgm,lat_weights) for _ in 1:N]
            lons,lats = [locs[i][1] for i in 1:N], [locs[i][2] for i in 1:N]
            locs_ocean = [TMI.deep_wetlocation(γ_lgm,lons[i],lats[i]) for i in 1:N]
            
            wis_lgm = [TMI.interpindex(locs[i], γ_lgm) for i in 1:N]
            wis_lgm_ocean = [TMI.interpindex(locs_ocean[i], γ_lgm) for i in 1:N]

            wis_mod = [TMI.interpindex(locs[i], γ_mod) for i in 1:N]
            wis_mod_ocean = [TMI.interpindex(locs_ocean[i], γ_mod) for i in 1:N]

            wis_gh19 = [TMI.interpindex(locs[i], γ_PI) for i in 1:N]
            wis_gh19_ocean = [TMI.interpindex(locs_ocean[i], γ_PI) for i in 1:N]

            # LGM observations
            y_lgm, _, _, _, _, _ = TMI.surface_synthetic_observations2(TMIversion_lgm, "θ", γ_lgm, N, σ=0.0, locs=locs, wis=wis_lgm)
            y_o_lgm, _, _, _, _, _ = TMI.ocean_synthetic_observations2(TMIversion_lgm, "θ", γ_lgm, N, σ=0.0, locs=locs_ocean, wis=wis_lgm_ocean)
            push!(mean_ys_lgm, mean(y_lgm))
            push!(std_ys_lgm, std(y_lgm))
            push!(mean_y_os_lgm, mean(y_o_lgm))
            push!(std_y_os_lgm, std(y_o_lgm))
            println("Surface y values: $y_lgm")
            println("Bottom y values: $y_o_lgm")
    
            # Modern observations
            y_mod, _, _, _, _, _ = TMI.surface_synthetic_observations2(TMIversion_mod, "θ", γ_mod, N, σ=0.0, locs=locs, wis=wis_mod)
            y_o_mod, _, _, _, _, _ = TMI.ocean_synthetic_observations2(TMIversion_mod, "θ", γ_mod, N, σ=0.0, locs=locs_ocean, wis=wis_mod_ocean)
            push!(mean_ys_mod, mean(y_mod))
            push!(std_ys_mod, std(y_mod))
            push!(mean_y_os_mod, mean(y_o_mod))
            push!(std_y_os_mod, std(y_o_mod))

            # GH19 obs.
            y_gh19, _, _, _, _, _ = TMI.surface_synthetic_observations3(TMIversion_PI, θtrue, γ_PI, N; σ=0.0, locs=locs, wis=wis_gh19)
            y_o_gh19, _, _, _, _, _  = TMI.ocean_synthetic_observations3(TMIversion_PI, θtrue,γ_PI, N; σ=0.0,locs=locs_ocean, wis=wis_gh19_ocean)
            push!(mean_ys_gh19, mean(y_gh19))
            push!(mean_y_os_gh19, mean(y_o_gh19))

        end
    
        return mean_ys_lgm, std_ys_lgm, mean_y_os_lgm, std_y_os_lgm, mean_ys_mod, std_ys_mod, mean_y_os_mod, std_y_os_mod, mean_ys_gh19,mean_y_os_gh19
end
    





function extract_deepest_valid_layer(field)
        lon_size, lat_size, depth_size = size(field)
        deepest_layer = similar(field[:, :, 1])  # Create a similar array for the result
        
        for i in 1:lon_size
            for j in 1:lat_size
                for k in depth_size:-1:1  # Iterate from the deepest to the shallowest
                    if !isnan(field[i, j, k])
                        deepest_layer[i, j] = field[i, j, k]
                        break
                    end
                end
            end
        end
        
        return deepest_layer
    end


    function compute_LGM_to_PI_diff(LGM_value, Modern_value, PI_value)
        Delta_LGM_to_modern = Modern_value - LGM_value
        Delta_LGM_to_PI = Delta_LGM_to_modern - (Modern_value - PI_value)
        return Delta_LGM_to_PI
    end



    TMIversion_lgm = "LGM_90x45x33_G14"
    A_lgm, Alu_lgm, γ_lgm, TMIfile_lgm, L_lgm, B_lgm = config_from_nc(TMIversion_lgm)
    TMIversion_mod = "modern_90x45x33_G14_v2"
    A_mod, Alu_mod, γ_mod, TMIfile_mod, L_mod, B_mod = config_from_nc(TMIversion_mod)


    #Load GH19 data and grid from GH11_GH12
    TMIversion_PI = "modern_180x90x33_GH11_GH12"
    A_PI, Alu_PI, γ_PI, TMIfile_PI, L_PI, B_PI = config_from_nc(TMIversion_PI)
    file_PI = "../../GH19.jl/data/Theta_EQ-1750.nc"
    θtrue, unit, longname = GH19.readtracer_snapshot(file_PI, "theta", 53)
    θtrue = Field(θtrue, γ_PI, :temperature, "C", "Potential Temperature")


    N = 20
    
    # Perform bootstrap resampling on LGM and modern data
    mean_ys_lgm, std_ys_lgm, mean_y_os_lgm, std_y_os_lgm, mean_ys_mod, std_ys_mod, mean_y_os_mod, std_y_os_mod, mean_ys_gh19, mean_y_os_gh19 =
    bootstrap_stats_same_locs(TMIversion_lgm, TMIversion_mod,TMIversion_PI, γ_lgm, γ_mod,γ_PI,θtrue, N, Bootst)
    


#======================================================================================================#
# What is the mean MOT and SST during LGM from random samples of N=20?

# Hottest, coldest and true mean  temp. during LGM on the scatter plot
y_all_lgm, _, _ = surface_synthetic_observations(TMIversion_lgm, "θ", γ_lgm,σ=0.0)
y_o_all_lgm_3d, _, _ = synthetic_observations(TMIversion_lgm, "θ", γ_lgm,σ=0.0)
y_o_all_lgm = extract_deepest_valid_layer(y_o_all_lgm_3d.tracer)

mean_y_all_lgm = round(mean(y_all_lgm, γ_lgm), digits=2)
mean_y_o_all_lgm = round(mean(y_o_all_lgm_3d), digits=2)

n = 100;
y_hot20_lgm = sort(y_all_lgm)[1:n]
y_o_hot20_lgm = sort(vec(y_o_all_lgm))[1:n]

y_cool20_lgm = sort(y_all_lgm,rev=false)[1:n]
y_o_cool20_lgm = sort(vec(y_o_all_lgm),rev=false)[1:n]



###==============================================================================
# What is the mean MOT and SST during modern time from random samples of N=20?

# Hottest, coldest and true mean  temp. during LGM on the scatter plot
y_all_mod, _, _ = surface_synthetic_observations(TMIversion_mod, "θ", γ_mod,σ=0.0)
y_o_all_mod_3d, _, _ = synthetic_observations(TMIversion_mod, "θ", γ_mod,σ=0.0)
y_o_all_mod = extract_deepest_valid_layer(y_o_all_mod_3d.tracer)

mean_y_all_mod = round(mean(y_all_mod, γ_mod), digits=2)
mean_y_o_all_mod = round(mean(y_o_all_mod_3d), digits=2)
n = 100;
y_hot20_mod = sort(y_all_mod)[1:n]
y_o_hot20_mod = sort(vec(y_o_all_mod))[1:n]

y_cool20_mod = sort(y_all_mod,rev=false)[1:n]
y_o_cool20_mod = sort(vec(y_o_all_mod),rev=false)[1:n]


# Compute Delta LGM to PI
delta_ys_lgm_pi = compute_LGM_to_PI_diff(mean_ys_lgm, mean_ys_mod, mean_ys_gh19)
delta_y_os_lgm_pi = compute_LGM_to_PI_diff(mean_y_os_lgm, mean_y_os_mod, mean_y_os_gh19)


scatter((delta_ys_lgm_pi), (delta_y_os_lgm_pi), 
        label="Mean (20 samples)",order=1,xlims=(-1,5), ylims=(-1,5),legend=:topleft, 
        xlabel="\$\\overline{\\Delta_{SST}}(^\\circ C)\$", 
        ylabel=" \$\\overline{\\Delta_{MOT}}(^\\circ C)\$",
        title="\$\\Delta\$Ocean Temp. (LGM vs. PI)", 
        markerstrokewidth=0.3,markersize=3, left_margin=2Plots.mm,
        legendfontsize=8,dpi=1000)


scatter!([mean_y_all_mod]-[mean_y_all_lgm], [mean_y_o_all_mod]-[mean_y_o_all_lgm], 
        label="True mean diff", color=:orangered4,
         markerstrokewidth=2, legendfontsize=8, markersize=10,order=2,
         )






plot!([-1, 5], [-1, 5], label="1:1", linestyle=:dash, linewidth=2, color=:black)

# scatter!([mean_y_all_mod], [mean_y_o_all_mod], 
#          label="True mean Mod. ($mean_y_all_mod (\$\\bar{\\theta}_s\$) and $mean_y_o_all_mod (\$\\bar{\\theta}\$))",
#          marker=(:diamond, 5), markerstrokewidth=2, opacity=0.9,
#         legendfontsize=8)

savefig("LGM_PI_Temp_diff_4_weighted.png",)


