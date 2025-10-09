import Pkg; Pkg.activate(".")

using TMI
using Revise
using LinearAlgebra, Statistics, SparseArrays
using Plots
# using GeoPythonPlot

function get_errors(θ̃_approx, θ̃_true, γ_compressed, depth_mask)
    rms(x, y) = sqrt(mean((x .- y).^2))
    abserr(x, y) = mean(abs.(x .- y))
    wet_points = γ_compressed.wet .* reshape(depth_mask, 1, 1, 33)
    model = θ̃_approx.tracer[wet_points]
    data = θ̃_true.tracer[wet_points]

    return (rms_err = rms(model, data), abs_err = abserr(model, data))
    
end

versions = ["LGM_90x45x33_G14", "LGM_90x45x33_OG18", "LGM_90x45x33_GPLS2", "LGM_90x45x33_G14A"]

CE15_MOT = 3.222263849739807
CE15_SST = 16.327818191261613

CE1990_MOT = 3.4289856231485647
CE1990_SST = 16.76771356694541

p = plot()

x = -2:0.1:3; y = fill(2.27, length(x))   
plot!(p, x, x, c = :black, label = nothing, linewidth = 3, linestyle = :dash)

y_upper = y .+ 0.46; y_lower = y .- 0.46
plot!(p, x, y_upper, fillrange=y_lower, fillalpha=0.3, label=nothing, color=:pink)
plot!(p, x, y, color=:pink, label="Seltzer 2024")

err_dicts = Dict()
for TMIversion in versions
    println(TMIversion)
    TMIversion_compressed =  TMIversion*"_compressed.nc"
    TMIfile_compressed = TMI.pkgdatadir("TMI_" * TMIversion_compressed)
    γ = Grid(TMIfile_compressed)
    θ = readfield(TMIfile_compressed,"θ",γ)
    SST = getsurfaceboundary(θ); surface_wet = γ.wet[:, :, 1]
    θ_ref = readfield(TMIfile_compressed,"θ_reference",γ)

    area = cellarea(γ)
    MOT_LGM = mean(θ)
    SST_LGM = sum((SST.tracer .* area.tracer)[surface_wet]) / sum(area.tracer[surface_wet])

    scatter!(p, [CE15_SST - SST_LGM], [CE15_MOT - MOT_LGM], 
             label = TMIversion * "_compressed", markersize = 10)

    compressed_depths = γ.depth
    err_dict = Dict()
    err_dict["surfacetobottom"] = get_errors(θ, θ_ref, γ, 3 .< compressed_depths)

    err_dict["surface"] = get_errors(θ, θ_ref, γ, compressed_depths .< 3)
    err_dict["0to1000"] = get_errors(θ, θ_ref, γ, 3 .<  compressed_depths .<= 1000)
    err_dict["1000tobottom"] = get_errors(θ, θ_ref, γ, 1000 .<  compressed_depths)
    err_dicts[TMIversion] = err_dict
end
xlabel!(p, "ΔSST")
ylabel!(p, "ΔMOT")

p

err_dicts["LGM_90x45x33_G14"]