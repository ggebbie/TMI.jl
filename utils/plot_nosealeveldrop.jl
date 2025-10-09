import Pkg; Pkg.activate(".")

using TMI
using Revise
using LinearAlgebra, Statistics, SparseArrays
using Plots
# using GeoPythonPlot


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

for TMIversion in versions
    suffix = "sealeveldrop"
    println(TMIversion)
    TMIversion_compressed =  TMIversion*"_$suffix.nc"
    TMIfile_compressed = TMI.pkgdatadir("TMI_" * TMIversion_compressed)
    γ = Grid(TMIfile_compressed)
    θ = readfield(TMIfile_compressed,"θ",γ)
    SST = getsurfaceboundary(θ); surface_wet = γ.wet[:, :, surfaceindex(γ)]
    # θ_ref = readfield(TMIfile_compressed,"θ_reference",γ)

    area = cellarea(γ)
    MOT_LGM = mean(θ)
    SST_LGM = sum((SST.tracer .* area.tracer)[surface_wet]) / sum(area.tracer[surface_wet])

    scatter!(p, [CE1990_SST - SST_LGM], [CE1990_MOT - MOT_LGM], 
             label = TMIversion * "_$suffix", markersize = 10)
end
xlabel!(p, "ΔSST")
ylabel!(p, "ΔMOT")

p