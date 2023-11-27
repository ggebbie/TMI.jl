#=
Code snippets to
- edit existing regions in region NetCDF file
- rename existing regions
- add new regions
- use a nearest neighbor approximation to generate regions in 90x45 file from 180x90

Don't run this whole script at once, comments should indicate how to run it 
=#
import Pkg; Pkg.activate(".")
using PyPlot, TMI, PyCall, DrWatson
using CSV, DataFrames
ccrs = pyimport("cartopy.crs")
close("all")

datapath = joinpath(dirname(projectdir()), "data/regions")
(!).(isdir(datapath)) && mkdir(datapath)

"""
function BC_to_CSV(b::BoundaryCondition, filepath)

    saves a BoundaryCondition as CSV for editing purposes
    paired with CSV_to_traceer
    converts land values to NaN first, because that makes hand-editing in
    an Excel equivalent easier 

    # Args
    - `b` : BoundaryCondition to write as a file
    - `filepath` : output file name, should be .csv 
"""
#save a boundary condition as a csv
#to make this easier to edit as a CSV, its helpful if the land values are NaN
function BC_to_CSV(b::BoundaryCondition, filepath::String)
    isfile(filepath) && error("file already exists!") 
    wet = 1 .- convert(Matrix{Float64}, b.wet) #turn the "true"s of the wet mask to false
    wet[wet .== 1] .= NaN #turn land into NaN 
    CSV.write(filepath, DataFrame(convert(Matrix{Float64}, b.tracer) .+ wet, :auto)) #add the tracer, and land = NaN mask 
end

"""
function CSV_to_tracer(filepath::String)

    converts a CSV file into a BitMatrix tracer to put into a BoundaryCondition
    sets the NaNs from BC_to_CSV to 0 
"""
function CSV_to_tracer(filepath::String, index = 1)
    tracer = Matrix(CSV.read(filepath, DataFrame))
    tracer[tracer .!= index] .= 0
    tracer[tracer .== index] .= 1
    return convert(BitMatrix, tracer) #return as BitMatrix 
end

"""
function plot_regions(TMIversion::String, TMIregion::String, savepath = nothing)

    generic plot of all regions in a specified file 
"""
function plot_regions(TMIversion::String, TMIregion::String, savepath = nothing)
    figure(figsize = (20,20))
    regionlist = TMI.return_regions(TMIversion, TMIregion) 
    N = length(regionlist)
    #_, _, γ, _, _, _ = config_from_nc(TMIversion)
    for (i, r) in enumerate(regionlist)
        ax = subplot(convert(Int64, ceil(N/5)), 5, i, projection = ccrs.PlateCarree())
        ax.coastlines()
        b = TMI.surfaceregion(TMIversion, r, TMIregion)
        pcolormesh(b.i, b.j, b.tracer', cmap = "cool")
        title(r)
    end
    suptitle(TMIversion)
    tight_layout()
    (!).(isnothing(savepath)) && savefig(savepath)
end

# ====== CODE TO EDIT EXISTING REGIONS ======= #
# write every BoundaryCondition in the default regions file as a CSV file
# NaNs will be land, 0s will be wet, nonregion, and 1s will be wet, region
TMIversionlist = ["modern_90x45x33_GH10_GH12", "modern_180x90x33_GH11_GH12"]
for TMIversion in TMIversionlist 
    for r in TMI.return_regions(TMIversion) 
        b = TMI.surfaceregion(TMIversion, r)
        BC_to_CSV(b, joinpath(datapath, TMIversion * "_" * r * ".csv"))
    end
end

TMIregionlist = joinpath.(dirname(projectdir()), "data", ["regions_90x45_edit", "regions_180x90_edit"] .* ".nc") # we will WRITE to these new files 
# then, edit any points in an spreadsheet editor
# then, write to new file in `TMIregionlist` using the following code 
for (TMIregion, TMIversion) in zip(TMIregionlist, TMIversionlist) 
    TMIfile = download_ncfile(TMIversion)
    γ = Grid(TMIfile)
    for r in TMI.regionlist()
        tracer = CSV_to_tracer(joinpath(datapath, TMIversion * "_" * r * ".csv"))
        b = TMI.surfaceregion(TMIversion, r)
        bnew = BoundaryCondition(tracer, b.i,b.j,b.k,b.dim,b.dimval,b.wet,b.name,b.longname,"none") #where T <: Real
        @show TMIregion
        TMI.write(TMIregion, bnew) 
    end
end

# ====== CODE TO RENAME A REGION ======= #
TMIregionlist_er = joinpath.(dirname(projectdir()), "data", ["regions_90x45_edit_rename", "regions_180x90_edit_rename"] .* ".nc") # we will WRITE to these new files 
[TMI.rename_region(TMIversion, :LAB, (:LIS, "Labrador and Irminger Seas"), output_path, TMIregion) for (output_path, TMIversion, TMIregion) in zip(TMIregionlist_er, TMIversionlist, TMIregionlist)]

# ====== CODE TO MAKE NEW REGIONS ======= # 
# Following snippet writes a CSV with 0s everywhere and NaNs for land
# Can then manually change gridcells to numbers corresponding to regions 
TMIversion, TMIregion = TMIversionlist[2], TMIregionlist[2]
CSVpathname = joinpath(datapath, "northatlantic_180x90.csv")
TMIfile = download_ncfile(TMIversion)
γ = Grid(TMIfile)
CSV.write(CSVpathname, DataFrame(zeros(3,1,γ,:none,"","").tracer, :auto))

# after editing the CSV file, write it to a NetCDF using following code
btemplate = TMI.surfaceregion(TMIversionlist[2], "GLOBAL")
# in CSV, 1 will correspond to LAB, 2 to ICE, etc...
shortnames =  Symbol.(["LAB", "IRM", "ICE", "GRE", "NOR", "CEL"])
longnames = ["Labrador Sea", "Irminger Sea", "Iceland Sea", "Greenland Sea", "Norwegian Sea", "Celtic Sea"]
# add regions to old regions file 
for (i, n) in enumerate(shortnames)
    mat = CSV_to_tracer(CSVpathname, i)
    bnew = BoundaryCondition(mat, btemplate.i, btemplate.j, btemplate.k, btemplate.dim, btemplate.dimval, btemplate.wet, n, longnames[i], "none")
    TMI.write(TMIregionlist_er[2], bnew)
end

# ====== CODE TO CONVERT NEW REGIONS TO 90x45  ======= #
wet = TMI.surfaceregion(TMIversionlist[1], "GLOBAL").wet
for n in shortnames
    b_fine = TMI.surfaceregion(TMIversionlist[2], n, TMIregionlist_er[2])
    b_coarse = TMI.decrease_res(b_fine, wet)
    TMI.write(TMIregionlist_er[1], b_coarse) 
end

# ===== MAKE PLOT OF REGIONS FOR BOTH FILES ===== # 
for (TMIversion, TMIregion) in zip(TMIversionlist, TMIregionlist_er)
    plot_regions(TMIversion, TMIregion, projectdir(TMIversion * ".png"))
end
