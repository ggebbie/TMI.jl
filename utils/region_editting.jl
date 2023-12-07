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
function return_regions(TMIversion, region_filepath = nothing)

    return list of regions for default or custom region file 
"""
function return_regions(TMIversion::String, region_filepath = nothing)
    Nx,Ny,Nz = gridsize(TMIversion)
    filename = "regions_"*Nx*"x"*Ny*".nc"
    pathname = isnothing(region_filepath) ? pkgdatadir(filename) : region_filepath
    nc = NCDataset(pathname)
    k = keys(nc)
    close(nc) 
    return [i for i in k if i ∉ ["lat", "lon", "depth"]]
end

"""
function remove_region(TMIversion::String, region::Symbol, output_file::String, region_filepath = nothing)
    
    doesn't seem to be any easy way to do this besides write a whole new file without that BC
"""
function remove_region(TMIversion::String, region::Symbol, output_file::String, region_filepath = nothing)
    Nx,Ny,Nz = gridsize(TMIversion)
    filename = "regions_"*Nx*"x"*Ny*".nc"
    pathname = isnothing(region_filepath) ? pkgdatadir(filename) : region_filepath
    nc = NCDataset(pathname, "r")
    key = keys(nc)
    bcs = [surfaceregion(TMIversion, k, region_filepath) for k in key]
    close(nc)
    [TMI.write(output_file, bc) for bc in bcs if bcs.name != region]
end


"""
function rename_region(TMIversion::String, region::Symbol, nameto::Tuple{Symbol, String}, output_file::String, region_filepath = nothing)

    doesn't seem to be any easy way to do this besides write a whole new file without that BC
"""
function rename_region(TMIversion::String, region::Symbol, nameto::Tuple{Symbol, String}, output_file::String, region_filepath = nothing)
    Nx,Ny,Nz = gridsize(TMIversion)
    filename = "regions_"*Nx*"x"*Ny*".nc"
    pathname = isnothing(region_filepath) ? pkgdatadir(filename) : region_filepath
    key = return_regions(TMIversion, region_filepath)
    bcs = [surfaceregion(TMIversion, k, region_filepath) for k in key]
    for bc in bcs
        if bc.name == region
            newbc = BoundaryCondition(bc.tracer, bc.i, bc.j, bc.k, bc.dim, bc.dimval, bc.wet, nameto[1], nameto[2], "none")
            TMI.write(output_file, newbc) 
        else
            TMI.write(output_file, bc)
        end
    end
    
end

"""
function decrease_res(b_fine::BoundaryCondition, wet)

    convert a BoundaryCondition in 180x90 to 90x45 using a nearest neighbor approach 
"""
function decrease_res(b_fine::BoundaryCondition, wet)
    i = StepRange(convert.(Int64, [b_fine.i[1], unique(diff(b_fine.i))[1], b_fine.i[end]])...)
    j = StepRange(convert.(Int64, [b_fine.j[1], unique(diff(b_fine.j))[1], b_fine.j[end]])...)
    itp = Interpolations.scale(interpolate(Matrix(b_fine.tracer), BSpline(Constant())), i, j)
    j_c = collect(-88:4:88)
    i_c = collect(2:4:358)
    return BoundaryCondition(convert(BitMatrix, [itp(x,y) for y in j_c, x in i_c]'), i_c, j_c, 0, b_fine.dim, b_fine.dimval, wet, b_fine.name, b_fine.longname, b_fine.units)
end

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
