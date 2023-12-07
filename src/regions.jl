"""
     function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    Read an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)
    Return a BoundaryCondition
"""
function surfaceregion(TMIversion::String,region::Union{String,Symbol}, region_filepath = nothing)::BoundaryCondition
    #get wet mask 
    TMIfile = download_ncfile(TMIversion)
    γ = Grid(TMIfile)
    wet = γ.wet[:, :, 1]

    
    file,Nx,Ny = download_regionfile(TMIversion::String, region_filepath) 
    
    # version 1: region masks were stored with "d" prefix
    # new version: regions defined by region name alone
    tracername = Symbol(region)

    #dsfc = ncread(file,tracername) # using NetCDF.jl
    ds = Dataset(file,"r") # using NCDatasets.jl
    v = ds[tracername]
    units = v.attrib["units"]
    longname = v.attrib["longname"]

    lon = ds["lon"][:]
    lat = ds["lat"][:]
    depth = ds["depth"][1]
    mask = Bool.(v[:,:]) # use BitMatrix to save some bits

    #name = v.attrib["name"] # error
    close(ds)
    
    #b = BoundaryCondition(mask,lon,lat,depth,3,1,trues(parse(Int64,Nx),parse(Int64,Ny)),Symbol(region),longname,units)
    b = BoundaryCondition(mask,lon,lat,depth,3,1,wet,Symbol(region),longname,units)
    return b
end

"""
function download_regionfile(TMIversion::String)

Return file name of regional masks.

Also download file from Google Drive, if not already done.

File name refers to the 2D size of domain. It is the same for modern and LGM configs and only depends on number of points in Nx and Ny directions.
"""
function download_regionfile(TMIversion::String, region_filepath = nothing)

    Nx,Ny,Nz = gridsize(TMIversion)
    filename = "regions_"*Nx*"x"*Ny*".nc"
    pathname = pkgdatadir(filename)

    if isnothing(region_filepath) 
        #make pkgdatadir() if it doesn't exist 
        !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 

        #if TMIfile doesn't exist, get GDrive url and download 
        if !isfile(pathname)
            println("read via GoogleDrive.jl")
            #- `url`: Google Drive URL for data
            url = regionurl(filename)
            google_download(url,pkgdatadir())
        end
        return pathname,Nx,Ny
    else
        return region_filepath, Nx, Ny 
    end
end

""" 
    function regionurl(TMIversion)
    placeholder function to give location (URL) of NetCDF Google Drive input
    in the future, consider a struct or Dict that describes all TMI versions.
# Arguments
- `file`: name of file to look for on Google Drive
# Output
- `regionurl`: location (URL) for download of regional mask file
"""
function regionurl(file)
    if file == "regions_90x45.nc"
        return "https://docs.google.com/uc?export=download&id=1JODNUH7KKTT8d80DPSWYHXUHJxIRxoPs"
    elseif file == "regions_180x90.nc"
        return  "https://docs.google.com/uc?export=download&id=1Ymea-b3sxCV9pPuOLPDFIKxY0aUH4BDr"
    else
        return nothing
    end
end

regionlist() =  ("GLOBAL","ANT","SUBANT",
    "NATL","NPAC","TROP","ARC",
    "MED","ROSS","WED","LAB","GIN",
    "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
    "TROPATL","TROPPAC","TROPIND")

regionnames() = Dict("GLOBAL" => "globally uniform",
    "ANT" => "Antarctic",
    "SUBANT" => "Subantarctic",
    "NATL" => "North Atlantic",
    "NPAC" => "North Pacific",
    "TROP" => "tropical and subtropical",
    "ARC" => "Arctic",
    "MED" => "Mediterranean",
    "ROSS" => "Ross Sea sector",
    "WED" => "Weddell Sea sector",
    "LAB" => "Labrador and Irminger Seas",
    "GIN" => "Greenland-Iceland-Norwegian Seas",
    "ADEL" => "Adélie Land sector",
    "SUBANTATL" => "Atlantic-sector Subantarctic",
    "SUBANTPAC" => "Pacific-sector Subantarctic",
    "SUBANTIND" => "Indian-sector Subantarctic",
    "TROPATL" => "tropical and subtropical Atlantic",
    "TROPPAC" => "tropical and subtropical Pacific",
    "TROPIND" => "tropical and subtropical Indian")

"""
function regions2nc(TMIversion,γ)

Read vectors from mat file, translate to 3D,
and save surface field to NetCDF file.

Consider deprecating this function.
"""
function regions_mat2nc(TMIversion,γ)

    # save in a different file.
    filenetcdf = pkgdatadir("regions_"*TMIversion*".nc")
    filemat = pkgdatadir("TMI_"*TMIversion*".mat")

    # region names
    # didn't figure out how to use an ordered dict, instead use a tuple
    list = regionlist()
    regionname = regionnames()
    
    matobj = matopen(filemat)
    if haskey(matobj,"d_all")
        d_all = read(matobj,"d_all")
        close(matobj)
    else
        return
    end

    # a kludge for now
    T = eltype(γ.lon)
    
    # iterate over all regions in d_all
    Izyx = cartesianindex(filemat)
    regions = Dict{String,Array{T,2}}()
    regionatts = Dict{String,Dict{String,String}}()
    
    for rr = 1:size(d_all,2)
        # 3D fields in zyx vector format
        # are changed to 3D xyz format
        d = tracerinit(d_all[:,rr],Izyx,γ.wet)

        # just save the surface 2D field
        push!(regions, list[rr] => d[:,:,1])
        
        push!(regionatts, list[rr] =>
         Dict("longname" => regionname[list[rr]], "units" => "[]"))
    end

    TMIgrids, TMIgridsatts = griddicts(γ)

    # iterate in regions Dictionary to write to NetCDF.
    for (varname,varvals) in regions
        dvarname = "d_"*varname
        nccreate(filenetcdf,dvarname,"lon",γ.lon,TMIgridsatts["lon"],"lat",γ.lat,TMIgridsatts["lat"],atts=regionatts[varname])
        println("write ",dvarname)
        ncwrite(varvals,filenetcdf,dvarname)

    end
end

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
