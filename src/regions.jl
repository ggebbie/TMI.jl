"""
     function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    Read an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)
    Return a BoundaryCondition
"""
function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    file = pkgdatadir("regions_"*TMIversion*".nc")
    tracername = "d_"*region

    # Didn't use readfiled because dsfc is 2d.
    dsfc = ncread(file,tracername)
    b = BoundaryCondition(dsfc,γ.lon,γ.lat,γ.depth[1],3,1,γ.wet[:,:,1])
    return b
end

"""
Read vectors from mat file, translate to 3D,
 and save surface field to NetCDF file.
"""
function regions2nc(TMIversion,γ)

    # save in a different file.
    filenetcdf = pkgdatadir("regions_"*TMIversion*".nc")
    filemat = pkgdatadir("TMI_"*TMIversion*".mat")

    # region names
    # didn't figure out how to use an ordered dict, instead use a tuple
    list = ("GLOBAL","ANT","SUBANT",
            "NATL","NPAC","TROP","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

    regionname = Dict("GLOBAL" => "globally uniform",
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
         Dict("longname" => regionname[list[rr]]*" surface region", "units" => "[]"))
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
