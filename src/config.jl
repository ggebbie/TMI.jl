# mostly unchanging functions that help configure TMI
"""
    function config_from_nc(TMIversion)
    Configure TMI environment from NetCDF input file.
# Arguments
- `TMIversion`: TMI version for water-mass/circulation model
# Output
- `A`: TMI steady-state water-mass matrix
- `Alu`: LU decomposition of A
- `γ`: TMI grid properties
- `TMIfile`: TMI file name
"""
function config_from_nc(TMIversion)

    TMIfile = download_ncfile(TMIversion)

    println("A")
    @time A = watermassmatrix(TMIfile)

    # LU factorization for efficient matrix inversions
    println("Alu")
    @time Alu = lu(A)

    γ = Grid(TMIfile)
    
    # would be good to make this optional
    println("L=")
    @time L = circulationmatrix(TMIfile,A,γ)

    println("B=")
    @time B = boundarymatrix(TMIfile,γ)
    
    return  A, Alu, γ, TMIfile, L, B

end

"""
    function download_ncfile(TMIversion::String)

    Download NetCDF file for given TMI version

# Arguments
- `TMIversion::String`

# Output
- `TMIfile::String`: TMI file name

# Side-effect
- download TMI input file if necessary
"""
function download_ncfile(TMIversion::String)

    #make datdir() if it doesn't exist 
    !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 
    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    #if TMIfile doesn't exist, get GDrive url and download 
    if !isfile(TMIfile)

        # add a workaround for large files
        if TMIversion == "modern_180x90x33_GH11_GH12"
            println("workaround for 2° x 2°")
            shellscript = pkgsrcdir("read_nc_modern_180x90x33_GH11_GH12.sh")
            run(`sh $shellscript`)
            mv(joinpath(pwd(),"TMI_"*TMIversion*".nc"),TMIfile)
        elseif  TMIversion == "nordic_201x115x46_B23"
            println("workaround for regional Nordic Seas file")
            shellscript = pkgsrcdir("read_nc_nordic_201x115x46_B23.sh")
            run(`sh $shellscript`)
            TMIfile = TMIversion*".nc"
            mv(joinpath(pwd(),"TMI_"*TMIversion*".nc"),TMIfile)
        else
            println("read via GoogleDrive.jl")
            #- `url`: Google Drive URL for data
            url = ncurl(TMIversion)
            google_download(url,pkgdatadir())
        end
    end
    return TMIfile
end


"""
function gridinit(TMIfile)

    Construct the Grid given a file name

# Arguments
- `TMIfile::String`: NetCDF file name for TMI version

# Output
- `γ::Grid`: TMI grid struct
"""
function Grid(TMIfile::String; A = watermassmatrix(TMIfile))
    # get properties of grid
    lon,lat,depth = gridprops(TMIfile)

    
    # make ocean mask
    wet = wetmask(TMIfile,length(lon),length(lat),length(depth))

    # make interior mask
    interior = interiormask(A,wet,length(lon),length(lat),length(depth))

    # do not store: compute on demand
    #R = linearindex(wet)
    #γ = Grid(lon,lat,depth,I,R,wet)

    return Grid(lon,lat,depth,wet,interior)
end

function wetmask(TMIfile,nx,ny,nz)
    # read Cartesian Index from file.
    I = cartesianindex(TMIfile)
    # make a mask
    # first choice: smaller but inconsistent with input grid
    #wet = falses(maximum(I)[1],maximum(I)[2],maximum(I)[3])
    wet = falses(nx,ny,nz)
    wet[I] .= 1
    return wet
end

function interiormask(A,wet,nx,ny,nz)
    interior = falses(nx,ny,nz)
    I = cartesianindex(wet)
    list = findall(.!(isone.(sum(abs.(A),dims=2))))
    interior[I[list]] .= true 
    return interior
end

"""
Configure TMI environment from original MATLAB output
"""
function config_from_mat(TMIversion)
    # repetitive and long
    TMIfile = download_matfile(TMIversion)
    
    # # make a sample field from zyx cartesian indices
    Izyx = cartesianindex(TMIfile)

    # # make a mask
    wet = falses(maximum(Izyx)[1],maximum(Izyx)[2],maximum(Izyx)[3])
    wet[Izyx] .= 1

    I = cartesianindex(wet)
    R = linearindex(wet)

    Azyx = watermassmatrix(TMIfile)

    A = matrix_zyx2xyz(TMIfile,Azyx,R)

    # # LU factorization for efficient matrix inversions
    Alu = lu(A)

    γ = Grid(TMIfile,A=A)

    # # get properties of grid
    # lat,lon,depth = gridprops(TMIfile)

    # γ = Grid(lon,lat,depth,I,R,wet)

    # need to make this optional
    L = circulationmatrix(TMIfile,γ)
    
    B = boundarymatrix(TMIfile,γ)
    
    # consider re-ordering this.
    # some output should be optional
    # return Izyx or I or neither?
    #return  A, Alu, γ, TMIfile, I, L, B
    return  A, Alu, γ, TMIfile, L, B
end

"""
    function download_matfile(TMIversion::String)

    Download MATLAB file for given TMI version

# Arguments
- `TMIversion::String`

# Output
- `TMIfile::String`: TMI file name

# Side-effect
- download TMI input file if necessary
"""
function download_matfile(TMIversion::String)

    #- `url`: Google Drive URL for data
    url = maturl(TMIversion)
    TMIfile = pkgdatadir("TMI_"*TMIversion*".mat")
    TMIfilegz = TMIfile*".gz"
    println(TMIfile)
    !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 
    #    !isfile(TMIfilegz) & !isfile(TMIfile) ? google_download(url,pkgdatadir()) : nothing

    if TMIversion == "modern_180x90x33_GH11_GH12"
        println("workaround for 2° x 2°")
        shellscript = pkgsrcdir("read_mat_modern_180x90x33_GH11_GH12.sh")
        run(`sh $shellscript`)
        mv(joinpath(pwd(),"TMI_"*TMIversion*".mat.gz"),TMIfilegz,force=true)
    elseif  TMIversion == "nordic_201x115x46_B23"
        println("workaround for regional Nordic Seas file")
        shellscript = pkgsrcdir("read_mat_nordic_201x115x46_B23.sh")
        run(`sh $shellscript`)
        mv(joinpath(pwd(),"TMI_"*TMIversion*".mat.gz"),TMIfilegz,force=true)
    else
        !isfile(TMIfilegz) & !isfile(TMIfile) && google_download(url,pkgdatadir())
    end
    
    # cloak mat file in gz to get Google Drive spam filter to shut down
    isfile(TMIfilegz) & !isfile(TMIfile) && run(`gunzip $TMIfilegz`) 
    return TMIfile
end

"""
    function cartesianindex(file)
    Read and assemble the grid coordinates
    according to the legacy MATLAB code (z,y,x order).
# Arguments
- `file`: TMI NetCDF file name
# Output
- `I`: TMI Cartesian index for wet points
"""
function cartesianindex(file::String)
    # make the Cartesian tracer grid
    if file[end-1:end] == "nc"

        it = convert(Vector{Int},ncread(file,"i"))
        jt = convert(Vector{Int},ncread(file,"j"))
        kt = convert(Vector{Int},ncread(file,"k"))

        I = CartesianIndex.(it,jt,kt)

    elseif file[end-2:end] == "mat"

        matobj = matopen(file)
        if haskey(matobj,"it")
            it=convert(Vector{Integer},vec(read(matobj,"it")))
            jt=convert(Vector{Integer},vec(read(matobj,"jt")))
            kt=convert(Vector{Integer},vec(read(matobj,"kt")))
        elseif haskey(matobj,"i")
            it = convert(Vector{Integer},vec(read(matobj,"i")))
            jt = convert(Vector{Integer},vec(read(matobj,"j")))
            kt = convert(Vector{Integer},vec(read(matobj,"k")))
        elseif haskey(matobj,"grd")
            #grd = read(matobj,"grd")
            it = convert(Vector{Integer},vec(read(matobj,"grd")["it"]))
            jt = convert(Vector{Integer},vec(read(matobj,"grd")["jt"]))
            kt = convert(Vector{Integer},vec(read(matobj,"grd")["kt"]))
        else
            error("grid index key not found")
        end
        close(matobj)
        I = CartesianIndex.(it,jt,kt) 
    end
    return I
end

"""
    function gridprops(file)
    Read and assemble the grid properties.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function gridprops(file)
    if file[end-1:end] == "nc"
        
        lon = convert(Vector{Float64},ncread(file,"lon"))
        lat = convert(Vector{Float64},ncread(file,"lat"))
        depth = convert(Vector{Float64},ncread(file,"depth"))

    elseif file[end-2:end] == "mat"

        matobj = matopen(file)
        if haskey(matobj,"grd")
            lon = convert(Vector{Float64},vec(read(matobj,"grd")["LON"]))
            lat = convert(Vector{Float64},vec(read(matobj,"grd")["LAT"]))
            depth = convert(Vector{Float64},vec(read(matobj,"grd")["DEPTH"]))
        else
            lon=convert(Vector{Float64},vec(read(matobj,"LON")))
            lat=convert(Vector{Float64},vec(read(matobj,"LAT")))
            depth=convert(Vector{Float64},vec(read(matobj,"DEPTH")))
        end
        close(matobj)
    end
    return lon,lat,depth
end

"""
    function watermassmatrix(file)
    Read and assemble the water-mass matrix.
# Arguments
- `file`: TMI NetCDF or MATLAB file name
# Output
- `A`: water-mass matrix
"""
function watermassmatrix(file)

    # consider adding a catch if A doesn't exist in file.

    if file[end-1:end] == "nc"
        # Int or Integer?
        i = convert(Vector{Int},ncread(file,"Arow"))
        j = convert(Vector{Int},ncread(file,"Acol"))
        m = ncread(file,"m")
        A = sparse(i,j,m)
    elseif file[end-2:end] == "mat"
        matobj = matopen(file)
        if haskey(matobj,"A")
            A=read(matobj,"A")
            close(matobj)
        else
            close(matobj)
            return nothing
        end

        # But MATLAB had zyx format and we need xyz format.
        # linearindices R not available so will do conversion in higher scope

    end
    return A
end

"""
        function matrix_zyx2xyz(TMIfile,Azyx,γ)
   
    Transfer zyx format water-mass matrix A to xyz format
# Arguments
- `Azyx`: water-mass matrix in zyx format
- `γ`: TMI grid
# Output
- `Axyz`: water-mass matrix in xyz format
"""
function matrix_zyx2xyz(file,Azyx,R)

    izyx, jzyx, mzyx = findnz(Azyx)
    Izyx = cartesianindex(file)
        
    # Julia accounting x,y,z
    ixyz = updatelinearindex(izyx,Izyx,R)
    jxyz = updatelinearindex(jzyx,Izyx,R)
    
    # use grid indices to switch i,j values
    Axyz = sparse(ixyz,jxyz,mzyx)
    return Axyz
end

"""
    function circulationmatrix(file,γ)
    Read and assemble the circulation matrix from MATLAB.
    Transfer to updated x,y,z version
# Arguments
- `file`: TMI MATLAB file name
- `γ`: TMI grid
# Output
- `L`: circulation matrix in xyz format
"""
function circulationmatrix(file,γ)

    if file[end-2:end] == "mat" 

        matobj = matopen(file)
        if haskey(matobj,"L")
            # Matlab output in zyx format
            Lzyx=read(matobj,"L")
            close(matobj)

            Izyx = cartesianindex(file)
            izyx, jzyx, Fzyx = findnz(Lzyx)
            # Julia accounting x,y,z
            ixyz = updatelinearindex(izyx,Izyx,γ.R)
            jxyz = updatelinearindex(jzyx,Izyx,γ.R)
            L = sparse(ixyz,jxyz,Fzyx)

        else
            close(matobj)
            return nothing
        end

    elseif file[end-1:end] == "nc"

        # based on function arguments, read from inefficient storage of L matrix.
        if haskey(NCDataset(file),"F")

            i = convert(Vector{Int},ncread(file,"Lrow"))
            j = convert(Vector{Int},ncread(file,"Lcol"))
            F = ncread(file,"F")
            L = sparse(i,j,F)
        else
            return nothing
        end
    end
    
    return L
end

"""
    function circulationmatrix(file,A,γ)
    Read and assemble the circulation matrix from the efficient storage of A and F₀ variables. 
# Arguments
- `file`: TMI MATLAB file name
- `A`: TMI water-mass matrix
- `γ`: TMI grid
# Output
- `L`: circulation matrix in xyz format
"""
function circulationmatrix(file,A,γ)

    file[end-1:end] !== "nc" && error("not a NetCDF file")

    # based on function arguments, read F₀ to efficiently reproduce L matrix.

    if haskey(NCDataset(file),"F₀")
        F₀ = ncread(file,"F₀")
        F₀vec = F₀[γ.wet]
    
        # For each row of A, multiply by F₀
        i, j, F = findnz(A)

        # careful, this loop can be really slow
        for nn in eachindex(i)
            F[nn] *= F₀vec[i[nn]]
        end

        L = sparse(i,j,F)
        return L
    else
        return nothing
    end
    
end

"""
        function boundarymatrix(file,γ)
    Read and assemble the boundary matrix from MATLAB.
    Transfer to updated x,y,z version
# Arguments
- `file`: TMI MATLAB file name
- `γ`: TMI grid
# Output
- `B`: boundary condition matrix
"""
function boundarymatrix(file,γ)

    if file[end-2:end] == "mat"

        matobj = matopen(file)
        if haskey(matobj,"B")
            Bzyx=read(matobj,"B")
            close(matobj)

            # matlab in zyx format.
            # consider using Azyx2xyz here.
            Izyx = cartesianindex(file)
            izyx, jzyx, Fzyx = findnz(Bzyx)
            # for B, rows are 3D grid space, columns are for the surface index. 
            # Julia accounting x,y,z
            Isfc = surfaceindex(Izyx)
            ixyz = updatelinearindex(izyx,Izyx,γ.R)
            jxyz = updatelinearindex(Isfc[jzyx],Izyx,γ.R)

            # assume surface at k = 1 (revisit for LGM problem)
            # give the full dimension of sparse matrix
            B = sparse(ixyz,jxyz,Fzyx,sum(γ.wet),sum(γ.wet[:,:,1]))
        else
            close(matobj)
            return nothing
        end

    elseif file[end-1:end] == "nc"

        # based on function arguments, read from inefficient storage of L matrix.
        if haskey(NCDataset(file),"b")

            i = convert(Vector{Int},ncread(file,"Brow"))
            j = convert(Vector{Int},ncread(file,"Bcol"))
            b = ncread(file,"b")
            B = sparse(i,j,b,sum(γ.wet),sum(γ.wet[:,:,1]))

        else
            return nothing
        end
    end
    return B
end

"""
    function updatelinearindex(izyx,Izyx,R)
    Linear index translated from z,y,x to x,y,z accounting

    get Izyx Cartesian index stored from legacy MATLAB code

# Arguments
- `izyx`: index of interest in z,y,x accounting
- `Izyx`: wet Cartesian Index for z,y,x
- `R`: Linear indices for x,y,z 
# Output
- `ixyz`: index of interest in x,y,z accounting
"""
updatelinearindex(izyx,Izyx,R) = R[Izyx[izyx]]

"""
     function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    Read an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)
    Return a BoundaryCondition
"""
function surfaceregion(TMIversion::String,region::String,γ::Grid)::BoundaryCondition

    file = pkgdatadir("TMI_"*TMIversion*".nc")
    tracername = "d_"*region

    # Didn't use readfiled because dsfc is 2d.
    dsfc = ncread(file,tracername)
    b = BoundaryCondition(dsfc,γ.lon,γ.lat,γ.depth[1],3,1,γ.wet[:,:,1])
    return b
end

""" 
    function ncurl(TMIversion)
    placeholder function to give location (URL) of NetCDF Google Drive input
    in the future, consider a struct or Dict that describes all TMI versions.
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `url`: location (URL) for download
"""
function ncurl(TMIname)
    if TMIname == "modern_90x45x33_GH10_GH12"
        url = "https://docs.google.com/uc?export=download&id=1Fn_cY-90_RDbBGh6kV0kpXmsvwdjp1Cd"
    elseif TMIname == "modern_180x90x33_GH11_GH12"
        url = "https://docs.google.com/uc?export=download&id=1-YEkB_YeQGqPRH6kauhBb2bi_BjVGt9b"
    elseif TMIname == "modern_90x45x33_unpub12"
        url = "https://docs.google.com/uc?export=download&id=1Kw_Mr7fiKqan0nx0dKvGHnSInP0hQ7AV"
    elseif TMIname == "modern_90x45x33_G14"
        url = "https://docs.google.com/uc?export=download&id=1aeE7EXA-vy3Cm_drt4qCFw4AlpYrdudk"
    elseif TMIname == "modern_90x45x33_G14_v2"
        url = "https://docs.google.com/uc?export=download&id=1Mwhv70soBX6-pYijU0ElNl0TZw0vSbXN"
    elseif TMIname == "LGM_90x45x33_G14"
        url = "https://docs.google.com/uc?export=download&id=1yoDi7_foBt3TVULCstlWnNLHFc2G47Fz"  
    elseif TMIname == "LGM_90x45x33_G14A"
        url = "https://docs.google.com/uc?export=download&id=1ADkDI3Fc3z4Vm75K5u6hx0Yu1P0iVnW1"
    elseif TMIname == "LGM_90x45x33_GPLS1"
        url = "https://docs.google.com/uc?export=download&id=1VOrZGUsO7lp21qw6Yw0dBqGlIy_ImdRW"
    elseif TMIname == "LGM_90x45x33_GPLS2"
        url = "https://docs.google.com/uc?export=download&id=1cOCrty9kvA2s3NoD1QZjnNehlVbP0rHP"
    elseif TMIname == "LGM_90x45x33_OG18"
        url = "https://docs.google.com/uc?export=download&id=19zccG1BSdspD9rti2OttsF2Dm4P2OLjt"
    else
        url = nothing
    end
end

""" 
    function maturl(TMIversion)
    Find *mat file here.
    placeholder function to give location (URL) of Google Drive input
    in the future, consider a struct or Dict that describes all TMI versions.
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
# Output
- `url`: location (URL) for download
"""
function maturl(TMIname)
    if TMIname == "modern_90x45x33_GH10_GH12"
        #        url = "https://docs.google.com/uc?export=download&id=1Z2knDctAmZHO2lcWTBCdR8zjkmbcyCGg"
        url = "https://docs.google.com/uc?export=download&id=11-b4L6D1bnDdIIgSg8SaiCkOuPk64pSY"
    elseif TMIname == "modern_180x90x33_GH11_GH12"
        url = "https://docs.google.com/uc?export=download&id=11zD1nOfT6V7G0qIHdjK2pDGHFk-ExXwU"
    elseif TMIname == "modern_90x45x33_unpub12"
        url = "https://docs.google.com/uc?export=download&id=1sqkjFCPxZT_2Bm9rsp0acyxxkBri9YAT"
    elseif TMIname == "modern_90x45x33_G14"
        url = "https://docs.google.com/uc?export=download&id=1dCrDe5VXrsXiOf04mbuHID7xc5Ymm8-z"
    elseif TMIname == "modern_90x45x33_G14_v2"
        url = "https://docs.google.com/uc?export=download&id=1Axaqn88HZv3i4rYn0g5dUztVXlx_TmGt"
    elseif TMIname == "LGM_90x45x33_G14"
                url = "https://docs.google.com/uc?export=download&id=1qWnL9SrcjFRt4TUQ0k7v0jEFrKtUYpev"
    elseif TMIname == "LGM_90x45x33_G14A"
                url = "https://docs.google.com/uc?export=download&id=1E7qYWdAnoz4YXFvJ0AkVcqpHDq8aFP4R"
    elseif TMIname == "LGM_90x45x33_GPLS1"
                url = "https://docs.google.com/uc?export=download&id=1nF5KjWE3n--NkOTgIymIhP6pkafJ7jl3"
    elseif TMIname == "LGM_90x45x33_GPLS2"
                url = "https://docs.google.com/uc?export=download&id=1CHiR2J60HPRrxXTvJupyWEHmggMkgnxV"
    elseif TMIname == "LGM_90x45x33_OG18"
                url = "https://docs.google.com/uc?export=download&id=1WKAWxDYDFls6C5YzO1sMLDlUhC8_tKV2"
    else
        url = nothing
    end
    return url
end

""" 
Save TMI configuration to NetCDF format for non-proprietary access
"""
function config2nc(TMIversion,A,γ,L,B)

    # make new netcdf file.
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    isfile(filenetcdf) && rm(filenetcdf)

    grid2nc(TMIversion,γ)
    
    matfields2nc(TMIversion,γ)

    !isnothing(A) && watermassmatrix2nc(TMIversion,A)

    !isnothing(L) && circulationmatrix2nc(TMIversion,L,γ)

    !isnothing(B) && boundarymatrix2nc(TMIversion,B)

    #= is this part of the config? Or should it go to
     a separate output? It is similar to the output fields above. Probably should be considered part of the config. =#
    regions2nc(TMIversion,γ)

    optim2nc(TMIversion)

end

"""
Save grid dictionaries of attributes for writing to NetCDF file
"""
function griddicts(γ)
    # update names and types in dictionary
    
    TMIgrids = Dict("lon" => γ.lon,
                    "lat" => γ.lat,
                    "depth" => γ.depth)
    
    TMIgridsatts = Dict("lon" => Dict("longname" => "Longitude", "units" => "°E"),
                        "lat" => Dict("longname" => "Latitude", "units" => "°N"),
                        "depth" => Dict("longname" => "depth", "units" => "m"))

    return TMIgrids, TMIgridsatts

end

"""
Read 3D fields from mat file and save to NetCDF file.
"""
function matfields2nc(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    filemat = pkgdatadir("TMI_"*TMIversion*".mat")
    vars = matread(filemat)

    TMIgrids, TMIgridsatts = griddicts(γ)

    T = eltype(γ.lon) # does the eltype of longitude have to equal the tracer eltype?
    #T =  Float64

    varlist = Dict("dP"=>"qPO₄","q"=>"qPO₄",
                   "Tobs"=>"θ","Tmod"=>"θ","Tlgm"=>"θ",
                   "Terr"=>"σθ",
                   "Sobs"=>"Sₚ","Smod"=>"Sₚ","Slgm"=>"Sₚ",
                   "Serr"=>"σSₚ",
                   "O18obs"=>"δ¹⁸Ow","O18mod"=>"δ¹⁸Ow","O18lgm"=>"δ¹⁸Ow",
                   "O18err"=>"σδ¹⁸Ow",
                   "Pobs"=>"PO₄","Pmod"=>"PO₄","Plgm"=>"PO₄",
                   "Perr" => "σPO₄",
                   "Nobs"=>"NO₃","Nmod"=>"NO₃","Nlgm"=>"NO₃",
                   "Nerr" => "σNO₃",
                   "Oobs"=>"O₂","Omod"=>"O₂","Olgm"=>"O₂",
                   "Oerr"=>"σO₂",
                   "C13obs"=>"δ¹³C","C13mod"=>"δ¹³C","C13lgm"=>"δ¹³C",
                   "C13err" =>  "σδ¹³C")

    # iterate over all possible variables listed above
    Izyx = cartesianindex(filemat)
    TMIfields = Dict{String,Array{T,3}}()
    for (kk,vv) in varlist
        haskey(vars,kk) ? push!(TMIfields, vv => tracerinit(vars[kk], Izyx, γ.wet)) : nothing
    end

    # also save fields that are stored in the x struct, if they exist
    if haskey(vars,"x")
        for (kk,vv) in varlist
            println(kk)
            haskey(vars["x"],kk) ? push!(TMIfields, vv => tracerinit(vars["x"][kk], Izyx, γ.wet)) : nothing
        end
    end
    
    TMIfieldsatts = fieldsatts()

    # iterate in TMIgrids Dictionary to write to NetCDF.
    for (varname,varvals) in TMIfields
        
        nccreate(filenetcdf,varname,"lon",γ.lon,TMIgridsatts["lon"],"lat",γ.lat,TMIgridsatts["lat"],"depth",γ.depth,TMIgridsatts["depth"],atts=TMIfieldsatts[varname])
        println("write ",varname)
        ncwrite(varvals,filenetcdf,varname)

    end
end

"""
All variable names and attributes.
Useful for writing NetCDF files.
"""
fieldsatts() = 
    Dict("θ" => Dict("longname" => "potential temperature", "units" => "°C"),
         "σθ" => Dict("longname" => "1σ standard error in potential temperature", "units" => "°C"),
         "Sₚ" => Dict("longname" => "practical salinity", "units" => "PSS-78"),
         "σSₚ" => Dict("longname" => "1σ standard error in practical salinity", "units" => "PSS-78"),
         "δ¹⁸Ow" => Dict("longname" => "oxygen-18 to oxygen-16 ratio in seawater", "units" => "‰ VSMOW"),
         "σδ¹⁸Ow" => Dict("longname" => "1σ standard error in oxygen-18 to oxygen-16 ratio in seawater", "units" => "‰ VSMOW"),
         "PO₄" => Dict("longname" => "phosphate", "units" => "μmol/kg"),
         "σPO₄" => Dict("longname" => "1σ standard error in phosphate", "units" => "μmol/kg"),
         "qPO₄" => Dict("longname" => "local source of phosphate", "units" => "μmol/kg"),
         "NO₃" => Dict("longname" => "nitrate", "units" => "μmol/kg"),
         "σNO₃" => Dict("longname" => "1σ standard error in nitrate", "units" => "μmol/kg"),
         "O₂" => Dict("longname" => "dissolved oxygen", "units" => "μmol/kg"),
         "σO₂" => Dict("longname" => "1σ standard error in dissolved oxygen", "units" => "μmol/kg"),
         "δ¹³C" => Dict("longname" => "carbon-13 to carbon-12 ratio in DIC", "units" => "‰ PDB"),
         "σδ¹³C" => Dict("longname" => "1σ standard error fin carbon-13 to carbon-12 ratio in DIC", "units" => "‰ PDB"),
         "F₀" => Dict("longname" => "normalized mass flux out of gridcell", "units" => "(kg seawater/yr)/(kg gridcell)"))

"""
Read vectors from mat file, translate to 3D,
 and save surface field to NetCDF file.
"""
function regions2nc(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
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

function watermassmatrix2nc(TMIversion,A)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    i, j, m = findnz(A)
    nelements = length(i)

    # add the circulation matrix: problem can't store sparse matrix.
    varname= "m"
    elementatts = Dict("longname" => "TMI sparse matrix element number")
    matts =  Dict("longname" => "TMI water-mass fraction (sparse matrix values)", "units" => "(kg source)/(kg total)")
    nccreate(filenetcdf,varname,"A_element",1:nelements,elementatts,atts=matts)
    println("write ",varname)
    ncwrite(m, filenetcdf,varname)
    
     varname= "Arow"
     destatts = Dict("longname" => "gridcell number of destination (row value)")
     nccreate(filenetcdf,varname,"A_element",1:nelements,elementatts,atts=destatts)
    println("write ",varname)
    ncwrite(i, filenetcdf,varname)

    varname= "Acol"
    sourceatts = Dict("longname" => "gridcell number of source (column value)")
    nccreate(filenetcdf,varname,"A_element",1:nelements,elementatts,atts=sourceatts)
    println("write ",varname)
    ncwrite(j, filenetcdf,varname)

end

"""
Save optimization parameters to NetCDF file)

Future considerations: split into 2 functions
1) read from mat
2) save to nc
"""
function optim2nc(TMIversion)

    filemat = pkgdatadir("TMI_"*TMIversion*".mat")
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")

    matobj = matopen(filemat)
    if haskey(matobj,"fval")
        J̃ = read(matobj,"fval")
        iteratts = Dict("longname" => "iteration number")
        Jatts =  Dict("longname" => "cost function value", "units" => "[]")

        varname = "J" # J̃ not output to screen properly
        println("write ",varname)
        nccreate(filenetcdf,varname,"iter",1:length(J̃),iteratts,atts=Jatts)
    end
    if haskey(matobj,"u")
        ũ = read(matobj,"u")
        iteratts = Dict("longname" => "control element number")
        uatts =  Dict("longname" => "control vector", "units" => "[]")

        varname = "ũ"
        println("write ",varname)
        nccreate(filenetcdf,varname,"control_element",1:length(ũ),iteratts,atts=uatts)
    end
    close(matobj)

end

"""
Save circulation matrix `L` to NetCDF file.
"""
function circulationmatrix2nc(TMIversion,L,γ)

    T = eltype(L)
    fullmatrix = false # more efficient to just save F₀, then modify A to get L 
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    if !fullmatrix
        F₀ = tracerinit(γ.wet,T)
        for nd ∈ eachindex(γ.I)
            # normalized mass flux out of gridcell is found on diagonal
            F₀[γ.I[nd]] = -L[nd,nd]
        end

        TMIgrids, TMIgridsatts = griddicts(γ)
        TMIfieldsatts = fieldsatts()
        varname = "F₀"
        nccreate(filenetcdf,varname,"lon",γ.lon,TMIgridsatts["lon"],"lat",γ.lat,TMIgridsatts["lat"],"depth",γ.depth,TMIgridsatts["depth"],atts=TMIfieldsatts[varname])
        println("write ",varname)
        ncwrite(F₀,filenetcdf,varname)

    else # fullmatrix = true, makes a bigger nc file

        i, j, F = findnz(L)
        nelements = length(i)

        # add the circulation matrix: problem can't store sparse matrix.
        varname = "F"
        elementatts = Dict("longname" => "TMI sparse matrix element number")
        Fatts =  Dict("longname" => "normalized mass flux (sparse matrix values)","units" => "(kg seawater/s)/(kg gridcell)")
        nccreate(filenetcdf,varname,"L_element",1:nelements,elementatts,atts=Fatts)
        println("write ",varname)
        ncwrite(F, filenetcdf,varname)
        
        varname= "Lrow"
        destatts = Dict("longname" => "gridcell number of destination (row value)")
        nccreate(filenetcdf,varname,"L_element",1:nelements,elementatts,atts=destatts)
        println("write ",varname)
        ncwrite(i, filenetcdf,varname)

        varname= "Lcol"
        sourceatts = Dict("longname" => "gridcell number of source (column value)")
        nccreate(filenetcdf,varname,"L_element",1:nelements,elementatts,atts=sourceatts)
        println("write ",varname)
        ncwrite(j, filenetcdf,varname)
    end
    
end

"""
Save boundary matrix for transient model to NetCDF file
"""
function boundarymatrix2nc(TMIversion,B)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    i, j, b = findnz(B)
    nelements = length(i)

    # add the circulation matrix: problem can't store sparse matrix.
    varname= "b"
    elementatts = Dict("longname" => "TMI boundary matrix element number")
    matts =  Dict("longname" => "Boundary matrix values", "units" => "[]")
    nccreate(filenetcdf,varname,"B_element",1:nelements,elementatts,atts=matts)
    println("write ",varname)
    ncwrite(b, filenetcdf,varname)

    varname= "Brow"
    destatts = Dict("longname" => "gridcell number of 3D field")
    nccreate(filenetcdf,varname,"B_element",1:nelements,elementatts,atts=destatts)
    println("write ",varname)
    ncwrite(i, filenetcdf,varname)

    varname= "Bcol"
    sourceatts = Dict("longname" => "gridcell number of surface boundary condition")
    nccreate(filenetcdf,varname,"B_element",1:nelements,elementatts,atts=sourceatts)
    println("write ",varname)
    ncwrite(j, filenetcdf,varname)

end

"""
Put grid properties (Cartesian index) into NetCDF file
"""
function grid2nc(TMIversion,γ)

    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")

    linearindexatts = Dict("longname" => "linear index")
    nfld = length(γ.I)

    iatts =  Dict("longname" => "Cartesian index for x-direction", "units" => "[]")
    jatts =  Dict("longname" => "Cartesian index for y-direction", "units" => "[]")
    katts =  Dict("longname" => "Cartesian index for z-direction", "units" => "[]")

    varname = "i"
    nccreate(filenetcdf,varname,"linearindex",1:nfld,linearindexatts,atts=iatts)
    ncwrite(lonindex(γ.I),filenetcdf,varname)
    
    varname = "j"
    nccreate(filenetcdf,varname,"linearindex",1:nfld,linearindexatts,atts=jatts)
    ncwrite(latindex(γ.I),filenetcdf,varname)
    
    varname = "k"
    nccreate(filenetcdf,varname,"linearindex",1:nfld,linearindexatts,atts=katts)
    ncwrite(depthindex(γ.I),filenetcdf,varname)
    
end

