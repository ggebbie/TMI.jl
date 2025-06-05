# mostly unchanging functions that help configure TMI
"""
    function config(TMIversion; compute_lu = true)
    Configure TMI environment from NetCDF input file.
# Arguments
- `TMIversion`: TMI version for water-mass/circulation model
# Output
- `A`: TMI steady-state water-mass matrix
- `Alu`: LU decomposition of A
- `γ`: TMI grid properties
- `TMIfile`: TMI file name
"""
function config(TMIversion; compute_lu = true)

    TMIfile = download_file(TMIversion)

    println("A")
    @time A = watermassmatrix(TMIfile)

    # LU factorization for efficient matrix inversions
    println("Alu")
    if compute_lu
        @time Alu = lu(A)
    else
        Alu = nothing
    end

    γ = Grid(TMIfile)
    
    # would be good to make this optional
    println("L=")
    @time L = circulationmatrix(TMIfile,A,γ)

    println("B=")
    @time B = boundarymatrix(TMIfile,γ)
    
    return  A, Alu, γ, TMIfile, L, B

end

"""
    function download_file(TMIversion::String)

    Download NetCDF file for given TMI version

# Arguments
- `TMIversion::String`

# Output
- `TMIfile::String`: TMI file name

# Side-effect
- download TMI input file if necessary
"""
function download_file(TMIversion::String)

    #make datdir() if it doesn't exist 
    !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 
    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    #if TMIfile doesn't exist, get GDrive url and download 
    if !isfile(TMIfile)

        # add a workaround for large files
        if  TMIversion == "nordic_201x115x46_B23" || TMIversion == "modern_180x90x33_GH11_GH12" 
            println("use `Downloads.download` for large files")
            #shellscript = pkgsrcdir("read_nc_nordic_201x115x46_B23.sh")
            #run(`sh $shellscript`)
            #mv(joinpath(pwd(),"TMI_"*TMIversion*".nc"),TMIfile)
            url = ncurl(TMIversion)
            Downloads.download(TMI.ncurl(TMIversion),TMIfile)
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
    it = convert(Vector{Int},ncread(file,"i"))
    jt = convert(Vector{Int},ncread(file,"j"))
    kt = convert(Vector{Int},ncread(file,"k"))
    return CartesianIndex.(it,jt,kt)
end

"""
    function axislabels(file)
    Read and assemble the grid properties.
# Arguments
- `file`: TMI NetCDF file name
# Output
- `grid`: TMI grid coordinates
"""
function axislabels(file::String)
        lon = convert(Vector{Float64},ncread(file,"lon"))
        lat = convert(Vector{Float64},ncread(file,"lat"))
        depth = convert(Vector{Float64},ncread(file,"depth"))
    return lon, lat, depth
end

"""
function gridsize(TMIversion)

Parse the TMIversion string for the grid size

Will break if the prefix contains underscores or 'x'
"""
function gridsize(TMIversion)
    underscores = findall('_',TMIversion)
    exes = findall('x',TMIversion)
    Nx = TMIversion[underscores[1]+1:exes[1]-1]
    Ny = TMIversion[exes[1]+1:exes[2]-1]
    Nz = TMIversion[exes[2]+1:underscores[2]-1]
    return Nx,Ny,Nz
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
    # Int or Integer?
    i = convert(Vector{Int},ncread(file,"Arow"))
    j = convert(Vector{Int},ncread(file,"Acol"))
    m = ncread(file,"m")
    A = sparse(i,j,m)
    return A
end

function watermassmatrix(file, halflife)

    #!haskey(NCDataset(file),"F₀") && error("watermassmatrix: no rates included with this version")
    Aradio = watermassmatrix(file) # placeholder
    F₀ = ncread(file,"F₀")

    # decay constant
    λ = log(2)/halflife
    diag_decay = (- λ ./ F₀)[γ.wet]
    remove_inf = x -> (isinf(x) ? 0.0 : x) # use pair substitution instead
    replace!(remove_inf, diag_decay)
    return Aradio -= Diagonal(diag_decay)
end

"""
    function circulationmatrix(file,γ)
    Read and assemble the circulation matrix from NetCDF.

# Arguments
- `file`: TMI MATLAB file name
- `γ`: TMI grid
# Output
- `L`: circulation matrix in xyz format
"""
function circulationmatrix(file, γ)
    # based on function arguments, read from inefficient storage of L matrix.
    if haskey(NCDataset(file),"F")
        i = convert(Vector{Int},ncread(file,"Lrow"))
        j = convert(Vector{Int},ncread(file,"Lcol"))
        F = ncread(file,"F")
        L = sparse(i,j,F)
        return L
    else
        return nothing
    end
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

    # based on function arguments, read from inefficient storage of L matrix.
    if haskey(NCDataset(file),"b")

        i = convert(Vector{Int},ncread(file,"Brow"))
        j = convert(Vector{Int},ncread(file,"Bcol"))
        b = ncread(file,"b")
        B = sparse(i,j,b,sum(γ.wet),sum(γ.wet[:,:,1]))
        return B
    else
        return nothing
    end
end

"""
     function mixedlayermatrix(A, γ, τ)

Read and assemble the circulation matrix from the efficient storage of A and F₀ variables. 

# Arguments
- `A`: TMI water-mass matrix
- `γ`: TMI grid
- `τ`: uniform residence timescale (years) for all mixed layer points 
# Output
- `Lmix`: circulation matrix in xyz format for mixed layer points
"""
function mixedlayermatrix(A, γ, τ)
    Lmix = spzeros(size(A))
    I = γ.I # coordinates of all wet points, precompute this for speed
    mixedlayer = mixedlayermask(A,γ)
    for r in  1:size(A,1)
        if mixedlayer[I[r]]
            Lmix[r,:] = - A[r,:] / τ
        end
    end
    return Lmix
end

"""
     function dirichletmatrix(γ, τ)

Dirichlet surface boundary matrix with uniform timescale.
Assumes that the Dirichlet boundary condition is zero.

# Arguments
- `γ`: TMI grid
- `τ`: uniform restoring timescale (years) for all boundary points 
# Output
- `Ldir`: circulation matrix in xyz format for boundary points
"""
function dirichletmatrix(γ::Grid, τ)
    nfield = sum(γ.wet)
    Ldir = spzeros(nfield, nfield)
    boundary = boundarymask(γ)
    I = γ.I
    for r in  1:nfield
        if boundary[I[r]]
            Ldir[r,r] = - 1.0 / τ
        end
    end
    return Ldir
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
        url = "https://argo.whoi.edu/jake/TMI_modern_180x90x33_GH11_GH12.nc"
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
    elseif  TMIname == "nordic_201x115x46_B23"
        url = "https://argo.whoi.edu/jake/TMI_nordic_201x115x46_B23.nc"
    else
        url = nothing
    end
end


""" 
Save TMI configuration to NetCDF format for non-proprietary access
"""
function config2nc(TMIversion,A,γ,L,B)

    # make new netcdf file.
    filenetcdf = pkgdatadir("TMI_"*TMIversion*".nc")
    isfile(filenetcdf) && rm(filenetcdf)
    
    matfields2nc(TMIversion,γ)

    matsource2nc(TMIversion,γ)

    grid2nc(TMIversion,γ)

    !isnothing(A) && watermassmatrix2nc(TMIversion,A)

    !isnothing(L) && circulationmatrix2nc(TMIversion,L,γ)

    !isnothing(B) && boundarymatrix2nc(TMIversion,B)

    #= is this part of the config? Or should it go to
     a separate output? It is similar to the output fields above. Probably should be considered part of the config. =#
    regions_mat2nc(TMIversion,γ)

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
    
    TMIgridsatts = gridatts()

    return TMIgrids, TMIgridsatts
end

gridatts() = Dict("lon" => Dict("longname" => "Longitude", "units" => "°E"),
    "lat" => Dict("longname" => "Latitude", "units" => "°N"),
    "depth" => Dict("longname" => "Depth", "units" => "m"))


"""
All variable names and attributes.
Useful for writing NetCDF files.
"""
fieldsatts() = 
    Dict("Θ" => Dict("longname" => "Conservative Temperature", "units" => "°C"),
         "θ" => Dict("longname" => "potential temperature", "units" => "°C"),
         "σθ" => Dict("longname" => "1σ standard error in potential temperature", "units" => "°C"),
         "S⋆" => Dict("longname" => "preformed salinity", "units" => "g/kg"),
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

# help for reading foreign files, translate to TMI convention
tracerdict() = 
    Dict("Θ" => :Θ,
        "CT" => :Θ,
        "θ" => :θ,
        "THETA" => :θ,
        "theta" => :θ,
        "σθ" => :σθ,
        "S⋆" => Symbol("S⋆"),
        "Sₚ" => :Sₚ,
        "Sp" => :Sₚ,
        "SALT" => :Sₚ,
        "σSₚ" => :σSₚ,
        "δ¹⁸Ow" => :δ¹⁸Ow,
        "σδ¹⁸Ow" => :σδ¹⁸Ow,
         "PO₄" => :PO₄,
         "σPO₄" => :σPO₄,
         "qPO₄" => :qPO₄,
         "NO₃" => :NO₃,
         "σNO₃" => :σNO₃,
         "O₂" => :O₂,
         "σO₂" => :σO₂,
         "δ¹³C" => :δ¹³C,
         "σδ¹³C" => :σδ¹³C,
         "F₀" => :F₀)

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

""" 
    function tracerinit(wet,vec,I)
          initialize tracer field on TMI grid
        perhaps better to have a tracer struct and constructor
# Arguments
- `wet`:: BitArray mask of ocean points
- `vec`:: vector of values at wet points
- `I`:: Cartesian Index for vector
# Output
- `field`:: 3d tracer field with NaN on dry points
"""
function tracerinit(vec,I,wet)

    # preallocate
    T = eltype(vec)
    field = Array{T}(undef,size(wet))
    fill!(field,zero(T)/zero(T))    

    #- a comprehension
    [field[I[n]]=vec[n] for n ∈ eachindex(I)]
    return field
end
function sourceinit(vec,I,γ)

    # preallocate
    T = eltype(vec)
    source = Array{T}(undef,size(γ.wet))
    fill!(source,zero(T)/zero(T))    

    if length(vec) == sum(γ.wet)
        #- a comprehension
        [source[I[n]]=vec[n] for n ∈ eachindex(I)]
    elseif length(vec) == sum(γ.interior)
        # get new index that drops boundaries
        inew = findall(γ.interior[I])
        Inew = I[inew]
        [source[Inew[n]]=vec[n] for n ∈ eachindex(Inew)]
    else
        error("length of source vector doesn't match grid")
    end
    return source
end

versionlist() = ("modern_90x45x33_GH10_GH12",
    "modern_180x90x33_GH11_GH12",
    "modern_90x45x33_unpub12",
    "modern_90x45x33_G14", 
    "modern_90x45x33_G14_v2",				  
    "LGM_90x45x33_G14",
    "LGM_90x45x33_G14A",				  
    "LGM_90x45x33_GPLS1",				  
    "LGM_90x45x33_GPLS2",				  
    "LGM_90x45x33_OG18",
    "nordic_201x115x46_B23")
