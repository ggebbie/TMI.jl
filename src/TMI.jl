module TMI

using Revise
using LinearAlgebra, SparseArrays, NetCDF, Downloads,
    GoogleDrive, Distances, DrWatson, GibbsSeaWater,  
    PyPlot, PyCall, Distributions, Optim,
    Interpolations, LineSearches, MAT

export config, config_from_mat, config_from_nc,
    vec2fld, fld2vec, surfaceindex,
    lonindex, latindex, depthindex,
    surfacepatch, section,
    layerthickness, cellarea, cellvolume,
    planview, dyeplot, plotextent, tracerinit,
    watermassmatrix, watermassdistribution,
    circulationmatrix, boundarymatrixXYZ,
    linearindex, nearestneighbor, updatelinearindex,
    nearestneighbormask, horizontaldistance,
    readtracer, cartesianindex, Γ,
    costfunction_obs, costfunction_obs!,
    costfunction, costfunction!,
    trackpathways, regeneratedphosphate, volumefilled,
    surfaceorigin, sample_observations, steadyclimatology,
    steady_inversion,
    interpweights, interpindex,
    wetlocation, iswet,
    control2state, control2state!,
    sparsedatamap, config2nc

#Python packages - initialize them to null globally
#const patch = PyNULL()
#const ccrs = PyNULL()

# following example at ClimatePlots.jl
const mpl = PyNULL()
const plt = PyNULL()
const cmocean = PyNULL()
const cartopy = PyNULL()

#Initialize all Python packages - install with conda through Julia
function __init__()
    #copy!(patch, pyimport_conda("matplotlib.patches", "patches"))
    #copy!(ccrs, pyimport_conda("cartopy.crs", "ccrs"))

    # following ClimatePlots.jl
    copy!(mpl, pyimport_conda("matplotlib", "matplotlib", "conda-forge"))
    #copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge"))
    #copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))
    copy!(cartopy, pyimport_conda("cartopy", "cartopy", "conda-forge"))

    println("Python libraries installed")
 end

struct grid
    lon::Vector{Float64}
    lat::Vector{Float64}
    depth::Vector{Float64}
    I::Vector{CartesianIndex{3}} # index
    R::Array{Int,3}
#    R::LinearIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}} 
    wet::BitArray{3}
end

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

    #- `url`: Google Drive URL for data
    url = ncurl(TMIversion)
    
    TMIfile = datadir("TMI_"*TMIversion*".nc")
    println(url)
    println(TMIfile)
    !isdir(datadir()) ? mkpath(datadir()) : nothing
    !isfile(TMIfile) ? google_download(url,datadir()) : nothing

    #ncdata = NetCDF.open(TMIfile) # necessary?

    # read Cartesian Index from file.
    I = cartesianindex(TMIfile)

    # make a mask
    wet = falses(maximum(I)[1],maximum(I)[2],maximum(I)[3])
    #wet = BitArray{3}(undef,maximum(I)[1],maximum(I)[2],maximum(I)[3])
    #fill!(wet,0)
    wet[I] .= 1

    R = linearindex(wet)

    println("A")
    @time A = watermassmatrix(TMIfile)

    # LU factorization for efficient matrix inversions
    println("Alu")
    @time Alu = lu(A)
    
    # get properties of grid
    lat,lon,depth = gridprops(TMIfile)
    γ = grid(lon,lat,depth,I,R,wet)

    # would be good to make this optional
    println("L=")
    @time L = circulationmatrix(TMIfile,A,γ)

    println("B=")
    @time B = boundarymatrix(TMIfile,γ)
    
    return  A, Alu, γ, TMIfile, L, B

end

"""
Configure TMI environment from original MATLAB output
"""
function config_from_mat(TMIversion)

    #- `url`: Google Drive URL for data
    url = maturl(TMIversion)
    TMIfile = datadir("TMI_"*TMIversion*".mat")
    TMIfilegz = TMIfile*".gz"
    println(TMIfile)
    !isdir(datadir()) ? mkpath(datadir()) : nothing
    !isfile(TMIfilegz) & !isfile(TMIfile) ? google_download(url,datadir()) : nothing

    # cloak mat file in gz to get Google Drive spam filter to shut down
    isfile(TMIfilegz) & !isfile(TMIfile) ? run(`gunzip $TMIfilegz`) : nothing
    
    # move this to runtests.jl to see if it is read correctly?
    # Azyx = watermassmatrix(TMIfile) 

    # # make a sample field from zyx cartesian indices
    Izyx = cartesianindex(TMIfile)

    # # make a mask
    wet = BitArray{3}(undef,maximum(Izyx)[1],maximum(Izyx)[2],maximum(Izyx)[3])
    fill!(wet,0)
    wet[Izyx] .= 1

    # # consistent with tracer definition?
    # #wet = .!isnan.(c)

    I = cartesianindex(wet)

    R = linearindex(wet)

    Azyx = watermassmatrix(TMIfile)
    A = Azyx2xyz(TMIfile,Azyx,R)
    # #R = R[ wet ] # eliminate land points

    # # LU factorization for efficient matrix inversions
    Alu = lu(A)

    # get properties of grid
    lat,lon,depth = gridprops(TMIfile)

    γ = grid(lon,lat,depth,I,R,wet)

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
        haskey(matobj,"it") ? it=convert(Vector{Integer},vec(read(matobj,"it"))) : it = convert(Vector{Integer},vec(read(matobj,"i")))
        haskey(matobj,"jt") ? jt=convert(Vector{Integer},vec(read(matobj,"jt"))) : jt=convert(Vector{Integer},vec(read(matobj,"j")))
        haskey(matobj,"kt") ? kt=convert(Vector{Integer},vec(read(matobj,"kt"))) : kt=convert(Vector{Integer},vec(read(matobj,"k")))
        close(matobj)
        I = CartesianIndex.(it,jt,kt) # should this be reversed?
    end
    return I
end

"""
    function cartesianindex(wet)
    Read and assemble the grid coordinates
    according to a 3D tracer in x,y,z order
# Arguments
- `wet`: BitArray logical mask for wet points
# Output
- `I`: 3D Cartesian indices
"""
cartesianindex(wet::BitArray{3}) = findall(wet)

"""
    function linearindex(wet)
    Read and assemble the grid coordinates.
# Arguments
- `wet`: 3D mask for wet points
# Output
- `R`: array of linear indices, but not a LinearIndices type
"""
function linearindex(wet)
    R = Array{Int64,3}(undef,size(wet))
    fill!(R,0)
    # R = Array{Union{Int64,Nothing},3}(nothing,size(wet))
    R[wet]=1:sum(wet)
    # R = LinearIndices((1:maximum(it),1:maximum(jt),1:maximum(kt)));
    # R = LinearIndices((it,jt,kt));
    #Rwet = R[γ.wet]
    return R
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
        
        lat = convert(Vector{Float64},ncread(file,"lat"))
        lon = convert(Vector{Float64},ncread(file,"lon"))
        depth = convert(Vector{Float64},ncread(file,"depth"))

    elseif file[end-2:end] == "mat"
        
        matobj = matopen(file)
        lon=convert(Vector{Float64},vec(read(matobj,"LON")))
        lat=convert(Vector{Float64},vec(read(matobj,"LAT")))
        depth=convert(Vector{Float64},vec(read(matobj,"DEPTH")))
        close(matobj)

    end
    
    return lat,lon,depth
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
    if file[end-1:end] == "nc"
        # Int or Integer?
        i = convert(Vector{Int},ncread(file,"Arow"))
        j = convert(Vector{Int},ncread(file,"Acol"))
        m = ncread(file,"m")
        A = sparse(i,j,m)
    elseif file[end-2:end] == "mat"
        matobj = matopen(file)
        A=read(matobj,"A")
        close(matobj)

        # But MATLAB had zyx format and we need xyz format.
        # linearindices R not available so will do conversion in higher scope

    end
    return A
end

"""
        function Azyx2xyz(TMIfile,Azyx,γ)
   
    Transfer zyx format water-mass matrix A to xyz format
# Arguments
- `Azyx`: water-mass matrix in zyx format
- `γ`: TMI grid
# Output
- `Axyz`: water-mass matrix in xyz format
"""
function Azyx2xyz(file,Azyx,R)

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
        # Matlab output in zyx format
        Lzyx=read(matobj,"L")
        close(matobj)

        Izyx = cartesianindex(file)
        izyx, jzyx, Fzyx = findnz(Lzyx)
        # Julia accounting x,y,z
        ixyz = updatelinearindex(izyx,Izyx,γ.R)
        jxyz = updatelinearindex(jzyx,Izyx,γ.R)
        L = sparse(ixyz,jxyz,Fzyx)

    elseif file[end-1:end] == "nc"

        # based on function arguments, read from inefficient storage of L matrix.
        i = convert(Vector{Int},ncread(file,"Lrow"))
        j = convert(Vector{Int},ncread(file,"Lcol"))
        F = ncread(file,"F")
        L = sparse(i,j,F)
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

    elseif file[end-1:end] == "nc"

        # based on function arguments, read from inefficient storage of L matrix.
        i = convert(Vector{Int},ncread(file,"Brow"))
        j = convert(Vector{Int},ncread(file,"Bcol"))
        b = ncread(file,"b")
        B = sparse(i,j,b,sum(γ.wet),sum(γ.wet[:,:,1]))
    end

    return B
end

"""
    function updatelinearindex(izyx,Izyx,R)
    Linear index translated from z,y,x to x,y,z accounting
# Arguments
- `izyx`: index of interest in z,y,x accounting
- `Izyx`: wet Cartesian Index for z,y,x
- `R`: Linear indices for x,y,z 
# Output
- `ixyz`: index of interest in x,y,z accounting
"""
function updatelinearindex(izyx,Izyx,R)
    # get Izyx Cartesian index stored from legacy MATLAB code
    ixyz = R[Izyx[izyx]]
    return ixyz
end
                                                  
"""
    function readtracer(file,tracername)
    Read a tracer field from NetCDF.
# Arguments
- `file`: TMI NetCDF file name
- `tracername`: name of tracer
# Output
- `c`: 3D tracer field
"""
function readtracer(file,tracername)
    c = ncread(file,tracername)
    return c
end

"""
    Horizontal area of grid cell
"""
function cellarea(γ)
    dx = zonalgriddist(γ)
    dy = haversine((γ.lon[1],γ.lat[1])
                  ,(γ.lon[1],γ.lat[2]))

    area = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(area,0.0)

    # to calculate area everywhere
    #[area[i,j] = dx[j] * dy for i ∈ eachindex(γ.lon) for j ∈ eachindex(γ.lat)]

    # to calculate sea surface area
    I = γ.I
    [area[I[ii][1],I[ii][2]] = dx[I[ii][2]] * dy for ii ∈ eachindex(I) if I[ii][3] == 1]

    return area
end

"""
    Volume of each grid cell.
"""
function cellvolume(γ)
    dz = layerthickness(γ)
    area = cellarea(γ)
    volume = Array{Float64,3}(undef,length(γ.lon),length(γ.lat),length(γ.depth))
    fill!(volume,0.0)

    # for volume everywhere
    # [volume[i,j,k] = area[i,j] * dz[k] for i ∈ eachindex(γ.lon) for j ∈ eachindex(γ.lat) for k ∈ eachindex(γ.depth)]

    # for ocean volume only
    I = γ.I
    [volume[I[ii]] = area[I[ii][1],I[ii][2]] * dz[I[ii][3]] for ii ∈ eachindex(I)]
    return volume
end

function layerthickness(γ::grid)
    zface= (γ.depth[1:end-1].+γ.depth[2:end])./2;
    dz = ([zface[1] ; diff(zface); 500]);
    return dz
end

function zonalgriddist(γ::grid)
    dx = similar(γ.lat)
    for j in eachindex(γ.lat)
        dx[j] = haversine((γ.lon[1],γ.lat[j])
                         ,(γ.lon[2],γ.lat[j]))
    end
    return dx
end

"""
    function vec2fld
    Transfer a vector to a 3D field with accounting for ocean bathymetry
# Arguments
- `vector`: field in vector form (no land points)
- `I`: cartesian indices of ocean points
# Output
- `field`: field in 3d form including land points (NaN)
"""
function vec2fld(vector::Vector{T},I::Vector{CartesianIndex{3}}) where T<:Real

    # choose NaN for now, zero better? nothing better?
    fillvalue = zero(T)/zero(T) # NaN32 or NaN64

    nx = maximum(I)[1]
    ny = maximum(I)[2]
    nz = maximum(I)[3]

    # faster instead to allocate as undef and then fill! ?
    field = fillvalue .* zeros(nx,ny,nz)

    # a comprehension
    [field[I[n]]=vector[n] for n ∈ eachindex(I)]
    return field
end

"""
    function fld2vec
    Transfer 3D field with accounting for ocean bathymetry to a vector without land points.
    This is done more easily with a BitArray mask, i.e., vector = field[mask].
    This function may be removed in the future.
# Arguments
- `field`: field in 3d form including land points (NaN)
- `I`: cartesian indices of ocean points
# Output
- `vector`: field in vector form (no land points)
"""
function fld2vec(field::Array{T,3},I::Vector{CartesianIndex{3}}) where T<:Real
    vector = Vector{T}(undef,length(I))
    #- a comprehension
     [vector[n] = field[I[n]] for n ∈ eachindex(I)];
     return vector
 end

"""
    function surfacepatch
    Make a surface boundary condition
    with a rectangular patch
# Arguments
- `lonbox`: longitudes of box edges
- `latbox`: latitudes of box edges
- `γ`: TMI.grid
# Output
- `d`: vector that describes surface patch
"""
function surfacepatch(lonbox,latbox,γ::grid)

    # ternary operator to handle longitudinal wraparound
    lonbox[1] ≤ 0 ? lonbox[1] += 360 : nothing
    lonbox[2] ≤ 0 ? lonbox[2] += 360 : nothing

    # define the surface boundary condition

    # preallocate
    d = tracerinit(γ.wet,Float64)

    # can you add a logical to a Float64? yes, it's 1.0
    [d[i,j,1] =  latbox[1] ≤ γ.lat[j] ≤ latbox[2] && lonbox[1] ≤ γ.lon[i] ≤ lonbox[2] for i in eachindex(γ.lon) for j in eachindex(γ.lat)] 

    # old method for vectors
        #nfield = length(γ.I) # number of ocean points
    #d = zeros(Int,nfield) # preallocate
    #[d[n]=1 for n ∈ 1:nfield if γ.I[n][3]==1 && latbox[1] ≤ γ.lat[γ.I[n][2]] ≤ latbox[2]
    #     && lonbox[1] ≤ γ.lon[γ.I[n][1]] ≤ lonbox[2] ]
    return d
end

"""
    function nearestneighbormask
    Make a 3D tracer field that is 1 at location 
    of nearest neighbor, 0 elsewhere
# Arguments
- `loc`: location in a 3-tuple (lon,lat,depth)
- `γ`: TMI.grid
# Output
- `δ`: nearest neighbor mask 3D field
"""
function nearestneighbormask(loc,γ::grid,N=1)

    Inn, Rnn = nearestneighbor(loc,γ,N)

    # preallocate
    δ = falses(size(γ.wet))
    #Array{BitArray,3}(undef,size(γ.wet))
    #fill!(δ,zero(Bool))
    δ[Inn] = 1 
    return δ
end

"""
    function nearestneighbor(loc,γ)
    return the Cartesian index and linear index 
    of the nearest N neighbors
# Arguments
- `loc`: 3-tuple of lon,lat,depth location
- `γ`: TMI.grid
# Output
- `Inn`: Cartesian indices of nearest neighbor
#- `Rnn`: linear indices of nearest neighbor, Removed from code
"""
function nearestneighbor(loc,γ,N=1)

    xydist = horizontaldistance(loc[1:2],γ)

    if N==1
        ijdist,ijmin = findmin(xydist[γ.wet[:,:,1]])
        kdist,kmin = findmin(abs.(loc[3] .- γ.depth))
    elseif N > 1
        ijmin = sortperm(xydist[γ.wet[:,:,1]])
        kmin = sortperm(abs.(loc[3] .- γ.depth))
    end

    if N == 1
        Inn = CartesianIndex.(γ.I[ijmin][1],γ.I[ijmin][2],kmin)
    elseif N > 1
        Inn = Vector{CartesianIndex}(undef,N)
        cN2 = ceil(Integer,N/2)
        for ii in 1:cN2
            # translate ijmin into imin, jmin
            Inn[ii] = CartesianIndex.(γ.I[ijmin[ii]][1],γ.I[ijmin[ii]][2],kmin[1])
        end
        for ii in 1:floor(Integer,N/2)
            Inn[cN2+ii] = CartesianIndex.(γ.I[ijmin[ii]][1],γ.I[ijmin[ii]][2],kmin[2])
        end        
    end
    
    return Inn
end

"""
    function horizontaldistance(loc,γ)
    return the Cartesian index and linear index 
    of the nearest N neighbors
# Arguments
- `loc`: 3-tuple of lon,lat,depth location
- `γ`: TMI.grid
# Output
- `hordist`: horizontal distance to nearest tracer grid points
"""
function horizontaldistance(loc,γ::grid)

    # hordist will have same type as lon,lat,depth
    T = eltype(γ.lon)
    
    # pre-allocate horizontal distance
    hordist = Matrix{T}(undef,length(γ.lon),length(γ.lat))
    # will give NaN with appropriate precision
    fill!(hordist,zero(T)/zero(T))
    
    # calculate haversine horizontal distance on sphere
    [hordist[γ.I[ii]] = haversine((loc[1],loc[2]),                  (γ.lon[γ.I[ii][1]],γ.lat[γ.I[ii][2]]))
       for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
    return hordist
end

"""
    function section
    View latitude-depth slice of field
# Arguments
- `c`: 3d tracer field
- `lon`: longitude of section
- `γ`: TMI.grid
# Output
- `csection`: 2d slice of field
"""
function section(c,lon,γ)

    isec = findall(==(lon),γ.lon)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    csection= dropdims(view(c,isec,:,:),dims=1)
    return csection
end

function planview(c,depth,γ)

    isec = findall(==(depth),γ.depth)

    # use view so that a new array is not allocated
    # note: if cfld changes, so does csection (automatically)
    cplan = dropdims(view(c,:,:,isec),dims=3)
    return cplan
end

"""
    function plotextent
    Generate image showing user-specified ROI
# Arguments
- `latbox`: in format [lat_start, lat_stop]
- `lonbox`: in format [lon_start, lon_stop]

"""
function plotextent(latbox, lonbox)
    
#    ccrs = pyimport("cartopy.crs")
    lower_left = [minimum(lonbox), minimum(latbox)] #array of lower left of box

    #calc width and height of box
    w = maximum(lonbox) - minimum(lonbox)
    h = maximum(latbox) - minimum(latbox)

    #init GeoAxes
    fig = figure()
    ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())

    #plot rectangle
    ax.add_patch(TMI.mpl.patches.Rectangle(xy=lower_left,
                                 width=w, height=h,
                                 facecolor="blue",
                                 alpha=0.2,
                                 transform=TMI.cartopy.crs.PlateCarree()))
    #define extent of figure
    pad = 10 #how many deg lat and lon to show outside of bbox
    pad_add = [-pad, pad] #add this to latbox and lonbox
    padded_lat = latbox + pad_add
    padded_lon = lonbox + pad_add
    ext = vcat(padded_lon, padded_lat) #make into one vector
    ax.set_extent(ext)

    # using cartopy 0.18 and NaturalEarth is missing
    ax.coastlines() #show coastlines

    #add gridlines
    gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
    gl.top_labels = false
    gl.right_labels = false

    ax.set_title("User-defined surface patch")
end

"""
    function dyeplot
    Plot of dye in ocean
# Arguments
- `lat`: latitude arrays
- `depth`: depth array
- `vals`: lat x depth value array
- `lims`: contour levels
"""
function dyeplot(lat, depth, vals, lims)

    #calc fignum - based on current number of figures
    figure()
    contourf(lat, depth, vals, lims) 
    gca().set_title("Meridional dye concentration")
end

"""
    function depthindex(I) 
    Get the k-index (depth level) from the Cartesian index
"""
function depthindex(I)
    T = eltype(I[1])
    k = Vector{T}(undef,length(I))
    [k[n]=I[n][3] for n ∈ eachindex(I)]
    return k
end

"""
    function lonindex(I) 
    Get the i-index (lon index) from the Cartesian index
"""
function lonindex(I)
    T = eltype(I[1])
    i = Vector{T}(undef,length(I))
    [i[n]=I[n][1] for n ∈ eachindex(I)]
    return i
end

"""
    function latindex(I) 
    Get the j-index (latitude index) from the Cartesian index
"""
function latindex(I)
    T = eltype(I[1])
    j = Vector{T}(undef,length(I))
    [j[n]=I[n][2] for n ∈ eachindex(I)]
    return j
end

"""
    function surfaceindex(I) 
    Get the vector-index where depth level == 1 and it is ocean.
"""
function surfaceindex(I)
    Isfc = findall(depthindex(I) .==1)
    return Isfc
end

""" 
    function tracerinit(wet,ltype=Float64)
      initialize tracer field on TMI grid
    perhaps better to have a tracer struct and constructor
# Arguments
- `wet`::BitArray mask of ocean points
- `ltype`:: optional type argument, default=Float64
# Output
- `d`:: 3d tracer field with NaN on dry points
"""
function tracerinit(wet,ltype=Float64)
    # preallocate
    d = Array{ltype}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    d[wet] .= zero(ltype)
    d[.!wet] .= zero(ltype)/zero(ltype) # NaNs with right type
    return d
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

""" 
    function control2state(tracer2D,γ)
    turn 2D surface field into 3D field with zeroes below surface    
# Arguments
- `tracer2D`:: 2D surface tracer field
- `wet`::BitArray mask of ocean points
# Output
- `tracer3D`:: 3d tracer field with NaN on dry points
"""
function control2state(tracer2D::Matrix{T},wet) where T<: Real
    # preallocate
    tracer3D = Array{T}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer3D[wet] .= zero(T)
    tracer3D[.!wet] .= zero(T)/zero(T)
    tracer3D[:,:,1] = tracer2D
    return tracer3D
end

""" 
    function control2state(u,γ)
    turn surface control vector into 3D field with zeroes below surface    
# Arguments
- `u`:: surface control vector
- `wet`::BitArray mask of ocean points
# Output
- `tracer3D`:: 3d tracer field with NaN on dry points
"""
function control2state(u::Vector{T},wet) where T<: Real
    # preallocate
    tracer3D = Array{T}(undef,size(wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer3D[wet] .= zero(T)
    tracer3D[.!wet] .= zero(T)/zero(T)
    tracer3D[:,:,1][wet[:,:,1]] = u
    return tracer3D
end

""" 
    function control2state!(c,u,γ)
    Add surface control vector to existing 3D field 
# Arguments
- `c`:: state field, 3d tracer field with NaN on dry points, modified by function
- `u`:: surface control vector
- `wet`::BitArray mask of ocean points
"""
function control2state!(c::Array{T,3},u::Vector{T},γ) where T<: Real
    #c[:,:,1][wet[:,:,1]] .+= u # doesn't work
#    [c[γ.I[ii][1],γ.I[ii][2],γ.I[ii][3]] += u[ii] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
    list = surfaceindex(γ.I)
    [c[γ.I[ii]] += u[list[ii]] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
end

""" 
    function control2state!(c,u,γ)
    Add surface control vector to existing 3D field 
# Arguments
- `c`:: state field, 3d tracer field with NaN on dry points, modified by function
- `u`:: surface control vector
- `wet`::BitArray mask of ocean points
"""
function control2state!(c::Vector{T},u::Vector{T},γ) where T<: Real
    list = surfaceindex(γ.I)
    [c[ii] += u[list[ii]] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
end

function state2obs(cvec,wis,γ)
    # interpolate onto data points
    N = length(wis)
    sumwis = Vector{Float64}(undef,N)
    list = vcat(1:length(γ.lon),1)

    # perhaps the most clever line in TMI.jl?
    wetwrap = view(γ.wet,list,:,:)

    [sumwis[i] = wetwrap[wis[i]...] for i in eachindex(wis)]

    # reconstruct the observations
    ỹ = Vector{Float64}(undef,N)
    c̃ = tracerinit(γ.wet)
    c̃[γ.wet] = cvec
    replace!(c̃,NaN=>0.0)
    [ỹ[i] = c̃[wis[i]...]/sumwis[i] for i in 1:N]
    return ỹ
end
    
""" 
    function trackpathways(TMIversion,latbox,lonbox)
    Track the pathways of a user-defined water mass.
     Steps: (a) define the water mass by a rectangular surface patch dyed with passive tracer concentration of         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.
     See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `latbox`: min and max latitude of box
- `lonbox`: min and max longitude of box
- `γ`: TMI grid
# Output
- `c`: fraction of water from surface source
"""
function trackpathways(Alu,latbox,lonbox,γ)

    d = surfacepatch(lonbox,latbox,γ)

    # do matrix inversion to get quantity of dyed water throughout ocean:
    c = tracerinit(γ.wet); # pre-allocate c

    # make methods that make the "wet" index unnecessary
    c[γ.wet] = Alu\d[γ.wet] # equivalent but faster than `c = A\d`

    return c
end

""" 
    function watermassdistribution(TMIversion,latbox,lonbox)
    Track the pathways of a user-defined water mass.
     Steps: (a) define the water mass by an oceanographically-relevant surface patch dyed with passive tracer concentration of one
         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.
     See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `region`: name of pre-defined surface region
- `γ`: TMI grid
# Output
- `g`: water-mass fraction
"""
function watermassdistribution(TMIversion,Alu,region,γ)

    d = surfaceregion(TMIversion,region,γ)

    # do matrix inversion to get quantity of dyed water throughout ocean:
    g = tracerinit(γ.wet); # pre-allocate c

    # make methods that make the "wet" index unnecessary
    g[γ.wet] = Alu\d[γ.wet] # equivalent but faster than `c = A\d`

    return g
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
    elseif TMIname == "modern_180x90x33_GH10_GH12"
        url = "https://docs.google.com/uc?export=download&id=1-YEkB_YeQGqPRH6kauhBb2bi_BjVGt9b"
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
        url = "https://docs.google.com/uc?export=download&id=1qPRq7sonwdjkPhpMcAnvuv67OMZSBqgh"
    elseif TMIname == "modern_180x90x33_GH10_GH12"
        url = "https://docs.google.com/uc?export=download&id=11zD1nOfT6V7G0qIHdjK2pDGHFk-ExXwU"
    else
        url = nothing
    end
    return url
end

""" 
    function regeneratedphosphate(TMIversion,Alu,γ)
    Regenerated (i.e., accumulated, remineralized) phosphate
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `PO₄ᴿ`: regenerated phosphate
"""
function regeneratedphosphate(TMIversion,Alu,γ)

    inputfile = datadir("TMI_"*TMIversion*".nc")
        
    #A, Alu, γ, inputfile = config(TMIversion)
    qPO₄ = readtracer(inputfile,"qPO₄")

    # PO₄ᴿ = cumulative regenerated phosphate
    PO₄ᴿ = tracerinit(γ.wet); # pre-allocate 
    PO₄ᴿ[γ.wet] = -(Alu\qPO₄[γ.wet])
    return PO₄ᴿ
end

""" 
    function volumefilled(TMIversion)
    Find the ocean volume that has originated from each surface box.
     This is equivalent to solving a sensitivity problem:
     The total volume is V = vᵀ c , where v is the volume of each box 
     and c is the fraction of volume from a given source which
     satisfies the equation A c = d.                     
     Next, dV/d(d) = A⁻ᵀ v, and dV/d(d) is exactly the volume originating from each source.

     See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011. 
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI.grid
# Output
- `volume`: global ocean volume filled by a surface region
"""
function volumefilled(TMIversion,Alu,γ)

    #A, Alu, γ = config(TMIversion)
    
    v = cellvolume(γ)
    area = cellarea(γ)
    
    # effectively take inverse of transpose A matrix.
    dVdd = tracerinit(γ.wet); # pre-allocate c
    dVdd[γ.wet] = Alu'\v[γ.wet]

    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    volume = zeros(Float64,length(γ.lon),length(γ.lat))
    #volume = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    #fill!(volume,0.0)

    # this step could use a function with γ.I argument
    [volume[I[ii][1],I[ii][2]] = dVdd[I[ii]] / area[I[ii][1],I[ii][2]] for ii ∈ eachindex(I) if I[ii][3] == 1]

    return volume
end

""" 
    function surfaceorigin(TMIversion,loc)
     Find the surface origin of water for some interior box 
     This is equivalent to solving a sensitivity problem:
     The mass fraction at a location `loc` of interest is 
    `c[loc] = δᵀ c`, where `δ` samples the location of the global mass-fraction variable, c.
    Then the sensitivity of `c[loc]` is: d(c[loc])/d(d) = A⁻ᵀ δ.
    The derivative is solved using the constraint: Ac = d.
    The sensitivity is exactly the mass fraction originating from each source.      
    This problem is mathematically similar to determining how the ocean is filled.
# Arguments
- `loc`: location (lon,lat,depth) of location of interest
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `origin`: surface map of fraction of source water for a given location
"""
function surfaceorigin(loc,Alu,γ)

    #A, Alu, γ = config(TMIversion)

    ctmp = tracerinit(γ.wet)
    δ = interpweights(loc,γ)
    
    # Find nearest neighbor on grid
    # set δ = 1 at grid cell of interest
    #δ = nearestneighbormask(loc,γ)
    # Note: ctrue[γ.wet]'*δ[γ.wet] returns interpolated value

    dvlocdd = tracerinit(γ.wet); # pre-allocate c
    dvlocdd[γ.wet] = Alu'\δ[γ.wet]

    # origin is defined at sea surface
    origin = view(dvlocdd,:,:,1)
    return origin
end

"""
function interpindex(loc,γ)
    Weights for linear interpolation.
    The derivative of linear interpolation is needed in sensitivity studies.
    ReverseDiff.jl could find this quantity automatically.
    Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.
# Arguments
- `c`: a temporary tracer field, would be nice to make it unnecessary
- `loc`: (lon,lat,depth) tuple of a location of interest
- `γ`: TMI grid
# Output
- `δ`: weights on a 3D tracer field grid
"""
function interpindex(loc,γ)

    # Handle longitudinal periodic condition (i.e., wraparound)
    lon = vcat(copy(γ.lon),γ.lon[1]+360.)
    list = vcat(1:length(γ.lon),1)
    nodes = (lon,γ.lat,γ.depth)

    # eliminate need to pass tracer value
    wis = Interpolations.weightedindexes((Interpolations.value_weights,),((Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))),nodes, loc)

    # issue, some of weighted points may be NaNs in tracer field
    # handle this in the Interpolations.jl routines
    # may involve chaging Gridded(Linear()) above
    return wis
end

"""
function interpweights(loc,γ)
    Weights for linear interpolation.
    The derivative of linear interpolation is needed in sensitivity studies.
    ReverseDiff.jl could find this quantity automatically.
    Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.
# Arguments
- `loc`: (lon,lat,depth) tuple of a location of interest
- `γ`: TMI grid
# Output
- `δ`: weights on a 3D tracer field grid
"""
function interpweights(loc,γ)

    # handle wraparound
    # repeated (unnecessarily?) in interpindex
    lon = vcat(copy(γ.lon),γ.lon[1]+360)
    list = vcat(1:length(γ.lon),1)

    wis = interpindex(loc,γ)

    # translate to weights via
    #http://juliamath.github.io/Interpolations.jl/latest/devdocs/
    δ = tracerinit(γ.wet)

    # changes in δwrap i=91 are translated back to δ i=1
    δwrap = view(δ,list,:,:)
    for ii = 1:2
        for jj = 1:2
            for kk = 1:2
                δwrap[wis[1].istart+ii-1,wis[2].istart+jj-1,wis[3].istart+kk-1] +=
                wis[1].weights[ii]*wis[2].weights[jj]*wis[3].weights[kk]
            end
        end
    end

    # if some adjacent points are dry, then re-normalize to keep this interpolation as an average.
    # The hope is that the interpolation is stable with this approach, but other side effects are likely.
    # Note that this should be handled earlier, like in the interpindex section. For this reason, there could be an inconsistency in the global map function.
    if iszero(sum(filter(!isnan,δ)))
        δ = nothing
    elseif sum(filter(!isnan,δ)) < 1.0
        δ ./= sum(filter(!isnan,δ))
    end
    
    return δ
end

"""
function steadyclimatology(u₀,Alu,y,d₀,W⁻,γ)
     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u₀
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)
    subject to Aỹ = d₀ + Γ u₀.                 
    W⁻ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u₀`:
- `Alu`:
- `d₀`: first guess of boundary conditions and interior sources
- `y`: observations on 3D grid
- `W⁻`: weighting matrix best chosen as inverse error covariance matrix
- `fg!`: compute cost function and gradient in place
- `γ`: grid
"""
function steadyclimatology(u₀,Alu,d₀,y,W⁻,fg!,γ)

    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(),Optim.Options(show_trace=true, iterations = 5))

    return out    
end

"""
function sparsedatamap(u₀,Alu,y,d₀,W⁻,γ)
     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u₀
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)
    subject to Aỹ = d₀ + Γ u₀.                 
    W⁻ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u₀`:
- `Alu`:
- `d₀`: first guess of boundary conditions and interior sources
- `y`: observations on 3D grid
- `W⁻`: weighting matrix best chosen as inverse error covariance matrix
- `fg!`: compute cost function and gradient in place
- `γ`: grid
"""
#function sparsedatamap(u₀,fg!)
function sparsedatamap(u₀,Alu,d₀,y,W⁻,wis,locs,Q⁻,γ)

    # ### added this
     fg!(F,G,x) = costfunction!(F,G,x,Alu,d₀,y,W⁻,wis,locs,Q⁻,γ)
    
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(linesearch = LineSearches.BackTracking()),Optim.Options(show_trace=true, iterations = 5))
#    out = optimize(Optim.only_fg!(fg!), u₀, GradientDescent(),Optim.Options(show_trace=true, iterations = 5))

    return out    
end

""" 
    function sample_observations(TMIversion,variable)
    Synthetic observations that are a contaminated version of real observations
    This version: gridded observations
# Arguments
- `TMIversion::String`: version of TMI water-mass/circulation model
- `variable::String`: variable name to use as template
# Output
- `y`: contaminated observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `θtrue`: real observations, 3D field
"""
function sample_observations(TMIversion,variable,γ)

    inputfile = datadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    θtrue = readtracer(inputfile,variable)
    σθ = readtracer(inputfile,"σ"*variable)

    ntrue = tracerinit(γ.wet)
    ntrue[γ.wet] = rand(Normal(),length(σθ[γ.wet])) .* σθ[γ.wet]

    y = θtrue .+ ntrue

    # get cost function (J) based on model misfit
    # here the data-model misfit is weighted by the expected error

    # weighting matrix
    #W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    W⁻ = (1/sum(γ.wet)) .* Diagonal(1 ./σθ[γ.wet].^2)
    return y, W⁻, θtrue
end
 
""" 
    function sample_observations(TMIversion,variable,locs)
    Synthetic observations that are a contaminated version of real observations
    This version: observations with random (uniform) spatial sampling
# Arguments
- `TMIversion::String`: version of TMI water-mass/circulation model
- `variable::String`: variable name to use as template
- `N`: number of observations
# Output
- `y`: contaminated observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `ytrue`: uncontaminated observations, 3D field
- `locs`: 3-tuples of locations for observations
- `wis`: weighted indices for interpolation to locs sites
"""
function sample_observations(TMIversion,variable,γ,N)

    inputfile = datadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    
    θtrue = readtracer(inputfile,variable)
    replace!(θtrue,NaN=>0.0)
    σθ = readtracer(inputfile,"σ"*variable)
    replace!(σθ,NaN=>0.0)

    # get random locations that are wet (ocean)
    locs = Vector{Tuple{Float64,Float64,Float64}}(undef,N)
    [locs[i] = wetlocation(γ) for i in 1:N]

    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    # look at total weight, < 1 if there are land points
    # later make sure total weight = 1 for proper average
    sumwis = Vector{Float64}(undef,N)
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)
    [sumwis[i] = wetwrap[wis[i]...] for i in 1:N]

    # sample the true field at these random locations
    ytrue = Vector{Float64}(undef,N)
    replace!(θtrue,NaN=>0.0)
    θwrap = view(θtrue,list,:,:)
    [ytrue[i] = θwrap[wis[i]...]/sumwis[i] for i in 1:N]

    # interpolate the standard deviation of expected error
    σtrue = Vector{Float64}(undef,N)
    σwrap = view(σθ,list,:,:)
    [σtrue[i] = σwrap[wis[i]...]/sumwis[i] for i in 1:N]

    ntrue = rand(Normal(),N) .* σtrue
    y = ytrue .+ ntrue

    # weighting matrix
    #W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    W⁻ = (1/N) .* Diagonal(1 ./σtrue.^2)
    return y, W⁻, ytrue, locs, wis
end
 
""" 
    function costfunction_obs(u,Alu,dfld,yfld,Wⁱ,γ)
    squared model-data misfit for gridded data
    controls are a vector input for Optim.jl
# Arguments
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `y`: observations on grid
- `d`: model constraints
- `Wⁱ`: inverse of W weighting matrix for observations
- `wet`: BitArray ocean mask
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
function costfunction_obs(u::Vector{T},Alu,dfld::Array{T,3},yfld::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},γ::grid) where T <: Real
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    # below surface = 0 % no internal sinks or sources.
    #u = tracerinit(γ.wet[:,:,1],T)
    #u[γ.wet[:,:,1]] = uvec

    d = dfld[γ.wet]
    y = view(yfld,γ.wet)
    
    #ỹ = tracerinit(γ.wet,T)
    #n = tracerinit(γ.wet,T)
    #dJdn = tracerinit(γ.wet,T)
    dJdd = tracerinit(γ.wet,T)
    
    # first-guess reconstruction of observations
    #Δd = d + Γ(u,wet)
    #ỹ[wet] =  Alu\Δd[wet]
    #n = y .- ỹ
    #d[γ.wet] = Alu\d[γ.wet]

    # use in-place functions to make this more performant
    control2state!(d,u,γ) # d stores Δd
    ldiv!(Alu,d) # d stores -ỹ
    d .-= y # d stores n
    J = d'* (Wⁱ * d)
    # move this to its own function
    #J += u[wet[:,:,1]]'* (Qⁱ * u[wet[:,:,1]])

    #dJdn[wet] = 2Wⁱ*[wet]
    #dJdd[wet] = Alu'\dJdn[wet]
    gd = 2*(Wⁱ*d)
    ldiv!(Alu',gd)

    # "transpose" of control2state! operation
    gJ = gd[surfaceindex(γ.I)]

    return J, gJ
end

""" 
    function costfunction_obs(u,Alu,y,d,Wⁱ,wet)
    squared model-data misfit for gridded data
    controls are a vector input for Optim.jl
# Arguments
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `d`: model constraints
- `y`: observations on grid
- `Wⁱ`: inverse of W weighting matrix for observations
- `γ`: grid
"""
function costfunction_obs!(J,gJ,u::Vector{T},Alu,dfld::Array{T,3},yfld::Array{T,3},Wⁱ::Diagonal{T, Vector{T}},γ::grid) where T <: Real

    d = dfld[γ.wet]
    y = view(yfld,γ.wet)
    control2state!(d,u,γ) # d stores Δd
    ldiv!(Alu,d)
    d .-= y # .- d

    if gJ != nothing
        gd = Vector{T}(undef,length(d))
        gd = 2Wⁱ*d
        ldiv!(Alu',gd)

        list = surfaceindex(γ.I)
        [gJ[ii] = gd[list[ii]] for ii in 1:length(list)]
        #gJ = gd[surfaceindex(γ.I)]
        #[gJ[γ.R[ii]] = gd[ii] for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]

        #pick out I[3]==1
        #gJ = -dJdd[:,:,1][γ.wet[:,:,1]]
    end
    
    if J !=nothing
        return  d'* (Wⁱ * d)       
    end
end

""" 
    function costfunction_obs(u,Alu,dfld,yfld,Wⁱ,wis,locs,γ)
    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
# Arguments
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `y`: pointwise observations
- `d`: model constraints
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation 
- `locs`: data locations
- `γ`: grid
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
function costfunction_obs(u::Vector{T},Alu,dfld::Array{T,3},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,γ::grid) where T <: Real

    d = dfld[γ.wet] # couldn't use view b.c. of problem with function below

    # use in-place functions: more performant
    control2state!(d,u,γ) # d stores Δd
    ldiv!(Alu,d) # d stores c̃

    ỹ = state2obs(d,wis,γ)
    ỹ .-= y # stores n, data-model misfit
    J = ỹ'* (Wⁱ * ỹ)
    #println(J)
    gỹ = 2*(Wⁱ*ỹ)

    #gd = Array{T,3}(undef,size(dfld))
    #gd = Vector{T}(undef,sum(γ.wet))
    gd = zeros(T,sum(γ.wet))
    
    # transpose of "E" operation in state2obs
    for ii in eachindex(y)
        # interpweights repeats some calculations
        gd .+= gỹ[ii] * interpweights(locs[ii],γ)[γ.wet]
    end
    # do Eᵀ gỹ 
    ldiv!(Alu',gd)
    list = surfaceindex(γ.I)
    gJ = Vector{T}(undef,sum(γ.wet[:,:,1]))
    [gJ[ii] = gd[list[ii]] for ii in eachindex(list)]
    return J, gJ
end

""" 
    function costfunction_obs!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,locs,γ)
    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
# Arguments
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `dfld`: model constraints
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `γ`: grid
"""
function costfunction_obs!(J,gJ,u::Vector{T},Alu,dfld::Array{T,3},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,γ::grid) where T <: Real

    d = dfld[γ.wet] # couldn't use view b.c. of problem with function below

    # use in-place functions: more performant
    control2state!(d,u,γ) # d stores Δd
    ldiv!(Alu,d) # d stores c̃

    ỹ = state2obs(d,wis,γ)
    ỹ .-= y # stores n, data-model misfit

    if gJ != nothing    
        gỹ = 2*(Wⁱ*ỹ)

        #gd = Array{T,3}(undef,size(dfld))
        #gd = Vector{T}(undef,sum(γ.wet))
        gd = zeros(T,sum(γ.wet))
        for ii in eachindex(y)
            # interpweights repeats some calculations
            gd .+= gỹ[ii] * interpweights(locs[ii],γ)[γ.wet]
        end
        # do Eᵀ gỹ 
        ldiv!(Alu',gd)
        list = surfaceindex(γ.I)
        [gJ[ii] = gd[list[ii]] for ii in eachindex(list)]

    end

    if J != nothing
        return  ỹ'* (Wⁱ * ỹ)
    end
end

""" 
    function costfunction(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,Q⁻,γ)
    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue: couldn't figure out how to nest with costfunction_obs!
# Arguments
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `dfld`: model constraints
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `Q⁻`: weights for control vector
- `γ`: grid
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
function costfunction(u::Vector{T},Alu,dfld::Array{T,3},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::grid) where T <: Real

    d = dfld[γ.wet] # couldn't use view b.c. of problem with function below
    list = surfaceindex(γ.I)
    gJ = Vector{T}(undef,size(u))
    [gJ[ii] = 2*(Q⁻*u[ii]) for ii in eachindex(list)]
    Jcontrol = u'*(Q⁻*u)

    # use in-place functions: more performant
    control2state!(d,u,γ) # d stores Δd
    ldiv!(Alu,d) # d stores c̃

    ỹ = state2obs(d,wis,γ)
    ỹ .-= y # stores n, data-model misfit

    gỹ = 2*(Wⁱ*ỹ)

    #gd = Array{T,3}(undef,size(dfld))
    #gd = Vector{T}(undef,sum(γ.wet))
    gd = zeros(T,sum(γ.wet))
    for ii in eachindex(y)
        # interpweights repeats some calculations
        gd .+= gỹ[ii] * interpweights(locs[ii],γ)[γ.wet]
    end
    # do Eᵀ gỹ 
    ldiv!(Alu',gd)
    list = surfaceindex(γ.I)
    #[gJ[ii] = gd[list[ii]] + 2*(Q⁻*u[ii]) for ii in eachindex(list)]
    [gJ[ii] += gd[list[ii]] for ii in eachindex(list)]

    J = ỹ'* (Wⁱ * ỹ) + Jcontrol
    return J, gJ
end

""" 
    function costfunction!(J,gJ,u,Alu,dfld,yfld,Wⁱ,wis,Q⁻,γ)
    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue: couldn't figure out how to nest with costfunction_obs!
# Arguments
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `dfld`: model constraints
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `Q⁻`: weights for control vector
- `γ`: grid
"""
function costfunction!(J,gJ,u::Vector{T},Alu,dfld::Array{T,3},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::grid) where T <: Real

    d = dfld[γ.wet] # couldn't use view b.c. of problem with function below

    if gJ != nothing
        #gJ = 2*(Q⁻*u)
    end
    if J != nothing

    end

    # use in-place functions: more performant
    control2state!(d,u,γ) # d stores Δd
    ldiv!(Alu,d) # d stores c̃

    ỹ = state2obs(d,wis,γ)
    ỹ .-= y # stores n, data-model misfit

    if gJ != nothing    
        gỹ = 2*(Wⁱ*ỹ)

        gd = zeros(T,sum(γ.wet))
        for ii in eachindex(y)
            # interpweights repeats some calculations
            gd .+= gỹ[ii] * interpweights(locs[ii],γ)[γ.wet]
        end
        # do Eᵀ gỹ 
        ldiv!(Alu',gd)
        list = surfaceindex(γ.I)

        [gJ[ii] = 2*(Q⁻*u[ii]) for ii in eachindex(list)]
        [gJ[ii] += gd[list[ii]] for ii in eachindex(list)]

    end

    if J != nothing
        Jcontrol = u'*(Q⁻*u)
        return  ỹ'* (Wⁱ * ỹ) + Jcontrol
    end
end

""" 
    function steady_inversion(u,Alu,d,γ.wet)
    invert for a steady-state tracer distribution
# Arguments
- `u`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `d`: model constraints
- `wet`: BitArray ocean mask
# Output
- `c`: steady-state tracer distribution
"""
function steady_inversion(uvec::Vector{T},Alu,d::Array{T,3},wet::BitArray{3}) where T <: Real
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    # below surface = 0 % no internal sinks or sources.
    u = tracerinit(wet[:,:,1],T)
    u[wet[:,:,1]] = uvec

    c = tracerinit(wet,T)
    n = tracerinit(wet,T)
    
    # first-guess reconstruction of observations
    Δd = d + control2state(u,wet)
    c[wet] =  Alu\Δd[wet]

    return c
end

"""
    function wetlocation(γ)
    Get (lon,lat,depth) tuples of wet locations.
    Allow a location to be wet if at least one out of 8 nearby gridpoints is wet.
    Certainly "wet" gridpoints could be defined more strictly.
# Arguments
- `γ`: TMI.grid
# Output
- `loc`: lon,lat,depth tuple
"""
function wetlocation(γ)

    confirmwet = false
    neighbors  = 8
    while !confirmwet
        loc = (rand(0.0:0.1:360.0),
               rand(-90.0:0.1:90.0),
               rand(0.0:1.0:5750.0))

        iswet(loc,γ) && return loc
        println("dry point, try again")
    end # if not, then start over.
end

function iswet(loc,γ,neighbors)
    # two approaches
    # approach 2
    # find 8 nearest neighbors
    Inn = nearestneighbor(loc,γ,neighbors)

    # are any of them wet?
    for ii = 1:neighbors
        if γ.wet[Inn[ii]]
            return true
        end
    end
    return false
end

function iswet(loc,γ)
    # two approaches
    # approach 1

    # 1 = very strict
    # 0 = all points
    wetness = 0.5
    
    wis = interpindex(loc,γ)

    # handle wraparound
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)

    # are any of them wet?
    # interpolate ones and zeros on to this loc.
    # if there is land nearby, the interpolated value
    # will be greater than 0.
    # this criterion only requires on land point nearby,
    # where nearby is one of the 8 corners of the cube that contains loc
    return wetwrap[wis...] > wetness
end

""" 
Save TMI configuration to NetCDF format for non-proprietary access
"""
function config2nc(TMIversion,A,γ,L,B)

    # make new netcdf file.
    filenetcdf = datadir("TMI_"*TMIversion*".nc")
    isfile(filenetcdf) && rm(filenetcdf)

    grid2nc(TMIversion,γ)
    
    matfields2nc(TMIversion,γ)

    watermassmatrix2nc(TMIversion,A)

    circulationmatrix2nc(TMIversion,L,γ)

    boundarymatrix2nc(TMIversion,B)
    
    regions2nc(TMIversion,γ)

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

    filenetcdf = datadir("TMI_"*TMIversion*".nc")
    filemat = datadir("TMI_"*TMIversion*".mat")
    vars = matread(filemat)

    TMIgrids, TMIgridsatts = griddicts(γ)

    T = eltype(γ.lon) # does the eltype of longitude have to equal the tracer eltype?
    #T =  Float64

    varlist = Dict("dP" => "qPO₄",
                   "Tobs" => "θ",
                   "Terr" => "σθ",
                   "Sobs" => "Sp",
                   "Serr" => "σSp",
                   "O18obs" => "δ¹⁸Ow",
                   "O18err" => "σδ¹⁸Ow",
                   "Pobs" => "PO₄",
                   "Perr" => "σPO₄",
                   "Nobs" => "NO₃",
                   "Nerr" => "σNO₃",
                   "Oobs" =>  "O₂",
                   "Oerr" =>  "σO₂",
                   "C13obs" =>  "δ¹³C",
                   "C13err" =>  "σδ¹³C")

    # iterate over all possible variables listed above
    Izyx = cartesianindex(filemat)
    TMIfields = Dict{String,Array{T,3}}()
    for (kk,vv) in varlist
        haskey(vars,kk) ? push!(TMIfields, vv => tracerinit(vars[kk], Izyx, γ.wet)) : nothing
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
         "Sp" => Dict("longname" => "practical salinity", "units" => "PSS-78"),
         "σSp" => Dict("longname" => "1σ standard error in practical salinity", "units" => "PSS-78"),
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
         "F₀" => Dict("longname" => "normalized mass flux out of gridcell", "units" => "(kg seawater/s)/(kg gridcell)"))


"""
Read vectors from mat file, translate to 3D,
 and save surface field to NetCDF file.
"""
function regions2nc(TMIversion,γ)

    filenetcdf = datadir("TMI_"*TMIversion*".nc")
    filemat = datadir("TMI_"*TMIversion*".mat")

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
    d_all = read(matobj,"d_all")
    close(matobj)

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

    filenetcdf = datadir("TMI_"*TMIversion*".nc")
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
Save circulation matrix `L` to NetCDF file.
"""
function circulationmatrix2nc(TMIversion,L,γ)

    T = eltype(L)
    fullmatrix = false # more efficient to just save F₀, then modify A to get L 
    filenetcdf = datadir("TMI_"*TMIversion*".nc")
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

    filenetcdf = datadir("TMI_"*TMIversion*".nc")
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

    filenetcdf = datadir("TMI_"*TMIversion*".nc")

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
Read an oceanographically-relevant surface region from NetCDF file. (Also could be read from mat file.)
"""
function surfaceregion(TMIversion,region,γ)

    file = datadir("TMI_"*TMIversion*".nc")
    T = eltype(γ.lon)
    tracername = "d_"*region
    dsfc = ncread(file,tracername)
    println("sumdsfc",sum(filter(!isnan,dsfc)))
    println("maxdsfc",maximum(filter(!isnan,dsfc)))
    # expand dsfc to cover 3D.
    # will not use control2state, because the control
    # make include non-surface regions in future.

    # preallocate
    d = Array{T}(undef,size(γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    d[γ.wet] .= zero(T)
    d[.!γ.wet] .= zero(T)/zero(T)
    d[:,:,1] = dsfc
    println(count(!isnan,d))
    return d
end

end
