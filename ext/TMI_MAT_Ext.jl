module TMI_MAT_Ext

using TMI 
using MAT

import TMI: config, download_file, cartesianindex
import TMI: axislabels, watermassmatrix, boundarymatrix

export config
export download_file
export cartesianindex
export axislabels
export watermassmatrix
export boundarymatrix 

"""
Configure TMI environment from original MATLAB output
"""
function config(TMIversion, filetype)

    filetype ≠ "mat" && error("config: matlab files only")
    # repetitive and long
    TMIfile = download_file(TMIversion, filetype)
    
    # # make a sample field from zyx cartesian indices
    Izyx = cartesianindex(TMIfile, filetype)

    # # make a mask
    wet = falses(maximum(Izyx)[1],maximum(Izyx)[2],maximum(Izyx)[3])
    wet[Izyx] .= 1

    I = cartesianindex(wet)
    R = linearindex(wet)

    Azyx = watermassmatrix(TMIfile, filetype)

    A = matrix_zyx2xyz(TMIfile,Azyx,R)

    # # LU factorization for efficient matrix inversions
    Alu = lu(A)

    γ = Grid(TMIfile,A=A)

    # need to make this optional
    L = circulationmatrix(TMIfile, γ, filetype)
    
    B = boundarymatrix(TMIfile, γ, filetype)
    
    return  A, Alu, γ, TMIfile, L, B
end

"""
    function download_file(TMIversion::String, filetype)

    Download MATLAB file for given TMI version

# Arguments
- `TMIversion::String`

# Output
- `TMIfile::String`: TMI file name

# Side-effect
- download TMI input file if necessary
"""
function download_file(TMIversion::String, filetype)

    filetype ≠ "mat" && error("download_file: matlab files only")

    #- `url`: Google Drive URL for data
    url = maturl(TMIversion)
    TMIfile = pkgdatadir("TMI_"*TMIversion*".mat")
    TMIfilegz = TMIfile*".gz"
    println(TMIfile)
    !isdir(pkgdatadir()) && mkpath(pkgdatadir()) 

    if TMIversion == "modern_180x90x33_GH11_GH12"
        println("workaround for 2° x 2°")
        shellscript = pkgsrcdir("read_mat_modern_180x90x33_GH11_GH12.sh")
        run(`sh $shellscript`)
        mv(joinpath(pwd(),"TMI_"*TMIversion*".mat.gz"),TMIfilegz,force=true)
    elseif  TMIversion == "nordic_201x115x46_B23"
        println("workaround for regional Nordic Seas file")
        # warning: may not work due to changing Google API
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

function cartesianindex(file::String, filetype)

    filetype ≠ "mat" && error("cartesianindex: matlab files only")

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
    return CartesianIndex.(it,jt,kt) 
end

function axislabels(file::String, filetype)

    filetype ≠ "mat" && error("axislabels: matlab files only")

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
    return lon, lat, depth
end

function watermassmatrix(file, filetype)
    filetype ≠ "mat" && error("watermassmatrix: matlab files only")

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
    return A
end

function circulationmatrix(file, γ, filetype)
    filetype ≠ "mat" && error("circulatonmatrix: matlab files only")

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
        return L
    else
        close(matobj)
        return nothing
    end
end

function boundarymatrix(file, γ, filetype)

    filetype ≠ "mat" && error("boundarymatrix: matlab files only")
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
        return B
    else
        close(matobj)
        return nothing
    end
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
    Izyx = cartesianindex(file, "mat")
        
    # Julia accounting x,y,z
    ixyz = updatelinearindex(izyx,Izyx,R)
    jxyz = updatelinearindex(jzyx,Izyx,R)
    
    # use grid indices to switch i,j values
    Axyz = sparse(ixyz,jxyz,mzyx)
    return Axyz
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
    function mat2ncfield

    Rename MATLAB variables to NetCDF variables
"""
mat2ncfield() = Dict("CT"=>"Θ","Tobs"=>"θ","Tmod"=>"θ","Tlgm"=>"θ",
    "Terr"=>"σθ",
    "Sstar"=>"S⋆",
    "Sobs"=>"Sₚ","Smod"=>"Sₚ","Slgm"=>"Sₚ",
    "Sp" => "Sp", # help for some out-of-date input files
    "Serr"=>"σSₚ",
    "O18obs"=>"δ¹⁸Ow","O18mod"=>"δ¹⁸Ow","O18lgm"=>"δ¹⁸Ow",
    "O18err"=>"σδ¹⁸Ow",
    "Pobs"=>"PO₄","Pmod"=>"PO₄","Plgm"=>"PO₄","P"=>"PO₄",
    "Perr" => "σPO₄",
    "Nobs"=>"NO₃","Nmod"=>"NO₃","Nlgm"=>"NO₃","N"=>"NO₃",
    "Nerr" => "σNO₃",
    "Oobs"=>"O₂","Omod"=>"O₂","Olgm"=>"O₂","O"=>"O₂",
    "Oerr"=>"σO₂",
    "C13obs"=>"δ¹³C","C13mod"=>"δ¹³C","C13lgm"=>"δ¹³C",
    "C13err" =>  "σδ¹³C")

mat2ncsource() = Dict("dP"=>"qPO₄","q"=>"qPO₄")

function matvarnames(filemat)
    matobj = matopen(filemat)
    varnames = keys(matobj)
    xvarnames = nothing
    if haskey(matobj,"x")
        xvarnames = keys(read(matobj,"x"))
    end
    close(matobj)
    return varnames, xvarnames
end

"""
Read 3D fields from mat file and save to NetCDF file.
"""
function matfields2nc(TMIversion,γ)

    netcdffile = pkgdatadir("TMI_"*TMIversion*".nc")
    matfile = pkgdatadir("TMI_"*TMIversion*".mat")
    varnames, xvarnames = matvarnames(matfile)
    Izyx = cartesianindex(matfile, "mat") # use scope to call right cartesianindex?

    for (kk,vv) in mat2ncfield()
        if kk in varnames || kk in xvarnames #haskey(vars,kk)
            field = readfield(matfile,kk,γ,Izyx)
            writefield(netcdffile,field)
        end
    end
end

"""
Read 3D source field from mat file and save to NetCDF file.
"""
function matsource2nc(TMIversion,γ)

    netcdffile = pkgdatadir("TMI_"*TMIversion*".nc")
    matfile = pkgdatadir("TMI_"*TMIversion*".mat")
    varnames, xvarnames = matvarnames(matfile)
    Izyx = cartesianindex(matfile, "mat") # test this

    for (kk,vv) in mat2ncsource()
        if kk in varnames || kk in xvarnames 
            source = readsource(matfile,kk,γ,Izyx)
            writesource(netcdffile,source)
        end
    end
end

end # module 
