# TMI solutions were originally saved in mat format. Here convert to NetCDF using 
using Revise, MAT, LinearAlgebra, SparseArrays, TMI, DrWatson

TMIversion = "modern_4x4x33_GH10_GH12"

A, Alu, γ, TMIfile = config(TMIversion)

TMIfile = config_from_mat(TMIversion)

on download(url,inputdir)
 
# load circulation file
Afile = "/home/gebbie/Dropbox/mcode/TMI_v8/A_4deg_2010.mat"
Avars = matread(Afile)

# TMI_modern_4x4x33_GH2010_GH2012.mat
https://drive.google.com/file/d/15O3QBzeSUIsuF7lZMbUixzZuEwDiaAdo/view?usp=sharing


tracerfile = "/home/gebbie/Dropbox/mcode/TMI_v8/tracerobs_4deg_33lev_TMI.mat"
tracervars = matread(tracerfile)

it = convert(Array{Int,1},vec(Avars["i"]));
jt = convert(Array{Int,1},vec(Avars["j"]));
kt = convert(Array{Int,1},vec(Avars["k"]));
Itmi = CartesianIndex.(it,jt,kt)

function initfld(vector::Array{Float64,2},I::Array{CartesianIndex{3},1})
# function for changing TMI vector to field

    #convert2vec!(vector)
    vector = convert(Array{Float64,1},vec(vector))
    println(size(vector))
    nv = length(vector)

    nx = Tuple(I[1])[1]
    ny = Tuple(I[1])[2]
    nz = Tuple(I[1])[3]
    
    for nn = 2:nv
        nx = max(nx,Tuple(I[nn])[1])
        ny = max(ny,Tuple(I[nn])[2])
        nz = max(nz,Tuple(I[nn])[3])
    end

    field = NaN.*zeros(nx,ny,nz)
    #- a comprehension
    [field[I[n]]=vector[n] for n ∈ 1:nv]

    return field
end

function convert2vec!(vector::Array{Float64,2})
    vector = convert(Array{Float64,1},vec(vector))
    return vector
end

function convert2intvec!(vector::Array{Float64,2})
    convert(Array{Int64,1},vec(vector))
end
function convert2int!(vector::Array{Float64,2})
    convert(Array{Int64,2},vector)
end

# update names and types in dictionary
TMIgrids = Dict("lon" => convert2intvec!(Avars["LON"]),
                "lat" => convert2intvec!(Avars["LAT"]),
                "depth" => convert2intvec!(Avars["DEPTH"]))
#                "Itmi" => Itmi,
#                 "A" => Avars["A"])
                 
TMIfields = Dict("ΔPO₄" => initfld(Avars["dP"],Itmi),
                "θ" => initfld(tracervars["Tobs"],Itmi),
                "σθ" => initfld(tracervars["Terr"],Itmi),
                "Sp" => initfld(tracervars["Sobs"],Itmi),
                "σSp" => initfld(tracervars["Serr"],Itmi),
                "δ¹⁸Ow" => initfld(tracervars["O18obs"],Itmi),
                "σδ¹⁸Ow" => initfld(tracervars["O18err"],Itmi),
                "PO₄" => initfld(tracervars["Pobs"],Itmi),
                "σPO₄" => initfld(tracervars["Perr"],Itmi),
                "NO₃" => initfld(tracervars["Nobs"],Itmi),
                "σNO₃" => initfld(tracervars["Nerr"],Itmi),
                "O₂" => initfld(tracervars["Oobs"],Itmi),
                "σO₂" => initfld(tracervars["Oerr"],Itmi),
                "δ¹³C" => initfld(tracervars["C13obs"],Itmi),
                "σδ¹³C" => initfld(tracervars["C13err"],Itmi))

## write to NetCDF
filenetcdf = "TMI_4deg_2010.nc"

TMIgridsatts = Dict("lon" => Dict("longname" => "Longitude", "units" => "°E"),
                    "lat" => Dict("longname" => "Latitude", "units" => "°N"),
                    "depth" => Dict("longname" => "depth", "units" => "m"))

TMIfieldsatts = Dict("θ" => Dict("longname" => "potential temperature", "units" => "°C"),
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
                     "σδ¹³C" => Dict("longname" => "1σ standard error fin carbon-13 to carbon-12 ratio in DIC", "units" => "‰ PDB"))
            

# make new netcdf file.
isfile(filenetcdf) && rm(filenetcdf)

# iterate in TMIgrids Dictionary to write to NetCDF.
for (varname,varvals) in TMIfields
  println(varname)
  nccreate(filenetcdf,varname,"lon",lon,TMIgridsatts["lon"],"lat",lat,TMIgridsatts["lat"],"depth",depth,TMIgridsatts["depth"],atts=TMIfieldsatts[varname])
  ncwrite(varvals, filenetcdf,varname)
end

isrc = A_CSC[2]
idst = A_CSC[1]
ival = A_CSC[3]

# add the circulation matrix: problem can't store sparse matrix.
varname= "m"
faceatts = Dict("longname" => "TMI grid face number")
matts =  Dict("longname" => "TMI water-mass (sparse) matrix values","version"=> "2010")
nccreate(filenetcdf,varname,"face",1:nface,faceatts,atts=matts)
ncwrite(ival, filenetcdf,varname)

varname= "i"
destatts = Dict("longname" => "gridcell number of destination (row value)")
nccreate(filenetcdf,varname,"face",1:nface,faceatts,atts=destatts)
ncwrite(idst, filenetcdf,varname)

varname= "j"
sourceatts = Dict("longname" => "gridcell number of source (column value)")
nccreate(filenetcdf,varname,"face",1:nface,faceatts,atts=sourceatts)
ncwrite(isrc, filenetcdf,varname)

varname = "xgrid"
xgridatts = Dict("longname" => "gridcell Cartesian x-location")
traceratts = Dict("longname" => "gridcell index")
nccreate(filenetcdf,varname,"index number",1:nv,traceratts,atts=xgridatts)
ncwrite(it, filenetcdf,varname)

varname = "ygrid"
ygridatts = Dict("longname" => "gridcell Cartesian y-location")
traceratts = Dict("longname" => "gridcell index")
nccreate(filenetcdf,varname,"index number",1:nv,traceratts,atts=ygridatts)
ncwrite(jt, filenetcdf,varname)

varname = "zgrid"
zgridatts = Dict("longname" => "gridcell Cartesian z-location")
traceratts = Dict("longname" => "gridcell index")
nccreate(filenetcdf,varname,"index number",1:nv,traceratts,atts=zgridatts)
ncwrite(kt, filenetcdf,varname)


##############################################################
## end of 4-jan-2021 work
##############################################################


function extract(d)
           expr = quote end
           for (k, v) in d
               push!(expr.args, :($(Symbol(k)) = $v))
           end
           eval(expr)
           return
       end

extract(Dict_TMI)

 dry = isnan.(po4)
 wet = .!dry

# get po4fld
 @time po4fld = NaN.*zeros(nx,ny,nz);
 @time tmivec2fld!(po4fld,po4,indextmi);

# get po4vector: slower process.
 @time po4vec = NaN.*zeros(nv);
 @time tmifld2vec!(po4vec,po4fld,indextmi);


function tmivec2fld!(field::Array{Real,3},vector::Array{Real,1},index::Array{CartesianIndex{3},1})

    nv = length(index)
    #- a comprehension
    [field[index[n]]=vector[n] for n ∈ 1:nv]

    return field
end

function tmifld2vec!(vector::Array{Real,1},field::Array{Real,3},index::Array{CartesianIndex{3},1})

    nv = length(index)
    
    #- a comprehension
    [vector[n] = field[indextmi[n]] for n ∈ 1:nv];
    return vector
end


Alu = lu(A)
nfield = size(A,1)



# vector to field via a comprehension (slow: 0.1 s)
# @time po4fld_test = NaN.*zeros(nx,ny,nz);
# @time [po4fld_test[indextmi[n]]=po4[n] for n ∈ 1:nv];

# field to vector via a comprehension
# @time po4 = [po4fld[indextmi[n]] for n ∈ 1:nv];
#- also works -# @time [po4[n] = po4fld[indextmi[n]] for n ∈ 1:nv];

#file = matopen(Afile)
#for (key, values) in Avars
#    println(key); #print(value)
    #stringer = $key=value
#    eval($key=value)
#end


