# TMI solutions were originally saved in mat format. Here convert to NetCDF using 
using Revise, MAT, LinearAlgebra, SparseArrays, TMI, DrWatson, NCDatasets, Test

TMIversion = "modern_90x45x33_GH10_GH12"
TMIversion = "modern_180x90x33_GH10_GH12"

# original NetCDF version
#A, Alu, γ, TMIfile = config(TMIversion)

@time A, Alu, γ, TMIfile, L, B = config_from_mat(TMIversion);

config2nc(TMIversion,A,γ,L,B)

filenetcdf = datadir("TMI_"*TMIversion*".nc")
NCDataset(filenetcdf)

# read from NetCDF?
@time A2, Alu2, γ2, TMIfile2, L2, B2 = config_from_nc(TMIversion);

@test A2 == A
@test L2 == L
@test B2 == B
