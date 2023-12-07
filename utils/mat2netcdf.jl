# TMI solutions were originally saved in mat format. Here convert to NetCDF using 
using Revise, MAT, LinearAlgebra, SparseArrays, TMI, DrWatson, NCDatasets, Test

# choose one version
TMIversion = "modern_90x45x33_GH10_GH12"
TMIversion = "modern_180x90x33_GH11_GH12"
TMIversion = "modern_90x45x33_unpub12"
TMIversion = "modern_90x45x33_G14_v2"
TMIversion = "modern_90x45x33_G14"
TMIversion = "LGM_90x45x33_G14"
TMIversion = "LGM_90x45x33_G14A"
TMIversion = "LGM_90x45x33_GPLS1"
TMIversion = "LGM_90x45x33_GPLS2"
TMIversion = "LGM_90x45x33_OG18"
TMIversion = "nordic_201x115x46_B23"

# original NetCDF version
#A, Alu, γ, TMIfile = config(TMIversion)

@time A, Alu, γ, TMIfile, L, B = config_from_mat(TMIversion);
#@time A2, Alu2, γ2, TMIfile2, L2, B2 = config_from_mat(TMIversion2);

config2nc(TMIversion,A,γ,L,B)

filenetcdf = datadir("TMI_"*TMIversion*".nc")
NCDataset(filenetcdf)

# read from NetCDF?
@time A2, Alu2, γ2, TMIfile2, L2, B2 = config_from_nc(TMIversion);

# Do NetCDF and mat-files agree?
@test A2 == A
@test L2 == L
@test B2 == B
