#=
% Create a TMI.Grid struct from foreign model output
=#

import Pkg; Pkg.activate("../scripts") # a more full-featured environment

using Revise
using TMI
using NCDatasets
using OrderedCollections

gamma =  Grid(TMI.pkgdatadir("simple_MITgcm.nc"), "THETA", "XC", "YC", "Z", "maskC") 

# scratch work: shows how the Grid constructor works
mitgcm_grid = Dataset(TMI.pkgdatadir("mitgcm_grid.nc"))
simple_mitgcm = Dataset(TMI.pkgdatadir("simple_MITgcm.nc"))

# for testing purposes, make a minimal NetCDF file
# ncgen(TMI.pkgdatadir("simple_MITgcm.nc")) -- gives helpful hints
ds = NCDataset("MITgcm_test.nc","c")
ds.dim["Z"] = 32
ds.dim["YC"] = 128
ds.dim["XC"] = 64

ncTHETA = defVar(ds,"THETA", Float32, ("XC", "YC", "Z"), attrib = OrderedDict(
    "_FillValue"                => Float32(NaN),
    "standard_name"             => "THETA",
    "long_name"                 => "Potential Temperature",
    "units"                     => "degC",
    "coordinates"               => "Depth PHrefC drF dxF dyF hFacC iter maskC rA rhoRef",
))

ncXC = defVar(ds,"XC", Float32, ("XC",), attrib = OrderedDict(
    "_FillValue"                => Float32(NaN),
    "standard_name"             => "longitude",
    "long_name"                 => "longitude",
    "units"                     => "degrees_east",
    "coordinate"                => "YC XC",
    "axis"                      => "X",
))

ncYC = defVar(ds,"YC", Float32, ("YC",), attrib = OrderedDict(
    "_FillValue"                => Float32(NaN),
    "standard_name"             => "latitude",
    "long_name"                 => "latitude",
    "units"                     => "degrees_north",
    "coordinate"                => "YC XC",
    "axis"                      => "Y",
))

ncZ = defVar(ds,"Z", Float32, ("Z",), attrib = OrderedDict(
    "_FillValue"                => Float32(NaN),
    "standard_name"             => "depth",
    "long_name"                 => "vertical coordinate of cell center",
    "units"                     => "m",
    "positive"                  => "down",
    "axis"                      => "Z",
))

ncmaskC = defVar(ds,"maskC", Int8, ("XC", "YC", "Z"), attrib = OrderedDict(
    "standard_name"             => "sea_binary_mask_at_t_location",
    "long_name"                 => "mask denoting wet point at center",
    "dtype"                     => "bool",
))

ncTHETA[:] = simple_mitgcm["THETA"][:,:,:,begin]
ncmaskC[:] = simple_mitgcm["maskC"]
ncXC[:] = simple_mitgcm["XC"]
ncYC[:] = simple_mitgcm["YC"]
ncZ[:] = simple_mitgcm["Z"]
close(ds)

# use potential temperature as a template
# should double-check consistency with other tracers (but is not done here)
# problem: Union Missing type
c = convert(Array{Float32,3},simple_mitgcm["THETA"][:,:,:,begin]) # take one snapshot
lon = convert(Vector{Float32},simple_mitgcm["XC"]) 
lat = convert(Vector{Float32},simple_mitgcm["YC"])
depth = - convert(Vector{Float32},simple_mitgcm["Z"]) # flip sign for actual "depth"
wet = Bool.(simple_mitgcm["maskC"])::BitArray # very slow! (couple of secs), use `convert` instead?
