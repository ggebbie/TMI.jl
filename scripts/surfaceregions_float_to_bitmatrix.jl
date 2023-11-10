using Revise, TMI

TMIversion = "modern_90x45x33_GH10_GH12"
# version 2 (BitMatrix) file with regional masks
regionsfile = TMI.pkgdatadir("regions_90x45.nc")

A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

#region = "ANT"

for region in TMI.regionlist()
    # read v1 of regions from NetCDF file: used Floating point numbers for mask
    b = TMI.surfaceregion(TMIversion,region,γ,v1=true)

    # change b to a BitArray
    mask = (b.tracer .==1 .&& b.wet)

    @assert sum(mask) ≤ sum(b.wet)
    @assert sum(mask) == sum(b.tracer .==1)

    bnew = BoundaryCondition(mask,b.i,b.j,b.k,b.dim,b.dimval,b.wet,b.name,b.longname,"none") #where T <: Real
    TMI.write(TMI.pkgdatadir(regionsfile),bnew)
end
