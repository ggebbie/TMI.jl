# version 2 (BitMatrix) file with regional masks
using Revise, TMI

for TMIversion in ["modern_90x45x33_GH10_GH12","modern_180x90x33_GH11_GH12"]

    TMIversion[8:9] == "90" ? gridsize = "90x45" : gridsize = "180x90"
    regionsfile = TMI.pkgdatadir("regions_"*gridsize*".nc")

    println(regionsfile)
    
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
end

# open a file and check contents.
ds = NCDataset(TMI.pkgdatadir("regions_90x45.nc")) # looks ok although ncview has non-zero and non-one artifacts
ds = NCDataset(TMI.pkgdatadir("regions_180x90.nc")) # looks ok although ncview has non-zero and non-one artifacts

btest = TMI.surfaceregion(TMIversion,:ANT,γ)

