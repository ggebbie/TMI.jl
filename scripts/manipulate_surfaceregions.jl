using Revise, TMI

TMIversion = "modern_90x45x33_GH10_GH12"

regionsfile = "regions_90x45.nc"

A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

region = "ANT"
b = TMI.surfaceregion(TMIversion,region,γ)

# change b to a BitArray

mask = (b.tracer .==1 .&& b.wet)

sum(mask) ≤ sum(b.wet)
sum(mask) == sum(b.tracer .==1)

bnew2 = BoundaryCondition(mask,b.i,b.j,b.k,b.dim,b.dimval,b.wet,Symbol(region),String(region)*" mask","none") #where T <: Real

TMI.write(,bnew2)
