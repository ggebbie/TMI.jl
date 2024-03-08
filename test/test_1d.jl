#using Revise
#using TMI
ngrid = (50) # number of grid cells
xmax = 1000.0 # domain size 
lon = collect(range(0.0,1000.0,length=ngrid[1]))
tracer = collect(1.0.-lon./xmax)

axes = (lon,)
wet = trues(ngrid)
interior = copy(wet)
interior[begin] = false
interior[end] = false

wrap = (false,)
Δ = [CartesianIndex(1,),CartesianIndex(-1,)]
γ = Grid(axes,wet,interior,wrap,Δ)
n = neighbors(γ)
m0 = massfractions_isotropic(γ)
c = Field(tracer,
    γ,
    :c,
    "linear equilibrated tracer",
    "μmol/kg")

y = (c = c,)
w = (c = 0.01,)

m = massfractions(y, w)
