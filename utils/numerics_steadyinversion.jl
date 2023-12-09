# get observations at surface
# set them as surface boundary condition
yPO₄ = readfield(TMIfile,"PO₄",γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.

# choice: BoundaryCondition or NamedTuple(BoundaryCondition)
#bPO₄ = getsurfaceboundary(yPO₄)
bPO₄ = getsurfaceboundary(yPO₄)

## preformed phosphate
PO₄pre = steadyinversion(Alu,bPO₄,γ)

## read phosphate source
qPO₄ = readsource(TMIfile,"qPO₄",γ)

# zero boundary condition, choose one line of next two
#b₀ = zerosurfaceboundary(γ)
b₀ = zerosurfaceboundary(γ)

# remineralized phosphate
PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)

# total (observed) phosphate
PO₄total = PO₄ᴿ + PO₄pre

## compute total phosphate directly
PO₄direct = steadyinversion(Alu,bPO₄,γ,q=qPO₄)
