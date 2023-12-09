## oxygen distribution, just be sure it runs
yO₂ = readfield(TMIfile,"O₂",γ)

bO₂ = getsurfaceboundary(yO₂)

O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)
