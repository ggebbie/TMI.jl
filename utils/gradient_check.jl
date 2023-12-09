
# check with forward differences
ϵ = 1e-3
ii = rand(1:length(uvec))
println("Location for test =",ii)
δu = copy(uvec); δu[ii] += ϵ
∇f_finite = (f(δu) - f(uvec))/ϵ

fg!(J₀,gJ₀,(uvec+δu)./2) # J̃₀ is not overwritten
∇f = gJ₀[ii]
 
# error less than 10 percent?
println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
