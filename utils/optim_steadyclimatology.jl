F,G = costfunction_gridded_obs(uvec,Alu,b,u,y,W⁻,γ)
fg!(F,G,x) = costfunction_gridded_obs!(F,G,x,Alu,b,u,y,W⁻,γ)
fg(x) = costfunction_gridded_obs(x,Alu,b,u,y,W⁻,γ)
f(x) = fg(x)[1]
J₀,gJ₀ = fg(uvec)
iterations = 10
out = steadyclimatology(uvec,fg!,iterations)
