#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

copy-paste from ex8
1d copy-paste from TVCR#bh-1dA/1dA_separate.jl 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#
import Pkg; Pkg.activate(".")

#using Revise
using LinearAlgebra
using TMI
#using GeoPythonPlot # will load optional extension
using PyPlot
using OrdinaryDiffEq
close("all") 
#using Unitful 
#TMIversion = versionlist()[6] # G14 has no remote mass fractions
#A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

umin = 0.1 
umax = 1000
κ = 10000
Δx = 10 
discrete_x = 0:Δx:2000
N = length(discrete_x) 
function makeA_periodic(u, κ, N) 
    #A = Matrix{Any}(undef, length(discrete_x) + 2, length(discrete_x) + 2)
    #A .= 0 * s^-1
    w₁ = u / Δx + κ / Δx^2
    w₂ = -(u / Δx + κ / Δx^2 + κ / Δx^2 )
    w₃ = κ / Δx^2
    display(fill(w₂, N + 2))
    A = diagm(fill(w₂, N + 2))
    [A[i, i+1] = w₃ for i in 1:N + 1]
    [A[i, i-1] = w₁ for i in 2:N + 2]
    A[1, end - 1] = w₁
    A[end, 2] = w₃
    A = A
    return A
end

A05 = makeA_periodic(umax, 0, N) # u = 5cm/s
Aκ = makeA_periodic(0, κ, N) #decrease κ to be more observable 
Acombined = A05 + Aκ
@show (A05 + Aκ) == makeA_periodic(umax, κ, N) #nice

c₀ = vcat(fill(0, 101), fill(1,10), fill(0,92))
function func!(du, u, p, t)
    du .= p * u
end
figure()
for (i, (A, tit)) in enumerate(zip((A05, Aκ, Acombined), L"\frac{\partial c}{\partial t} = " .*["Aᵘ", L"A^\kappa", L"(A^u + A^\kappa)"] .* "c"))
    subplot(3,1,i)
    prob = ODEProblem(func!, c₀, (0.0, 1.0),A, saveat = 0.0:0.1:1.0)
    sol = solve(prob, Tsit5())
    cmap = get_cmap()
    [plot(discrete_x, sol.u[i][2:end-1], color = cmap(i / length(sol.t))) for i in 1:length(sol.u)]
    ylabel("c [unitless]")
    title(tit)
end
xlabel("Distance [km]")
tight_layout()
#savefig(plotsdir("1dAs_split.png"))


# get observations at surface
# set them as surface boundary condition
y = (θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
    P★ = preformedphosphate(TMIversion,Alu,γ),
    δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

@time m̃ = massfractions(y);
Ã = watermassmatrix(m̃, γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(y.θ)

## reconstruct temperature
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ã,bθ,γ)

# compare to c.θ
Base.maximum(y.θ - θ̃)
Base.minimum(y.θ - θ̃)

