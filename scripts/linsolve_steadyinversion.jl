import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using FiniteDiff
using Test
using Optim, LineSearches, BenchmarkTools
using LinearSolve
percent_difference(x, y) = @. 100 * (x - y) / y

TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

cobs = (θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
#     δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
#     P★ = preformedphosphate(TMIversion,Alu,γ),
    # δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)
c0 = cobs
m0 = massfractions_isotropic(γ)

u₀ = map(v -> getsurfaceboundary(v), cobs)
q₀ = map(v -> nothing, cobs)

function classic_steady_inversion(m, b, q, γ)
    A = watermassmatrix(m, γ)
    Alu = lu(A)
    c = steadyinversion(Alu, b, q, γ)
    return c
end

using LinearSolve


function linsolve_steady_inversion!(linsolve,b,γ::Grid{T};q=nothing,r=1.0) where T
    # preallocate Field for equation constraints
    d = zeros(γ,b.name,b.longname,b.units)
    
    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end
    linsolve.b .= vec(d)
    solve!(linsolve)
    c = zeros(d.γ,b.name,b.longname,b.units)
    c.tracer[wet(c)] .= linsolve.u
    return c
end
function linsolve_steady_inversion(m, b, q, γ)
    A = watermassmatrix(m, γ)
    x = zeros(size(A, 1))
    prob = LinearProblem(A, x)
    linsolve = init(prob);

    c_nt = map(b, q) do b_i, q_i
            linsolve_steady_inversion!(linsolve, b_i, γ; q = q_i, r = 1.0)
            end
    return c_nt
end
