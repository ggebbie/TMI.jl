import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using FiniteDiff
using Test
using BenchmarkTools
percent_difference(x, y) = @. 100 * (x - y) / y

TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

@btime lu(A)

findnz(Alu.L * Alu.U)

using SparseArrays

findnz(Alu)
# set them as surface boundary condition


cobs = (θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
#     δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
#     P★ = preformedphosphate(TMIversion,Alu,γ),
    # δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)
c0 = cobs

m0 = massfractions_isotropic(γ)
w = (θ =  0.01, S = 0.001,) #δ¹⁸O = 0.05, P★ = 0.05, δ¹³C★ = 0.05)

u₀ = map(v -> getsurfaceboundary(v), cobs)
q₀ = map(v -> nothing, cobs)

# Controls and priors -------------------------------------------------------
ub = deepcopy(u₀)
uq = deepcopy(q₀)

Qᵤ = map((v, wi) -> Diagonal(one.(vec(v)) .* 1/wi), ub,w )
Qₛ = map(v -> nothing, uq)
Qₘ = Diagonal(one.(vec(m0)))

controls = ControlParameters(; u = ub, q = uq, m = m0,
                                u₀ = u₀, q₀ = q₀, m₀ = m0, 
                                Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)

# Objective and gradient check for full field inversion ---------------------------------------------
W_full = map(v -> Diagonal(one.(vec(v))), c0)                                # full-field, unit weights
c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)                    # full-field obs


fg!(F, G, x) = optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=nothing)
x0 = vcat([vec(u₀), vec(q₀), randn(length(vec(m0)))]...)
lb = -Inf
ub = +Inf

using Optim, LineSearches, BenchmarkTools

fg!(NaN, zero.(x0), x0)

@btime fg!(NaN, zero.(x0), x0)


result_opt_fg = Optim.optimize(
    Optim.only_fg!(fg!),
    lb, ub, x0,
    Fminbox(LBFGS(;m = 5, 
            alphaguess = LineSearches.InitialHagerZhang(α0=NaN), 
            linesearch = LineSearches.HagerZhang())
            ),
    Optim.Options(f_abstol = 1e-12, g_tol = 1e-12,
                  iterations = 2500, store_trace = false,
                  show_trace = false, show_warnings = false));

uopt, qopt, xopt =unvec(controls, result_opt_fg.minimizer)
mopt = softmax_massfractions(xopt; α = 5.0)


@btime deepcopy(q₀)

vec(q₀)
vec(qopt)


vec(u₀)
vec(uopt)
