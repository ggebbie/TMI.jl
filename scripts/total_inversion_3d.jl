import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
# using FiniteDiff
using Test
using Optim, LineSearches, BenchmarkTools
using LinearSolve, IncompleteLU

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
w = (θ =  0.01, S = 0.001,) #δ¹⁸O = 0.05, P★ = 0.05, δ¹³C★ = 0.05)

u₀ = map(v -> getsurfaceboundary(v), cobs)
q₀ = map(v -> nothing, cobs)


# Controls and priors -------------------------------------------------------
ub = deepcopy(u₀)
uq = deepcopy(q₀)

Qᵤ = map((v, wi) -> Diagonal(one.(vec(v)) .* 1/wi), ub,w )
Qₛ = map(v -> nothing, uq)
Qₘ = Diagonal(one.(vec(m0)))

controls = ControlParameters(; γ = γ, ub = ub, uq = uq, m = m0,
                                u₀ = u₀, q₀ = q₀, m₀ = m0, 
                                Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)

# Objective and gradient check for full field inversion ---------------------------------------------
W_full = map(v -> Diagonal(one.(vec(v))), c0)                                # full-field, unit weights
c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)                    # full-field obs


function constrained_global_optimization(controls::ControlParameters, c_obs, γ::Grid; 
                                        x0::Union{Nothing,Vector{T}} = nothing, lowerb = -Inf, upperb = +Inf) where T

    function fg!(F, G, x)
        optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=nothing)
    end

    if isnothing(x0) #a naive first guess
        x0 = vcat([vec(controls.u₀), vec(controls.q₀), randn(length(vec(controls.m₀)))]...)
    end

    result_opt_fg = Optim.optimize(
        Optim.only_fg!(fg!), x0,
        LBFGS(;m = 10, 
                alphaguess = LineSearches.InitialHagerZhang(α0=NaN), 
                linesearch = LineSearches.HagerZhang()),
        Optim.Options(f_abstol = 1e-12, g_tol = 1e-8,
                    iterations = 5, store_trace = false,
                    show_trace = true, show_warnings = false));
    return result_opt_fg
end

x0 = vcat([vec(controls.u₀), vec(controls.q₀), randn(length(vec(controls.m₀)))]...)
# constrained_global_optimization(controls, c_obs, γ; x0 = x0)

# @btime constrained_global_optimization(controls, c_obs, γ; x0 = x0)