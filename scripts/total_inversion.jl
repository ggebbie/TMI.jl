import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using FiniteDiff
using Test

percent_difference(x, y) = @. 100 * (x - y) / y

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
m0 = (west = m0[1], east = m0[2])
c = Field(tracer,
    γ,
    :c,
    "linear equilibrated tracer",
    "μmol/kg")

A = watermassmatrix(m0, γ)
Alu = lu(A)

dim = 1
b = (west = TMI.getboundarycondition(c, dim, 1, γ),
        east = TMI.getboundarycondition(c, 1, ngrid[dim], γ))

c̃ = steadyinversion(A,b,γ)
# Introduce a nonconservative variable with a source
qfield = 1.0e-2 * ones(ngrid)
# requires negative sign which is counterintutive (needs to be fixed)
q = TMI.Source(-qfield, γ, :q, "remineralized stuff", "μmol/kg", false)
c_noncons = steadyinversion(A,b,γ; q = q)
Δc = c_noncons - c


# Deep-copy priors so each control entry owns its storage.
u₀ = (c = deepcopy(b), c_q = deepcopy(b))
q₀ = (c = nothing, c_q = deepcopy(q),)
c0 = steadyinversion(Alu,u₀, q₀, γ)
Alu = lu(A)


# Controls and priors -------------------------------------------------------
ub = deepcopy(u₀)
uq = deepcopy(q₀)
Qᵤ = map(v -> Diagonal(one.(vec(v))), ub)
Qₛ = map(v -> isnothing(v) ? nothing : Diagonal(one.(vec(v))), uq)
Qₘ = Diagonal(one.(vec(m0)))

controls = ControlParameters(; u = ub, q = uq, m = m0,
                                u₀ = u₀, q₀ = q₀, m₀ = m0, 
                                Qᵤ = Qᵤ, Qₛ = Qₛ, Qₘ = Qₘ)

# Objective and gradient check for pointwise inversion ---------------------------------------------
obs_loc = [lon[2]]                                                           # single-point example
observed_vals = map(v -> observe(v, obs_loc, γ), c0)
W_pt = map(vals -> Diagonal(one.(vals)), observed_vals)
c_obs = map((vals, w) -> Observations(vals; locs = obs_loc, γ = γ, W = w), observed_vals, W_pt)

control_vector = randn(length(vec(controls)))
objective(x) = unconstrained_global_costfunction(x, controls, c_obs, γ)
Jg_finite = FiniteDiff.finite_difference_gradient(objective, control_vector, Val{:central})
J, Jg = unconstrained_global_costfunction(control_vector, controls, c_obs, γ; return_gradients = true)

@test all(abs.(percent_difference(Jg_finite, Jg)) .< 0.1) #all within 1 percent 

# Objective and gradient check for full field inversion ---------------------------------------------
W_full = map(v -> Diagonal(one.(vec(v))), c0)                                # full-field, unit weights
c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)                    # full-field obs

control_vector = randn(length(vec(controls)))

objective(x) = unconstrained_global_costfunction(x, controls, c_obs, γ)
Jg_finite = FiniteDiff.finite_difference_gradient(objective, control_vector, Val{:central})
J, Jg = unconstrained_global_costfunction(control_vector, controls, c_obs, γ; return_gradients = true)
@test all(abs.(percent_difference(Jg_finite, Jg)) .< 0.1) #all within 1 percent 

control_vector = randn(length(vec(controls)))
objective(x) = constrained_global_costfunction(x, controls, c_obs, γ; return_gradients = false)
Jg_finite = FiniteDiff.finite_difference_gradient(objective, control_vector, Val{:central})
J, Jg = constrained_global_costfunction(control_vector, controls, c_obs, γ; return_gradients = true)
@test all(abs.(percent_difference(Jg_finite, Jg)) .< 0.1) #all within 1 percent 



fg!(F, G, x) = optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=nothing)
x0 = vcat([vec(u₀), vec(q₀), randn(length(vec(m0)))]...)
lb = -Inf
ub = +Inf

_, _, xstart =unvec(controls, x0)


using Optim, LineSearches, BenchmarkTools
@btime result_opt_fg = Optim.optimize(
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

vec(q₀)
vec(qopt)


vec(u₀)
vec(uopt)
