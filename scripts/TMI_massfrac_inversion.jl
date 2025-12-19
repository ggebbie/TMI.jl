import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using Test
using LinearSolve, IncompleteLU
using SparseArrays
using Ipopt


TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

# Load observed tracer fields from the TMI NetCDF file.
# These will be the target for the inversion.
c0 = (
    θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
    PO₄ = readfield(TMIfile, "PO₄", γ),
    NO₃= readfield(TMIfile, "NO₃", γ),
    O₂= readfield(TMIfile, "O₂", γ),
    # δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

vec(c0.θ)
vec(c0.θ.tracer[γ.interior])

filter(c0.θ)

cobs = deepcopy(c0)

# Define tracer standard deviations (σ) for observations.
tracer_error = (θ = 0.1, S = 0.01, 
                δ¹⁸O = 0.2, 
                PO₄ = 0.05, 
                NO₃ = 1.0, 
                O₂ = 5.0
                )

# Square them to get variances (σ²), which are used in the cost function.
tracer_error_variance = map(x -> x^2, tracer_error)

# Create observation objects with their corresponding precision (inverse variance).
# The precision matrix W is 1/σ², weighting the cost function.
observation_precision = map((v, σ²) -> Diagonal( (1/σ²) .* one.(vec(v)) ), c0, tracer_error_variance)
c_obs = map((v, w) -> Observations(v; W = w), c0, observation_precision)

# Priors for surface boundary conditions (u₀) are taken from the observations.
u₀ = map(v -> getsurfaceboundary(v), cobs)

u_lower = (θ = -2.0, S = 0.0, δ¹⁸O = -Inf, PO₄ = 0.0, NO₃ = 0.0, O₂ = 0.0)
u_upper = (θ = 35.0, S = 45.0, δ¹⁸O = Inf, PO₄ = 50.0, NO₃ = Inf, O₂ = Inf)

# This bundles all priors, uncertainties, and bounds into a single object
# for the optimization routine, using the new direct construction API.
boundary_controls = BoundaryControls(
    u₀;
);

# Priors for interior sources (q₀) are set to nothing (no sources).
q₀ = (
    θ =  nothing,
    S = nothing,
    δ¹⁸O = nothing,
    PO₄ = getsource(readfield(TMIfile, "qPO₄", γ)), #zerosource(γ),
    NO₃= nothing,
    O₂= nothing,
)
q₀.PO₄.tracer .*= -1 #sources are -1 for some reason? 
source_dependencies = (NO₃ = (:PO₄, 15.5), O₂ = (:PO₄, -170.0))

source_controls = SourceControls(
    q₀;
    dependencies = source_dependencies, 
);

# Isotropic mass fractions as a prior estimate.
m0 = massfractions_isotropic(γ)
# m0 = TMI.mass_frac_from_watermassatrix(A, γ)
# msol = inverse_watermassmatrix(A, γ)


massfrac_controls = MassFracControls(
    m0;
    variance = 1.0,
    lower_bound = map(v -> 0.0, m0),
    upper_bound = map(v -> 1.0, m0),
    γ = γ
);

controls = TMI.Controls(γ;
    boundary = boundary_controls,
    source = source_controls,
    massfrac = massfrac_controls
);

function fg!(F, G, x)
    joint_global_cost!(F, G, x, controls, c_obs, γ; locs=nothing, debug = true)
end

x0 = vec(controls)
ncontrols = length(x0)
f, grad = cache_f_and_g!((F, G, x) -> fg!(F, G, x), ncontrols)
f(x0)
g(x, grad::Function) = (g = zero(x); grad(x0, g); return g)
g(x) = g(x, grad)

ngrad_check = 5 #check 15 random points 
nrand_idx = idx = unique(rand(1:length(x0), ngrad_check))
perc_diff = zeros(ngrad_check)
@inbounds for(i, idx) in enumerate(nrand_idx)
    perc_diff[i] = check_gradient(x0, idx, f, g)
end

!all(perc_diff .< 0.1) && @error "Gradients are not consistent with finite differences!";

n_constraints = sum(γ.interior)  # number of interior points since only constraint is mass conservation
g_L = g_U = ones(n_constraints)  # constraint lower bounds both = 1 for mass conservation
n_nonzero_constraint_jacobian = length(vec(m0)) # number non-zeros in sparse

eval_h(x, rows, cols, obj_factor, lambda, values) = (nothing) # Empty since using a Hessian approximation
eval_c!(x, g) = mass_conservation_constraints!(x, g; m0 = m0) #mass conservation constraint
eval_jac_c(x, rows, cols, vals) = mass_conservation_jacobian!(x, rows, cols, vals; m0 = m0) #jacobian of constraints

# Create Ipopt problem
prob = Ipopt.CreateIpoptProblem(
    ncontrols,                              # number of variables
    controls.lower_bound,           # lower bounds
    controls.upper_bound,           # upper bounds
    n_constraints,                              # number of constraints (0 if only box constraints)
    g_L,                     
    g_U,                      # constraint upper bounds
    n_nonzero_constraint_jacobian,                              # nnz in constraint Jacobian (0 if no constraints)
    0,                              # nnz in Hessian (0 = use approximation)
    f,
    eval_c!, 
    grad,
    eval_jac_c, 
    nothing
);

# Some Ipopt solver options 
n_iterations = 5
Ipopt.AddIpoptIntOption(prob, "max_iter", n_iterations) #iterations 
Ipopt.AddIpoptStrOption(prob, "hessian_approximation", "limited-memory") #enable LBFGS
Ipopt.AddIpoptIntOption(prob, "acceptable_iter", 30)

Ipopt.AddIpoptIntOption(prob, "limited_memory_max_history", 27) #LBFGS Hessian approximation parameter 
Ipopt.AddIpoptNumOption(prob, "acceptable_tol", 1e-6)
Ipopt.AddIpoptStrOption(prob, "mu_strategy", "adaptive") 
Ipopt.AddIpoptStrOption(prob, "adaptive_mu_globalization", "kkt-error") 
Ipopt.AddIpoptStrOption(prob, "nlp_scaling_method", "gradient-based") #scale the gradients to avoid one dominating the problem

# Some Ipopt output options 
Ipopt.AddIpoptIntOption(prob, "print_level", 5); #verbose print to terminal
Ipopt.AddIpoptStrOption(prob, "print_timing_statistics", "yes")
Ipopt.AddIpoptStrOption(prob, "output_file", "./ipopt_trace.log") #output to a log file (Ipopt has no callback! :( )

prob.x = deepcopy(x0) # Set starting point
status = Ipopt.IpoptSolve(prob) # Solve
sol_x = prob.x

#QUICK ANALYSIS
b_sol, q_sol, m_sol= unvec(controls, sol_x)

m_true = mass_frac_from_watermassatrix(A, γ)
Alu_true = Alu
Alu_0 = lu(watermassmatrix(m0, γ))
Alu_sol= lu(watermassmatrix(m_sol, γ))

get_vfilled(Alu) = volumefilled(NaN, Alu, γ)
vfilled_true = get_vfilled(Alu_true)
vfilled_0 = get_vfilled(Alu_0) 
vfilled_sol = get_vfilled(Alu_sol)

vfill_min, vfill_max = extrema(filter(!isnan, vfilled_true.tracer))

using PythonPlot
fig, ax = subplots(1, 3, figsize = (15, 5), sharey = true)
ax[0].pcolormesh(γ.lon, γ.lat, vfilled_true.tracer', vmin = vfill_min, vmax = vfill_max)
ax[1].pcolormesh(γ.lon, γ.lat, vfilled_0.tracer', vmin = vfill_min, vmax = vfill_max)
ax[2].pcolormesh(γ.lon, γ.lat, vfilled_sol.tracer', vmin = vfill_min, vmax = vfill_max)
fig


get_vfilled_distribution(vfilled) = sort(filter(!isnan, vfilled.tracer))

fig, axes = subplots()
axes.plot(get_vfilled_distribution(vfilled_true), label = "true volume filled")
axes.plot(get_vfilled_distribution(vfilled_0), label = "isotropic volume filled")
axes.plot(get_vfilled_distribution(vfilled_sol), label = "optimized volume filled ($n_iterations Iterations)")
axes.legend()
fig