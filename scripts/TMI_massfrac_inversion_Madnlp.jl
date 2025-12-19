import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using Test
using SparseArrays

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
)

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

source_dependencies = (NO₃ = (:PO₄, 15.5),O₂ = (:PO₄, -170.0))

q_lower = (PO₄ = 0.0, )
q_upper = (PO₄ = 0.5, )

source_controls = SourceControls(
    q₀;
    dependencies = source_dependencies, 
);

# Isotropic mass fractions as a prior estimate.
m0 = massfractions_isotropic(γ)

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

# This is the core objective and gradient function.
# It must return the objective value and modify the gradient `G` in-place.
function fg!(F_val_dummy, G, x)
    objective = joint_global_cost!(F_val_dummy, G, x, controls, c_obs, γ; locs=nothing, debug=false)
    return objective
end

x0 = vec(controls)
ncontrols = length(x0)
f, grad = cache_f_and_g!((F, G, x) -> fg!(F, G, x), ncontrols)
f(x0)
g(x, grad::Function) = (g = zero(x); grad(x0, g); return g)
g(x) = g(x, grad)

ngrad_check = 5 #check a few random points 
nrand_idx = idx = unique(rand(1:length(x0), ngrad_check))
perc_diff = zeros(ngrad_check)
@inbounds for(i, idx) in enumerate(nrand_idx)
    perc_diff[i] = check_gradient(x0, idx, f, g)
end

!all(perc_diff .< 0.1) && @error "Gradients are not consistent with finite differences!"

# ----------------------------------------------------
# --- MadNLP via NLPModels.jl wrapper definition ---
# ----------------------------------------------------
using NLPModels
using MadNLP

@inline function _same_x(a::AbstractVector, b::AbstractVector)
    @inbounds for i in eachindex(a, b)
        if a[i] != b[i]
            return false
        end
    end
    return true
end

mutable struct TMIMadModel <: NLPModels.AbstractNLPModel{Float64, Vector{Float64}}
    meta::NLPModels.NLPModelMeta{Float64, Vector{Float64}}
    counters::NLPModels.Counters
    fg!::Function
    cons!_impl::Function
    jac!_impl::Function
    jac_rows::Vector{Int}
    jac_cols::Vector{Int}
    xcache::Vector{Float64}
    fcache::Float64
    gcache::Vector{Float64}
    cache_valid::Bool
end

function TMIMadModel(x0::Vector{Float64}, lvar::Vector{Float64}, uvar::Vector{Float64},
                     lcon::Vector{Float64}, ucon::Vector{Float64}, nnzj::Int;
                     fg!::Function, cons!_impl::Function, jac!_impl::Function)

    jac_rows = Vector{Int}(undef, nnzj)
    jac_cols = Vector{Int}(undef, nnzj)
    jac!_impl(x0, jac_rows, jac_cols, nothing)

    # The NLPModelMeta constructor is sensitive to the version of NLPModels.jl.
    # By omitting the `lin` keyword argument, all constraints are treated as nonlinear.
    # This avoids a call to the unimplemented `jtprod_lin!` function.
    meta = NLPModels.NLPModelMeta(
        length(x0),
        x0 = x0,
        lvar = lvar,
        uvar = uvar,
        ncon = length(lcon),
        lcon = lcon,
        ucon = ucon,
        nnzj = nnzj,
        nnzh = 0,
        hess_available = false,
        minimize = true
    )

    return TMIMadModel(
        meta,
        NLPModels.Counters(),
        fg!,
        cons!_impl,
        jac!_impl,
        jac_rows,
        jac_cols,
        copy(x0),
        0.0,
        zeros(length(x0)),
        false
    )
end

function _update_cache!(nlp::TMIMadModel, x::AbstractVector{<:Real})
    if nlp.cache_valid && _same_x(nlp.xcache, x)
        return
    end
    nlp.fcache = nlp.fg!(NaN, nlp.gcache, x)
    copyto!(nlp.xcache, x)
    nlp.cache_valid = true
    return
end

# --- NLPModels API Implementation ---
# The :neval_* symbols are the modern API for incrementing counters.
function NLPModels.obj(nlp::TMIMadModel, x::AbstractVector)
    NLPModels.increment!(nlp, :neval_obj)
    _update_cache!(nlp, x)
    return nlp.fcache
end

function NLPModels.grad!(nlp::TMIMadModel, x::AbstractVector, g::AbstractVector)
    NLPModels.increment!(nlp, :neval_grad)
    _update_cache!(nlp, x)
    copyto!(g, nlp.gcache)
    return g
end

function NLPModels.objgrad!(nlp::TMIMadModel, x::AbstractVector, g::AbstractVector)
    NLPModels.increment!(nlp, :neval_obj)
    NLPModels.increment!(nlp, :neval_grad)
    _update_cache!(nlp, x)
    copyto!(g, nlp.gcache)
    return nlp.fcache, g
end

function NLPModels.cons!(nlp::TMIMadModel, x::AbstractVector, c::AbstractVector)
    NLPModels.increment!(nlp, :neval_cons)
    nlp.cons!_impl(x, c)
    return c
end

function NLPModels.jac_structure!(nlp::TMIMadModel, I::AbstractVector{T}, J::AbstractVector{T}) where {T}
    copyto!(I, nlp.jac_rows)
    copyto!(J, nlp.jac_cols)
    return
end

function NLPModels.jac_coord!(nlp::TMIMadModel, x::AbstractVector, vals::AbstractVector)
    NLPModels.increment!(nlp, :neval_jac)
    nlp.jac!_impl(x, nlp.jac_rows, nlp.jac_cols, vals)
    return vals
end

# This function for the Jacobian-transpose product of nonlinear constraints
# is required by the API. Since we are not flagging constraints as linear, 
# this function is responsible for the Jacobian-transpose-vector product.
import NLPModels: jtprod_nln!
function NLPModels.jtprod_nln!(nlp::TMIMadModel, x::AbstractVector, v::AbstractVector, Jtv::AbstractVector)
    # The parent function `jtprod!` may not zero the vector first, so we do it here for safety.
    fill!(Jtv, 0.0)

    # The Jacobian of the constraints is constant, so `x` is not used.
    # We compute J' * v.
    # The Jacobian values are all 1.0, and the structure is stored in nlp.jac_rows and nlp.jac_cols.
    for k in 1:length(nlp.jac_rows)
        i = nlp.jac_rows[k]
        j = nlp.jac_cols[k]
        # The Jacobian values are implicitly 1.0 for all specified entries.
        Jtv[j] += v[i]
    end

    return Jtv
end

# ---------------------------
# --- Problem Definition ----
# ---------------------------
n_constraints = sum(γ.interior)
lcon = ones(n_constraints)
ucon = ones(n_constraints)
nnzj = length(vec(m0))

cons!_impl(x, c) = mass_conservation_constraints!(x, c; m0 = m0)
jac!_impl(x, rows, cols, vals) = mass_conservation_jacobian!(x, rows, cols, vals; m0 = m0)

nlp = TMIMadModel(
    x0,
    controls.lower_bound,
    controls.upper_bound,
    lcon,
    ucon,
    nnzj;
    fg! = fg!,
    cons!_impl = cons!_impl,
    jac!_impl = jac!_impl
)

# -------------------
# --- Run Solver ----
# -------------------
n_iterations = 25
using MadNLPMumps
solver = MadNLP.MadNLPSolver(
    nlp;
    # termination
    max_iter = n_iterations,
    acceptable_tol = 1e-6,

    # quasi-Newton (L-BFGS)
    hessian_approximation = MadNLP.CompactLBFGS,
    quasi_newton_options = MadNLP.QuasiNewtonOptions(; max_history = 15),

    # Jacobian is linear/constant here
    jacobian_constant = true,
    blas_num_threads = 4,
    # output
    print_level = MadNLP.INFO,
    output_file = "./madnlp_trace.log",
    file_print_level = MadNLP.TRACE,
    # linear solver
    linear_solver = MumpsSolver,
)

# MadNLP.MadNLPMumps
results = MadNLP.solve!(solver)
sol_x = results.solution

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

# if results.status == MadNLP.SOLVE_SUCCEEDED
#     println("Solve succeeded!")
#     
# else
#     println("Solve failed with status: ", results.status)
# end