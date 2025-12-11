import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using Test
# using Optim, LineSearches, BenchmarkTools
using LinearSolve, IncompleteLU
using NCDatasets
using SparseArrays
# Load a pre-configured TMI model and grid.
# This sets up the transport matrix `A`, its LU factorization `Alu`,
# the grid `γ`, and other model parameters.
TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

# Load observed tracer fields from the TMI NetCDF file.
# These will be the target for the inversion.
cobs = (
    θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
    # δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
    # PO₄ = readfield(TMIfile, "PO₄", γ),
    # NO₃= readfield(TMIfile, "NO₃", γ),
    # O₂= readfield(TMIfile, "O₂", γ),
    # δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

c0 = deepcopy(cobs)
# Define tracer standard deviations (σ) for boundary conditions and observations.
# These represent our confidence in the prior values.
tracer_error = (θ = 0.1, S = 0.01, 
                # δ¹⁸O = 0.2, 
                # PO₄ = 0.05, 
                # NO₃ = 1.0, O₂ = 5.0
                )

# Square them to get variances (σ²), which are used in the cost function.
tracer_error_variance = map(x -> x^2, tracer_error)

# Create observation objects with their corresponding precision (inverse variance).
# The precision matrix W is 1/σ², weighting the cost function.
observation_precision = map((v, σ²) -> Diagonal( (1/σ²) .* one.(vec(v)) ), c0, tracer_error_variance)
c_obs = map((v, w) -> Observations(v; W = w), c0, observation_precision)


c0 = deepcopy(cobs)

# Priors for surface boundary conditions (u₀) are taken from the observations.
u₀ = map(v -> getsurfaceboundary(v), cobs)
u₀_scale = map(v -> 1.0, u₀) #optional scaling to precondition Hessian (no scaling is 1.0)

# Priors for interior sources (q₀) are set to nothing (no sources).
q₀ = (
    θ =  nothing,
    S = nothing,
    # δ¹⁸O = nothing,
    # PO₄ = zerosource(γ),
    # NO₃= nothing,
    # O₂= nothing,
    # δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

q₀_scale =  map((v, u) -> (isnothing(v) ? nothing : deepcopy(u)), q₀, u₀_scale)
# source_error_variance = (PO₄ = 1.0, )

# q₀ = merge(!)
# Define optimization bounds for the boundary control variables.
u_lower = (θ = -2.0, S = 0.0, 
        #   δ¹⁸O = -Inf, 
        # PO₄ = 0.0, 
        # NO₃ = 0.0, O₂ = 0.0
          )
u_upper = (θ = 35.0, S = 45.0, 
        #    δ¹⁸O = Inf, 
        # PO₄ = 50.0, 
        # NO₃ = Inf, O₂ = Inf
           )

# q_lower = (PO₄ = 0,
        #   δ¹⁸O = -Inf, PO₄ = 0.0, NO₃ = 0.0, O₂ = 0.0
        #   )
# q_upper = (PO₄ = 35.0,
        #    δ¹⁸O = Inf, PO₄ = Inf, NO₃ = Inf, O₂ = Inf
        #    )

# Isotropic mass fractions as a prior estimate.
m0 = massfractions_isotropic(γ)
# m0 = inverse_watermassmatrix(A, γ)



# This bundles all priors, uncertainties, and bounds into a single object
# for the optimization routine, using the new `setup_inversion` API.
controls = TMI.setup_inversion(γ;
    boundary = (
        prior = u₀, 
        #variances should take the shape of the controls 
        variance = tracer_error_variance,
        lower_bound = descale_parameter(u_lower, u₀_scale),
        upper_bound = descale_parameter(u_upper, u₀_scale)
    ),
    source = (prior = q₀, 
    # variance = source_error_variance, 
    # lower_bound = descale_parameter(q_lower, q₀_scale),
    # upper_bound = descale_parameter(q_upper, q₀_scale)
    ),

    mass_fraction = (prior = m0, variance = 1.0, lower_bound = map(v -> 0.0, m0), upper_bound = map(v -> 1.1, m0))
);



function cache_f_and_g!(f_and_g!, dimension)
    last_x= nothing
    last_f, last_grad = 0.0, zeros(Float64, dimension)

    function update!(x)
        if x != last_x
            # println("Updating cache")
            last_f = f_and_g!(NaN, last_grad, x)
            last_x = copy(x)
        end
    end

    function f(x)
        update!(x)
        return last_f
    end

    function grad(x, G)
        update!(x)
        G .= last_grad
    end

    return f, grad
end

# The `fg!` function calculates the objective (F) and gradient (G) for Optim.
function fg_scaled!(F, G, x; u₀_scale, q₀_scale)
    optim_fg_unconstrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=nothing, 
                                                u₀_scale = u₀_scale, q₀_scale = q₀_scale)
end

x0_scaled = vcat([vec(descale_parameter(controls.boundary.u₀, u₀_scale)), 
                  vec(descale_parameter(controls.source.q₀, q₀_scale)), 
                  vec(controls.massfrac.m₀)]...)


f, grad = cache_f_and_g!((F, G, x) -> fg_scaled!(F, G, x; u₀_scale = u₀_scale, q₀_scale = q₀_scale), length(x0_scaled))


f(x0_scaled)
grad(x0_scaled, zero(x0_scaled))

f(x0_scaled .+ 1e-16) 
# idx = 489800
# ee = 1e-3
# x01 = deepcopy(x0_scaled); x01[idx] += ee
# x02 = deepcopy(x0_scaled); x02[idx] -= ee

unique(vec(controls.boundary.u₀) ./ vec(descale_parameter(controls.boundary.u₀, u₀_scale)))

using Ipopt
# Ipopt.PrintAvailableOptions()
n = length(x0_scaled)

# Define constraints (empty if only box constraints)
function eval_g!(x, g; m0)
    nm = length(vec(m0))
    m_x = x[end-nm+1:end]
    m = unvec(m0, m_x)

    firstm = first(m)
    γ = firstm.γ
    T = eltype(firstm.fraction)
    nint = sum(γ.interior)
    ngrid = size(γ.wet)

    total = zeros(T, ngrid)
    Iint = cartesianindex(γ.interior)
    for m1 in m
        sum(m1.γ.interior) == nint || error("mass fraction interior grids do not match")
    end
    for I in Iint
        s = zero(T)
        for m1 in m
            m1.γ.wet[I] || continue
            s += m1.fraction[I]
        end
        total[I] = s
    end

    total = total[firstm.γ.interior]
    for i in eachindex(g)
        g[i] = total[i]
    end
end

gtmp = zeros(sum(γ.interior))
eval_g!(x0_scaled, gtmp; m0)
gtmp
# Correct signature: (x, rows, cols, values)

function eval_jac_g(x, rows, cols, values; m0)
    nm     = length(vec(m0))
    offset = length(x) - nm              # start of mass-fraction block in x

    firstm    = first(m0)
    interior  = firstm.γ.interior
    row_of_I  = zeros(Int, size(interior))
    r = 0
    for I in cartesianindex(interior)     # map each interior cell -> constraint row
        interior[I] || continue
        r += 1
        row_of_I[I] = r
    end

    col = offset
    nz  = 1
    oneT = one(eltype(x))
    for m1 in m0
        wet = m1.γ.wet
        for I in eachindex(wet)          # follow vec(m0) ordering (wet cells only)
            wet[I] || continue
            col += 1
            row = row_of_I[I]
            row == 0 && continue          # not an interior constraint -> skip
            if values === nothing
                rows[nz] = row
                cols[nz] = col
            else
                values[nz] = oneT         # ∂g_row/∂m_cell = 1
            end
            nz += 1
        end
    end
    return
end


function eval_h(x, rows, cols, obj_factor, lambda, values)
    # Empty - using Hessian approximation
end

function jac_structure_stats(m0)
    interior = first(m0).γ.interior
    row_of_I = zeros(Int, size(interior))
    r = 0
    for I in cartesianindex(interior)
        interior[I] || continue
        r += 1
        row_of_I[I] = r
    end

    nnz = 0
    for m1 in m0
        wet = m1.γ.wet
        for I in eachindex(wet)
            (wet[I] && row_of_I[I] != 0) && (nnz += 1)
        end
    end
    return (nconstr = r, nnz = nnz)
end

vec(m0)

sum(TMI.sum_massfractions(m0).tracer[γ.interior])

TMI.sum_massfractions(m0).tracer[m0.east.γ.interior]

sum(m0.east.γ.interior)
sum(m0.west.γ.interior)
sum(γ.interior)

jac_stats = jac_structure_stats(m0)
nconstr = jac_stats.nconstr
g_L = g_U = ones(nconstr)
nnz_jac_g = jac_stats.nnz # number non-zeros in sparse

# quick check that the jacobian fills exactly nnz entries
rows_chk = zeros(Int, nnz_jac_g)
cols_chk = similar(rows_chk)
vals_chk = zeros(Float64, nnz_jac_g)
eval_jac_g(x0_scaled, rows_chk, cols_chk, vals_chk; m0 = m0)

# Create Ipopt problem
prob = Ipopt.CreateIpoptProblem(
    n,                              # number of variables
    Float64.(controls.lower_bound),           # lower bounds
    Float64.(controls.upper_bound),           # upper bounds
    nconstr,                              # number of constraints (0 if only box constraints)
    g_L,                      # constraint lower bounds = 1 for mass conservation
    g_U,                      # constraint upper bounds
    nnz_jac_g,                              # nnz in constraint Jacobian (0 if no constraints)
    0,                              # nnz in Hessian (0 = use approximation)
    f,
   (x,g) -> eval_g!(x, g; m0 = m0), #mass conservation constraint vector
    grad,
   (x, rows, cols, vals) -> eval_jac_g(x, rows, cols, vals; m0 = m0), #jacobian of cosntraint
    nothing
)

# Set Ipopt options
# Ipopt options
# Set Ipopt options
Ipopt.AddIpoptIntOption(prob, "print_level", 5)
Ipopt.AddIpoptStrOption(prob, "print_timing_statistics", "yes")
Ipopt.AddIpoptIntOption(prob, "max_iter", 5)
Ipopt.AddIpoptNumOption(prob, "acceptable_tol", 1e-6)
# # Barrier / Hessian
Ipopt.AddIpoptStrOption(prob, "mu_strategy", "adaptive")
Ipopt.AddIpoptStrOption(prob, "hessian_approximation", "limited-memory")
Ipopt.AddIpoptIntOption(prob, "limited_memory_max_history", 10)

# # Scaling (MOST IMPORTANT for your problem!)
Ipopt.AddIpoptStrOption(prob, "nlp_scaling_method", "gradient-based")

# Set starting point
prob.x = deepcopy(x0_scaled)

# Solve
Ipopt.AddIpoptStrOption(prob, "output_file", "./ipopt_trace.log")
status = Ipopt.IpoptSolve(prob)

C

# Extract solution
sol_x = prob.x
sol_obj = prob.obj_val
A, B, C= unvec(controls, sol_x)

extrema(vec(TMI.sum_massfractions(C)))

println("Status: ", status)
println("Objective value: ", sol_obj)
