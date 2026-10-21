# METIS ordering: about 40% less fill than UMFPACK's default AMD ordering on TMI
# grids; `lu!` keeps the ordering of the factorization it refactors
const METIS_CONTROL = let control = SparseArrays.UMFPACK.get_umfpack_control(Float64, Int64)
    control[SparseArrays.UMFPACK.JL_UMFPACK_ORDERING] =
        SparseArrays.UMFPACK.LibSuiteSparse.UMFPACK_ORDERING_METIS
    control
end

"""
    function InversionWorkspace(inversion)

Preallocated storage for evaluating the cost function and its gradient.

# Arguments
- `inversion::Inversion`: inverse problem
# Output
- workspace with the controls `u` and their shadow `du`, and the UMFPACK
  factorization `F` of the water-mass matrix (METIS ordering), refactored in
  place at each control vector, with its shadow `dF`
"""
struct InversionWorkspace{I,U,L}
    inversion::I
    u::U
    du::U
    F::L
    dF::L
end
function InversionWorkspace(inversion::Inversion)
    u = deepcopy(inversion.u)
    F = lu(watermassmatrix(inversion.m₀, inversion.y.γ); control=METIS_CONTROL)
    return InversionWorkspace(inversion, u, Enzyme.make_zero(u), F, Enzyme.make_zero(F))
end

"""
    function costgradient!(g, workspace, uvec)

Cost function at `uvec`, with its gradient written to `g`, from one Enzyme
reverse pass through `costfunction`.

# Arguments
- `g`: gradient, overwritten
- `workspace::InversionWorkspace`
- `uvec`: control vector
# Output
- `J`: cost function
"""
function costgradient!(g, workspace::InversionWorkspace, uvec)
    fill!(g, 0)
    Enzyme.remake_zero!(workspace.du)
    return Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        Enzyme.Const(costfunction), Enzyme.Active, Enzyme.Duplicated(uvec, g),
        Enzyme.Const(workspace.inversion), Enzyme.Duplicated(workspace.u, workspace.du),
        Enzyme.Duplicated(workspace.F, workspace.dF))[2]
end

"""
    function gradient_check(inversion; n=5, rtol=1e-4, atol=1e-7, ε=1e-4, workspace)

Compare the Enzyme gradient at `uvec = 0` with centered finite differences at up
to `n` controls in each of the `b`, `q`, and `m` blocks, at the midpoints of `n`
equal segments of the block. Each block is compared on its own, so large
derivatives in one block cannot hide errors in another. The finite difference
is computed from the sampled tracers `ỹ`,
`J(δ) - J(-δ) = (ỹ₊ - ỹ₋)ᵀ Wⁱ (ỹ₊ + ỹ₋ - 2y)`, which avoids subtracting two
large values of `J`; the prior costs are even in `uvec` and cancel.

# Arguments
- `inversion::Inversion`: inverse problem
# Optional Arguments
- `n`: number of controls checked in each block
- `rtol`, `atol`: tolerances of the comparison within each block
- `ε`: finite-difference step
- `workspace`: `InversionWorkspace` to reuse
# Output
- `NamedTuple` keyed by block of `(indices, ∇J, ∇J_finite)`; errors if a block disagrees
"""
function gradient_check(inversion::Inversion; n::Integer=5, rtol=1e-4, atol=1e-7,
    ε=1e-4, workspace=InversionWorkspace(inversion))
    y = inversion.y
    uvec = zero(inversion.σ)
    g = similar(uvec)
    costgradient!(g, workspace, uvec)
    checks = map(inversion.ranges) do r
        k = min(n, length(r))
        indices = isempty(r) ? Int[] : r[ceil.(Int, (length(r) / k) .* ((1:k) .- 0.5))]
        ∇J_finite = map(indices) do i
            δ = zero(uvec)
            δ[i] = ε
            ỹ₊ = observe(forwardmodel(uvec + δ, inversion, workspace.u, workspace.F).c, y)
            ỹ₋ = observe(forwardmodel(uvec - δ, inversion, workspace.u, workspace.F).c, y)
            dot(ỹ₊ - ỹ₋, y.Wⁱ, ỹ₊ + ỹ₋ - 2vec(y)) / 2ε
        end
        (; indices, ∇J=g[indices], ∇J_finite)
    end
    for (block, check) in pairs(checks)
        isempty(check.indices) && continue
        for (i, g, g_finite) in zip(check.indices, check.∇J, check.∇J_finite)
            println("  $block control $i: Enzyme $g, finite difference $g_finite")
        end
        difference = norm(check.∇J - check.∇J_finite)
        scale = max(norm(check.∇J), norm(check.∇J_finite))
        println("  $block relative difference ", iszero(scale) ? difference : difference / scale)
        isapprox(check.∇J, check.∇J_finite; rtol, atol) ||
            error("gradient check failed for $block controls $(check.indices)")
    end
    return checks
end

"""
    function runinversion(inversion, optimizer; name, iterations, checkpoint_interval,
        checkpoint_directory, data_directory, number_of_gradient_checks=5)

Check the gradient, then minimize `J` with Ipopt subject to the limits on
`uvec` and mass conservation. Ipopt works on `z = uvec ./ σ`, starting at the
first guess `z = 0`. The observations go to
`checkpoint_directory/name_observations.nc`, every `checkpoint_interval`-th
iteration to `checkpoint_directory/name_iteration_NNNNNN.nc`, and the result to
`data_directory/name_final.nc`, all written by `save`.

# Arguments
- `inversion::Inversion`: inverse problem
- `optimizer::IpoptOptimizer`: Ipopt settings; the caller's copy is not modified
- `name`: experiment name for output files
- `iterations`: maximum number of Ipopt iterations
- `checkpoint_interval`: iterations between saved checkpoints
# Optional Arguments
- `checkpoint_directory`, `data_directory`: output directories
- `number_of_gradient_checks`: controls checked in each block before the run; 0 skips the check
# Output
- `file`: NetCDF file with the final iteration
"""
function runinversion(inversion::Inversion, optimizer::IpoptOptimizer; name::AbstractString,
    iterations::Integer, checkpoint_interval::Integer,
    checkpoint_directory::AbstractString=TMI.pkgdir("checkpoints", name),
    data_directory::AbstractString=TMI.pkgdatadir(), number_of_gradient_checks::Integer=5)
    σ = inversion.σ
    workspace = InversionWorkspace(inversion)
    uvec = similar(σ)
    # cost alone at Ipopt's trial points, most of which it rejects; gradient only when asked
    J(z, _=nothing) = (uvec .= σ .* z; costfunction(uvec, inversion, workspace.u, workspace.F))
    gJ!(gz, z, _=nothing) = (uvec .= σ .* z; costgradient!(gz, workspace, uvec); gz .*= σ; nothing)
    z₀ = zero(σ)
    (; cons, cons_j, cons_jac_prototype, lcons, ucons) = massconservation(inversion)
    # Optimization.jl requires an AD type; the gradient comes from gJ! and the Jacobian from cons_j
    f = OptimizationFunction(J, AutoEnzyme(); grad=gJ!, cons, cons_j, cons_jac_prototype)
    problem = OptimizationProblem(f, z₀; lb=inversion.lb ./ σ, ub=inversion.ub ./ σ, lcons, ucons)
    number_of_gradient_checks > 0 && gradient_check(inversion; n=number_of_gradient_checks, workspace)

    checkpoint(k) = joinpath(checkpoint_directory, "$(name)_iteration_$(lpad(k, 6, '0')).nc")
    save(joinpath(checkpoint_directory, "$(name)_observations.nc"), inversion)
    save(checkpoint(0), workspace, σ .* z₀; iteration=0)
    callback = (state, _) -> begin
        state.iter > 0 && state.iter % checkpoint_interval == 0 &&
            save(checkpoint(state.iter), workspace, σ .* state.u; iteration=state.iter)
        false
    end
    optimizer = deepcopy(optimizer)
    logfile = joinpath(checkpoint_directory, "ipopt_output.txt")
    optimizer.additional_options["output_file"] = logfile
    println("Starting Ipopt; writing log to $logfile")
    solution = solve(problem, optimizer; maxiters=iterations, verbose=5, callback)
    return save(joinpath(data_directory, "$(name)_final.nc"), workspace, σ .* solution.u;
        iteration=solution.stats.iterations)
end
