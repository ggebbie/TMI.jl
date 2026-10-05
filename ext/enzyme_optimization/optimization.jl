# METIS ordering: about 40% less fill than UMFPACK's default AMD ordering on TMI
# grids; `lu!` keeps the ordering of the factorization it refactors
const METIS_CONTROL = let control = SparseArrays.UMFPACK.get_umfpack_control(Float64, Int64)
    control[SparseArrays.UMFPACK.JL_UMFPACK_ORDERING] =
        SparseArrays.UMFPACK.LibSuiteSparse.UMFPACK_ORDERING_METIS
    control
end

"""
    EnzymeCostGradientCache(inversion)

Storage for the cost function and its gradient, so that Ipopt's separate cost
and gradient calls at the same control vector share one Enzyme reverse pass.

# Arguments
- `inversion::Inversion`: inverse problem

# Output
- cache with preallocated controls `u` and their shadow `du`, the UMFPACK
  factorization `F` of the water-mass matrix (METIS ordering), refactored in
  place at each control vector, and its shadow `dF`, the last control vector
  `uvec`, gradient `guvec`, and cost function `J`
"""
mutable struct EnzymeCostGradientCache{I,U,L,V,T}
    inversion::I
    u::U
    du::U
    F::L
    dF::L
    uvec::V
    guvec::V
    J::T
    valid::Bool
end
function EnzymeCostGradientCache(inversion::Inversion)
    u = deepcopy(inversion.u)
    F = lu(watermassmatrix(inversion.m₀, inversion.y.γ); control=METIS_CONTROL)
    return EnzymeCostGradientCache(inversion, u, Enzyme.make_zero(u), F, Enzyme.make_zero(F),
        similar(inversion.σ), zero(inversion.σ), 0.0, false)
end

"""
    costgradient!(cache, uvec)

Update the cost function `cache.J` and gradient `cache.guvec` at `uvec`, unless
the cache already holds them.

# Arguments
- `cache::EnzymeCostGradientCache`
- `uvec`: control vector

# Output
- `nothing`
"""
function costgradient!(cache::EnzymeCostGradientCache, uvec)
    cache.valid && isequal(uvec, cache.uvec) && return nothing
    fill!(cache.guvec, 0)
    Enzyme.remake_zero!(cache.du)
    cache.J = Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        Enzyme.Const(costfunction), Enzyme.Active, Enzyme.Duplicated(uvec, cache.guvec),
        Enzyme.Const(cache.inversion), Enzyme.Duplicated(cache.u, cache.du),
        Enzyme.Duplicated(cache.F, cache.dF))[2]
    copyto!(cache.uvec, uvec)
    cache.valid = true
    return nothing
end

"""
    cost!(cache, uvec)

Cost function at `uvec`: `cache.J` if the cache holds `uvec`, otherwise the cost
alone, without a reverse pass. Most of Ipopt's cost calls are at trial points
it rejects, which need no gradient.

# Arguments
- `cache::EnzymeCostGradientCache`
- `uvec`: control vector

# Output
- `J`: cost function
"""
function cost!(cache::EnzymeCostGradientCache, uvec)
    cache.valid && isequal(uvec, cache.uvec) && return cache.J
    return costfunction(uvec, cache.inversion, cache.u, cache.F)
end

"""
    gradient_check(inversion; n=5, rtol=1e-4, atol=1e-7, ε=1e-4, cache)

Compare the Enzyme gradient at `uvec = 0` with centered finite differences at up
to `n` controls in each of the `b`, `q`, and `m` blocks, at the midpoints of `n`
equal segments of the block. Each block is compared on its own, so large
derivatives in one block cannot hide errors in another. The finite difference
is computed from the sampled tracers `ỹ`,
`J(δ) - J(-δ) = (ỹ₊ - ỹ₋)ᵀ Wⁱ (ỹ₊ + ỹ₋ - 2y)`, which avoids subtracting two
large values of `J`; the control penalty is even in `uvec` and cancels.

# Arguments
- `inversion::Inversion`: inverse problem
- `n`: number of controls checked in each block
- `rtol`, `atol`: tolerances of the comparison within each block
- `ε`: finite-difference step
- `cache`: `EnzymeCostGradientCache` to reuse

# Output
- `NamedTuple` keyed by block of `(indices, ∇J, ∇J_finite)`; errors if a block disagrees
"""
function gradient_check(inversion::Inversion; n::Integer=5, rtol=1e-4, atol=1e-7,
    ε=1e-4, cache=EnzymeCostGradientCache(inversion))
    y = inversion.y
    uvec = zero(inversion.σ)
    costgradient!(cache, uvec)
    checks = map(inversion.ranges) do r
        k = min(n, length(r))
        indices = isempty(r) ? Int[] : r[ceil.(Int, (length(r) / k) .* ((1:k) .- 0.5))]
        ∇J_finite = map(indices) do i
            δ = zero(uvec)
            δ[i] = ε
            ỹ₊ = observe(steadyinversion(uvec + δ, inversion, cache.u, cache.F), y)
            ỹ₋ = observe(steadyinversion(uvec - δ, inversion, cache.u, cache.F), y)
            dot(ỹ₊ - ỹ₋, y.Wⁱ, ỹ₊ + ỹ₋ - 2vec(y)) / 2ε
        end
        (; indices, ∇J=cache.guvec[indices], ∇J_finite)
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
    runinversion(inversion, optimizer; name, iterations, checkpoint_interval,
        checkpoint_directory, data_directory, number_of_gradient_checks=5)

Check the gradient, then minimize `J` with Ipopt subject to the limits on
`uvec` and mass conservation. Ipopt works on `z = uvec ./ σ`, starting at the
first guess `z = 0`. Iterations are saved to `checkpoint_directory` every
`checkpoint_interval` iterations and the result to `data_directory`.

# Arguments
- `inversion::Inversion`: inverse problem
- `optimizer::IpoptOptimizer`: Ipopt settings; the caller's copy is not modified
- `name`: experiment name for output files
- `iterations`: maximum number of Ipopt iterations
- `checkpoint_interval`: iterations between saved checkpoints
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
    cache = EnzymeCostGradientCache(inversion)
    uvec = similar(σ)
    J(z, _=nothing) = (uvec .= σ .* z; cost!(cache, uvec))
    gJ!(gz, z, _=nothing) = (uvec .= σ .* z; costgradient!(cache, uvec); gz .= σ .* cache.guvec; nothing)
    z₀ = zero(σ)
    (; cons, cons_j, cons_jac_prototype, lcons, ucons) = massconservation(inversion)
    f = OptimizationFunction(J, AutoEnzyme(); grad=gJ!, cons, cons_j, cons_jac_prototype)
    problem = OptimizationProblem(f, z₀; lb=inversion.lb ./ σ, ub=inversion.ub ./ σ, lcons, ucons)
    number_of_gradient_checks > 0 && gradient_check(inversion; n=number_of_gradient_checks, cache)

    saveobservations(name, checkpoint_directory, inversion)
    saveinversion(name, checkpoint_directory, zero(σ), cache; iteration=0, J=J(z₀))
    callback = (state, Jₖ) -> begin
        state.iter > 0 && state.iter % checkpoint_interval == 0 &&
            saveinversion(name, checkpoint_directory, σ .* state.u, cache;
                iteration=state.iter, J=Jₖ)
        false
    end
    optimizer = deepcopy(optimizer)
    logfile = joinpath(checkpoint_directory, "ipopt_output.txt")
    optimizer.additional_options["output_file"] = logfile
    optimizer.additional_options["file_print_level"] = 5
    println("Starting Ipopt; writing log to $logfile")
    solution = solve(problem, optimizer; maxiters=iterations, verbose=5, callback)
    return saveinversion(name, data_directory, σ .* solution.u, cache;
        iteration=solution.stats.iterations, J=solution.objective, final=true)
end
