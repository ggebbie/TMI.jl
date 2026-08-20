"""
    EnzymeCostGradientCache(inversion)
    evaluate_objective_and_gradient!(cache, x)

Cache and evaluate one Enzyme objective-and-gradient result for Ipopt callbacks.

# Arguments
- `inversion`, `cache`, `x`: inputs, cache, and optimizer controls

# Output
- `result`: invalid cache or `nothing` after updating it
"""
mutable struct EnzymeCostGradientCache{I,W,DW,X,G,T,M}
    inversion::I
    workspace::W
    workspace_gradient::DW
    x::X
    gradient::G
    J::T
    valid::Bool
    mode::M
end
function EnzymeCostGradientCache(inversion::Inversion)
    x₀ = inversion.x₀
    controls = unvec(inversion, x₀)
    workspace = (
        boundary_values=copy(inversion.boundary_template_values),
        boundary_conditions=controls.boundary_conditions,
        source_adjustments=zero(inversion.control_priors.source),
        mass_fractions=controls.mass_fractions,
    )
    return EnzymeCostGradientCache(inversion, workspace, Enzyme.make_zero(workspace),
        similar(x₀), zero(x₀), zero(eltype(x₀)), false,
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal))
end
function evaluate_objective_and_gradient!(cache::EnzymeCostGradientCache, x)
    if !cache.valid || !isequal(x, cache.x)
        fill!(cache.gradient, zero(eltype(cache.gradient)))
        Enzyme.remake_zero!(cache.workspace_gradient)
        cache.J = Enzyme.autodiff(cache.mode, Enzyme.Const(costfunction),
            Enzyme.Active, Enzyme.Duplicated(x, cache.gradient),
            Enzyme.Const(cache.inversion),
            Enzyme.Duplicated(cache.workspace, cache.workspace_gradient))[2]
        copyto!(cache.x, x)
        cache.valid = true
    end
    return nothing
end
"""
    gradient_check(inversion; n=5, seed=1234, control_ranges=nothing,
        indices=nothing, rtol=1e-4, atol=1e-7, ε=1e-6)

Compare Enzyme derivatives with centered finite differences.

# Arguments
- `inversion`, `n`, `seed`: inversion and deterministic selection settings
- `control_ranges`, `indices`, `rtol`, `atol`, `ε`: selection and comparison settings

# Output
- `comparison`: selected indices and both derivative estimates
"""
function gradient_check(inversion::Inversion; n::Integer=5, seed::Integer=1234,
    control_ranges=nothing, indices=nothing, rtol::Real=1e-4,
    atol::Real=1e-7, ε::Real=1e-6,
    cache::EnzymeCostGradientCache=EnzymeCostGradientCache(inversion))
    x = inversion.x₀
    rng = MersenneTwister(seed)
    if !isnothing(indices)
        selected_indices = collect(indices)
    elseif !isnothing(control_ranges)
        selected_indices = [rand(rng, range) for range in control_ranges]
    else
        selected_indices = randperm(rng, length(x))[1:n]
    end
    println("Checking Enzyme gradient at $(length(selected_indices)) controls...")
    flush(stdout)
    evaluate_objective_and_gradient!(cache, x)
    analytic_gradient = cache.gradient
    finite_difference = similar(x, length(selected_indices))
    for (result_index, control_index) in enumerate(selected_indices)
        step = ε * max(one(eltype(x)), abs(x[control_index]))
        x_plus, x_minus = copy(x), copy(x)
        x_plus[control_index] += step
        x_minus[control_index] -= step
        finite_difference[result_index] = (costfunction(x_plus, inversion) -
            costfunction(x_minus, inversion)) / (2step)
        println("  control $control_index: Enzyme $(analytic_gradient[control_index]), " *
            "finite difference $(finite_difference[result_index])")
    end
    isapprox(analytic_gradient[selected_indices], finite_difference; rtol, atol) ||
        error("Enzyme gradient check failed at controls $selected_indices")
    println("Gradient check passed.")
    return (indices=selected_indices,
        analytic=analytic_gradient[selected_indices], finite_difference)
end

"""
    runinversion(inversion, optimizer, checkpointer; max_iterations,
        number_of_gradient_checks=5, random_seed=1234,
        normalize_controls=true)

Check the gradient and run the checkpointed Ipopt inversion.

# Arguments
- `inversion`, `optimizer`, `checkpointer`: scientific, solver, and output setup
- `max_iterations`, `number_of_gradient_checks`, `random_seed`: run limits and checks
- `normalize_controls`: use dimensionless prior departures as Ipopt controls

# Output
- `filename`: terminal inversion NetCDF path
"""
function runinversion(inversion::Inversion, optimizer::IpoptOptimizer,
    checkpointer::InversionCheckpointer; max_iterations::Integer,
    number_of_gradient_checks::Integer=5, random_seed::Integer=1234,
    normalize_controls::Bool=true)
    println("Preparing Enzyme/Ipopt problem...")
    flush(stdout)
    if normalize_controls
        control_offset = inversion.x_prior
        control_scale = Vector(vcat(inv.(sqrt.(diag(inversion.Q_b))),
            inv.(sqrt.(diag(inversion.Q_q))),
            inv.(sqrt.(diag(inversion.Q_m)))))
    else
        control_offset = zero(inversion.x₀)
        control_scale = ones(eltype(inversion.x₀), length(inversion.x₀))
    end
    optimizer_initial_controls =
        (inversion.x₀ - control_offset) ./ control_scale
    optimizer_lower_bounds =
        (inversion.bounds.lower - control_offset) ./ control_scale
    optimizer_upper_bounds =
        (inversion.bounds.upper - control_offset) ./ control_scale
    physical_controls = similar(inversion.x₀)
    cache = EnzymeCostGradientCache(inversion)
    objective(optimizer_controls, _=nothing) = begin
        @. physical_controls = control_offset + control_scale * optimizer_controls
        evaluate_objective_and_gradient!(cache, physical_controls)
        cache.J
    end
    gradient!(destination, optimizer_controls, _=nothing) = begin
        @. physical_controls = control_offset + control_scale * optimizer_controls
        evaluate_objective_and_gradient!(cache, physical_controls)
        @. destination = control_scale * cache.gradient
        nothing
    end
    conservation = mass_conservation_constraints(
        inversion, control_offset, control_scale)
    optimization_function = OptimizationFunction(objective, AutoEnzyme();
        grad=gradient!, conservation.callbacks...)
    problem = OptimizationProblem(optimization_function,
        optimizer_initial_controls, inversion;
        lb=optimizer_lower_bounds, ub=optimizer_upper_bounds,
        conservation.bounds...)
    gradient_check(inversion; cache,
        n=number_of_gradient_checks, seed=random_seed)
    mkpath(checkpointer.checkpoint_directory)
    mkpath(checkpointer.data_directory)
    saveobservations(checkpointer, inversion)
    initial_objective = objective(optimizer_initial_controls, inversion)
    saveinversion(checkpointer, inversion.x₀, inversion;
        iteration=0, objective=initial_objective)
    checkpoint_callback = (state, value) -> begin
        if state.iter > 0 && state.iter % checkpointer.interval == 0
            checkpoint_controls =
                control_offset + control_scale .* state.u
            saveinversion(checkpointer, checkpoint_controls, inversion;
                iteration=state.iter, objective=value)
        end
        false
    end
    output_file = joinpath(checkpointer.checkpoint_directory, "ipopt_output.txt")
    optimizer.additional_options["output_file"] = output_file
    optimizer.additional_options["file_print_level"] = 5
    println("Starting Ipopt; writing log to $output_file")
    solution = solve(problem, optimizer; maxiters=max_iterations,
        verbose=5, callback=checkpoint_callback)
    final_controls = control_offset + control_scale .* solution.u
    final_objective = objective(solution.u, inversion)
    return saveinversion(checkpointer, final_controls, inversion;
        iteration=solution.stats.iterations, objective=final_objective,
        terminal=true)
end
