using Enzyme, NCDatasets, Optimization, OptimizationIpopt, SparseArrays, Test, TMI
using LinearAlgebra: Diagonal, diag, dot, inv, issymmetric, lu
const TMIEnzymeOptimization = Base.get_extension(TMI, :TMIEnzymeOptimizationExt)
using .TMIEnzymeOptimization
function inversioninputs(n=5)
    coordinates = collect(range(-1.0, 1.0; length=n))
    interior = falses(n, n, n)
    interior[2:end-1, 2:end-1, 2:end-1] .= true
    neighbors = CartesianIndex.(((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))
    γ = Grid((180coordinates, 60coordinates, 500(coordinates .+ 1)), trues(n, n, n), interior, (false, false, false), collect(neighbors))
    fields = (active=Field(ones(n, n, n), γ, :active, "active", "unitless"), fixed=Field(fill(2.0, n, n, n), γ, :fixed, "fixed", "unitless"))
    boundary = map(getsurfaceboundary, fields)
    sources = (active_source=onesource(γ, :active_source, "active", "unitless"), fixed_source=2onesource(γ, :fixed_source, "fixed", "unitless"))
    stoichiometry = (active=(active_source=1.0,), fixed=(fixed_source=1.0,))
    fractions = massfractions_isotropic(γ)
    foreach(fraction -> replace!(fraction.fraction, NaN => 0.0), fractions)
    mass_fraction_uncertainty = unvec(fractions,
        collect(range(0.2, 0.8; length=length(vec(fractions)))))
    A = watermassmatrix(fractions, γ)
    recovered_fractions = massfractions(A, γ)
    @test keys(recovered_fractions) == (:north, :east, :south, :west, :up, :down)
    @test all(fraction -> fraction isa TMI.MassFraction, recovered_fractions)
    @test vec(recovered_fractions) ≈ vec(fractions)
    @test watermassmatrix(recovered_fractions, γ) ≈ A
    row_sign = ones(size(A, 1))
    row_sign[findall(γ.interior[γ.wet])] .= -1.0
    legacy_A = Diagonal(row_sign) * A
    @test vec(massfractions(legacy_A, γ)) ≈ vec(fractions)
    @test vec(one(boundary)) == vec(one(one(boundary))); @test all(iszero, vec(zero(sources)))
    @test vec(2one(fractions)) == fill(2.0, length(vec(fractions))); @test vec(unvec(fractions, vec(fractions))) == vec(fractions)
    reference = steadyinversion(lu(A), boundary, (active=sources.active_source, fixed=sources.fixed_source), γ)
    observations = Observations(map(field -> 0.95field, reference), 0.5, γ)
    inversion = Inversion(observations; boundary_conditions=(active=(prior=boundary.active, σ=0.5), fixed=boundary.fixed),
        sources=(active_source=(prior=zerosource(γ, :active_source, "active", "unitless"; logscale=true), σ=1.0, lb=log(0.25), ub=log(4.0)), fixed_source=sources.fixed_source), stoichiometry,
        mass_fractions=(prior=fractions, σ=mass_fraction_uncertainty,
            lb=0.0, ub=1.0))
    return (; γ, fields, boundary, sources, stoichiometry, fractions,
        mass_fraction_uncertainty, reference, observations, inversion)
end
@testset "Enzyme inversion" begin
    inputs = inversioninputs(); (; γ, inversion) = inputs
    predicted = steadyinversion(lu(watermassmatrix(inputs.fractions, γ)), inputs.boundary, inputs.sources, inputs.observations; stoichiometry=inputs.stoichiometry)
    @test keys(predicted) == keys(inputs.observations.y)
    @test vec(inputs.observations) == vec(map(field -> field.tracer[γ.interior], inputs.observations.y))
    locations = [(γ.lon[i], γ.lat[i], γ.depth[i]) for i in (2, 2, 4)]
    sparse_values = observe(inputs.reference.fixed, locations, γ)
    sparse_values[2] += 2
    mixed = Observations((active=inputs.observations.y.active, fixed=sparse_values), (active=0.5, fixed=fill(0.5, 3)), γ; locations=(fixed=locations,))
    mixed_prediction = steadyinversion(lu(watermassmatrix(inputs.fractions, γ)), inputs.boundary, inputs.sources, mixed; stoichiometry=inputs.stoichiometry)
    @test length(mixed_prediction.active) == sum(γ.interior); @test length(mixed_prediction.fixed) == 3
    @test diag(inv(Matrix(inputs.observations.Q))) ≈ fill(0.25, 2sum(γ.interior))
    σ = Field(fill(2.0, size(γ.wet)), γ, :σ, "uncertainty", "unitless")
    L = Field(fill(1_000.0, size(γ.wet)), γ, :L, "scale", "km")
    correlated = Observations((active=inputs.reference.active,), (active=σ,), γ; L)
    @test correlated.Q isa SparseMatrixCSC && issymmetric(correlated.Q)
    supplied = Observations((active=sparse_values,), fill(0.5, 3), γ; locations=(active=locations,), Q=(active=Diagonal(fill(4.0, 3)),))
    @test supplied.Q ≈ Diagonal(fill(4.0, 3))
    nan_precision = Observations((active=sparse_values,), 0.5, γ; locations=(active=locations,), Q=(active=Diagonal(fill(NaN, 3)),))
    @test isnan(data_cost(ones(3), nan_precision.Q))
    source_range = inversion.control_ranges.source; controls = unvec(inversion, inversion.x₀)
    @test all(iszero, inversion.x₀[source_range]); @test inversion.bounds.lower[source_range] ≈ fill(log(0.25), length(source_range))
    @test inversion.bounds.upper[source_range] ≈ fill(log(4.0), length(source_range)); @test vec(controls.sources.active_source) == vec(inputs.sources.active_source)
    @test !controls.sources.active_source.logscale
    tripled = copy(inversion.x₀); tripled[source_range] .= log(3.0)
    q = unvec(inversion, tripled).sources.active_source
    @test vec(q) ≈ fill(3.0, length(q)); @test vec(-2.0q) ≈ fill(-6.0, length(q))
    control_scale = Vector(vcat(inv.(sqrt.(diag(inversion.Q_b))),
        inv.(sqrt.(diag(inversion.Q_q))), inv.(sqrt.(diag(inversion.Q_m)))))
    normalized_initial_controls = (inversion.x₀ - inversion.x_prior) ./ control_scale
    normalized_lower_bounds =
        (inversion.bounds.lower - inversion.x_prior) ./ control_scale
    normalized_upper_bounds =
        (inversion.bounds.upper - inversion.x_prior) ./ control_scale
    @test all(iszero, normalized_initial_controls)
    @test inversion.x_prior + control_scale .* normalized_initial_controls ≈ inversion.x₀
    @test all(isapprox.(inversion.x_prior +
        control_scale .* normalized_lower_bounds, inversion.bounds.lower;
        atol=eps()))
    @test all(isapprox.(inversion.x_prior +
        control_scale .* normalized_upper_bounds, inversion.bounds.upper;
        atol=eps()))
    @test control_scale[inversion.control_ranges.mass_fraction] ≈
        vec(inputs.mass_fraction_uncertainty)
    x = copy(inversion.x₀); x[first(inversion.control_ranges.boundary)] += 0.2
    x[first(source_range)] += 0.1; x[first(inversion.control_ranges.mass_fraction)] += 0.05
    deviation = x - inversion.x_prior; ranges = inversion.control_ranges
    @test boundary_prior_cost(deviation[ranges.boundary], inversion.Q_b) ≈ 0.2^2 / 0.5^2; @test iszero(boundary_smoothness_cost(deviation[ranges.boundary], inversion.Q_s))
    @test source_prior_cost(deviation[ranges.source], inversion.Q_q) ≈ 0.1^2
    @test mass_fraction_prior_cost(deviation[ranges.mass_fraction], inversion.Q_m) ≈
        (0.05 / first(vec(inputs.mass_fraction_uncertainty)))^2
    normalized_controls = fill(0.01, length(inversion.x₀))
    physical_controls = inversion.x_prior + control_scale .* normalized_controls
    @test (physical_controls - inversion.x_prior) ./ control_scale ≈ normalized_controls
    normalized_cache = TMIEnzymeOptimization.EnzymeCostGradientCache(inversion)
    TMIEnzymeOptimization.evaluate_objective_and_gradient!(
        normalized_cache, physical_controls)
    normalized_gradient = control_scale .* normalized_cache.gradient
    for control_index in first.(filter(!isempty, values(inversion.control_ranges)))
        step = 1.0e-6
        controls_plus, controls_minus =
            copy(normalized_controls), copy(normalized_controls)
        controls_plus[control_index] += step
        controls_minus[control_index] -= step
        finite_difference = (costfunction(inversion.x_prior +
            control_scale .* controls_plus, inversion) -
            costfunction(inversion.x_prior + control_scale .* controls_minus,
                inversion)) / (2step)
        @test normalized_gradient[control_index] ≈ finite_difference rtol=1.0e-4
    end
    # Conservation and every active gradient block remain valid.
    conservation = TMIEnzymeOptimization.mass_conservation_constraints(inversion); g = zeros(sum(γ.interior))
    conservation.callbacks.cons(g, inversion.x₀, inversion)
    @test g ≈ ones(length(g))
    physical_jacobian = copy(conservation.callbacks.cons_jac_prototype)
    conservation.callbacks.cons_j(physical_jacobian, inversion.x₀, inversion)
    normalized_conservation = TMIEnzymeOptimization.mass_conservation_constraints(
        inversion, inversion.x_prior, control_scale)
    normalized_jacobian = normalized_conservation.callbacks.cons_jac_prototype
    @test normalized_jacobian ≈ physical_jacobian * Diagonal(control_scale)
    @test length(unique(control_scale[inversion.control_ranges.mass_fraction])) > 1
    mass_control_index = first(inversion.control_ranges.mass_fraction)
    step = 1.0e-6
    controls_plus, controls_minus =
        copy(normalized_initial_controls), copy(normalized_initial_controls)
    controls_plus[mass_control_index] += step
    controls_minus[mass_control_index] -= step
    constraint_plus, constraint_minus = similar(g), similar(g)
    normalized_conservation.callbacks.cons(
        constraint_plus, controls_plus, inversion)
    normalized_conservation.callbacks.cons(
        constraint_minus, controls_minus, inversion)
    @test (constraint_plus - constraint_minus) / (2step) ≈
        Vector(normalized_jacobian[:, mass_control_index])
    check = TMI.gradient_check(inversion; control_ranges=filter(!isempty, values(inversion.control_ranges)), rtol=1e-4)
    @test length(check.indices) == 3
    mixed_boundary_controls = (active=(prior=inputs.boundary.active, σ=0.5), fixed=inputs.boundary.fixed)
    mixed_source_controls = (active_source=(prior=inputs.sources.active_source, σ=1.0, lb=-4.0, ub=4.0), fixed_source=inputs.sources.fixed_source)
    mixed_inversion = Inversion(mixed; boundary_conditions=mixed_boundary_controls, sources=mixed_source_controls, stoichiometry=inputs.stoichiometry,
        mass_fractions=(prior=inputs.fractions,
            σ=inputs.mass_fraction_uncertainty, lb=0.0, ub=1.0))
    raw_range = mixed_inversion.control_ranges.source
    @test mixed_inversion.bounds.lower[raw_range] ≈ fill(-5.0, length(raw_range)); @test mixed_inversion.bounds.upper[raw_range] ≈ fill(3.0, length(raw_range))
    raw_controls = copy(mixed_inversion.x₀); raw_controls[first(raw_range)] = -2.0
    sources = unvec(mixed_inversion, raw_controls).sources
    @test first(vec(sources.active_source)) == -1.0
    mixed_cache = TMIEnzymeOptimization.EnzymeCostGradientCache(mixed_inversion)
    TMIEnzymeOptimization.evaluate_objective_and_gradient!(mixed_cache, mixed_inversion.x₀); @test all(isfinite, mixed_cache.gradient)
    mixed_check = TMI.gradient_check(mixed_inversion; control_ranges=filter(!isempty, values(mixed_inversion.control_ranges)), rtol=1e-4); @test length(mixed_check.indices) == 3
    fixed_mass_inversion = Inversion(mixed; boundary_conditions=mixed_boundary_controls, sources=mixed_source_controls,
        stoichiometry=inputs.stoichiometry, mass_fractions=inputs.fractions)
    @test isempty(fixed_mass_inversion.control_ranges.mass_fraction); @test vec(unvec(fixed_mass_inversion, fixed_mass_inversion.x₀).mass_fractions) == vec(inputs.fractions)
    fixed_mass_check = TMI.gradient_check(fixed_mass_inversion;
        control_ranges=filter(!isempty, values(fixed_mass_inversion.control_ranges)), rtol=1e-4)
    @test length(fixed_mass_check.indices) == 2
    mass_only = Inversion(inputs.observations; boundary_conditions=inputs.boundary, sources=inputs.sources, stoichiometry=inputs.stoichiometry,
        mass_fractions=(prior=inputs.fractions, σ=1.0, lb=0.0, ub=1.0))
    mass_cache = TMIEnzymeOptimization.EnzymeCostGradientCache(mass_only)
    TMIEnzymeOptimization.evaluate_objective_and_gradient!(mass_cache, mass_only.x₀); @test all(isfinite, mass_cache.gradient)
    mass_check = TMI.gradient_check(mass_only; control_ranges=(mass_only.control_ranges.mass_fraction,), rtol=1e-4); @test length(mass_check.indices) == 1
    mktempdir(directory -> begin
        checkpointer = InversionCheckpointer("test", 2; checkpoint_directory=joinpath(directory, "checkpoints"),
            data_directory=joinpath(directory, "data"))
        optimizer = IpoptOptimizer(hessian_approximation="limited-memory")
        runinversion(mixed_inversion, optimizer, checkpointer; max_iterations=2, number_of_gradient_checks=0)
        observation_file = joinpath(checkpointer.checkpoint_directory, "test_observations.nc")
        NCDataset(dataset -> begin
            @test dataset["fixed_sample"][:] == sparse_values
            @test dataset["fixed"][:, :, :][nearestneighbor(first(locations), γ)] ≈ sum(sparse_values[1:2]) / 2
            @test dataset["fixed"].attrib["support"] == "sparse"
            @test !any(contains.(keys(dataset), "precision"))
        end, observation_file)
        file = joinpath(checkpointer.checkpoint_directory, "test_iteration_000000.nc")
        NCDataset(dataset -> begin
            @test dataset.attrib["iteration"] == 0
            @test haskey(dataset, "active_source")
            @test !any(startswith.(keys(dataset.attrib), "normalized_misfit"))
        end, file)
        @test isfile(joinpath(checkpointer.checkpoint_directory, "test_iteration_000002.nc")); @test isfile(joinpath(checkpointer.data_directory, "test_final.nc"))
        @test isfile(joinpath(checkpointer.checkpoint_directory, "ipopt_output.txt"))
        physical_checkpointer = InversionCheckpointer("physical", 1;
            checkpoint_directory=joinpath(directory, "physical_checkpoints"),
            data_directory=joinpath(directory, "physical_data"))
        physical_optimizer = IpoptOptimizer(hessian_approximation="limited-memory")
        runinversion(mixed_inversion, physical_optimizer, physical_checkpointer;
            max_iterations=1, number_of_gradient_checks=0,
            normalize_controls=false)
        physical_initial_file = joinpath(physical_checkpointer.checkpoint_directory,
            "physical_iteration_000000.nc")
        NCDataset(dataset -> @test(dataset.attrib["objective"] ≈
            costfunction(mixed_inversion.x₀, mixed_inversion)), physical_initial_file)
        @test isfile(joinpath(physical_checkpointer.data_directory,
            "physical_final.nc"))
    end)
end
