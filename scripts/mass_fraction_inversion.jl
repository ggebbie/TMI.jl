using Enzyme
using Optimization
using OptimizationIpopt
using TMI

const TMIEnzymeOptimization = Base.get_extension(TMI, :TMIEnzymeOptimizationExt)
using .TMIEnzymeOptimization

"""
    main(experiment_name; max_iterations, checkpoint_interval,
        number_of_gradient_checks, random_seed)

Invert directional water-mass fractions with fixed tracer controls.

# Arguments
- `experiment_name` and keyword arguments: artifact name and run limits

# Output
- `checkpointer::InversionCheckpointer`: configuration for the saved results
"""
function main(
    experiment_name = "mass_fraction_inversion";
    max_iterations = 10,
    checkpoint_interval = 1,
    number_of_gradient_checks = 5,
    random_seed = 1234,
)
    input_tmi_version = "modern_90x45x33_G14_v2"
    input_tracers = (
        θ = "θ", S = "Sp", δ¹⁸O = "δ¹⁸Ow", PO₄ = "PO₄", NO₃ = "NO₃", O₂ = "O₂",
    )
    _, _, γ, input_tmi_file, _, _ = config(input_tmi_version)
    y = map(variable -> readfield(input_tmi_file, variable, γ), input_tracers)
    b = map(getsurfaceboundary, y)
    qPO₄ = readsource(input_tmi_file, "qPO₄", γ)
    checkpointer = InversionCheckpointer(experiment_name, checkpoint_interval)
    σ =
        (θ = 0.1, S = 0.01, δ¹⁸O = 0.2, PO₄ = 0.05, NO₃ = 1.0, O₂ = 5.0)
    observations = Observations(y, σ, γ)

    inversion = Inversion(
        observations;
        boundary_conditions=b,
        sources=(; qPO₄),
        stoichiometry=(
            PO₄ = (qPO₄ = -1.0,),
            NO₃ = (qPO₄ = -15.5,),
            O₂ = (qPO₄ = 170.0,),
        ),
        mass_fractions=(prior=massfractions_isotropic(γ),
            σ=1.0, lb=0.0, ub=1.0),
    )
    optimizer = IpoptOptimizer(
        hessian_approximation = "limited-memory",
        limited_memory_max_history = 27,
        acceptable_tol = 1.0e-6,
        mu_strategy = "adaptive",
        adaptive_mu_globalization = "kkt-error",
        nlp_scaling_method = "gradient-based",
    )
    runinversion(
        inversion,
        optimizer,
        checkpointer;
        max_iterations=max_iterations,
        number_of_gradient_checks=number_of_gradient_checks,
        random_seed=random_seed,
        normalize_controls=true,
    )
    return checkpointer
end

main()
