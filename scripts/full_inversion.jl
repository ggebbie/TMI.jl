using Enzyme
using Optimization
using OptimizationIpopt
using TMI

const TMIEnzymeOptimization = Base.get_extension(TMI, :TMIEnzymeOptimizationExt)
using .TMIEnzymeOptimization

"""
    woceobservations()

Read the gridded WOCE tracer observations, standard deviations, and horizontal
decorrelation scale used by the full inversion example.

# Arguments
This function reads the paths configured in its body and takes no arguments.

# Output
- `observations::Observations`: WOCE `y`, `σ`, and `L` on the modern TMI grid
"""
function woceobservations()
    data_directory = normpath(TMI.pkgdir("..", "regridding_WOCE_for_TMI", "data"))
    observations_file = joinpath(data_directory,
        "TMI_gridded_WOCE_4x4_Variables_33_levels.nc")
    uncertainties_file = joinpath(data_directory,
        "TMI_gridded_WOCE_Errors_4x4_Variables_33_levels.nc")
    woce_tracers = (
        θ="potential_temperature", S="salinity", δ¹⁸O="d18o",
        PO₄="phosphate", NO₃="nitrate", O₂="oxygen",
    )
    γ = TMI.Grid(TMI.pkgdatadir("TMI_modern_90x45x33_G14_v2.nc"))
    observations = map(name -> readfield(observations_file, name, γ), woce_tracers)
    σ = map((name, observation) -> readfield(
        uncertainties_file,
        name,
        γ;
        name=Symbol("σ", observation.name),
    ), woce_tracers, observations)
    L = readfield(uncertainties_file, "decorrelation_length", γ)
    return Observations(observations, σ, γ; L)
end

"""
    main(experiment_name; max_iterations, checkpoint_interval,
        number_of_gradient_checks, random_seed)

Run the full WOCE inversion with Enzyme gradients and Ipopt.

# Arguments
- `experiment_name` and keyword arguments: artifact name and run limits

# Output
- `checkpointer::InversionCheckpointer`: configuration for the saved results
"""
function main(
    experiment_name = "full_inversion";
    max_iterations = 6000,
    checkpoint_interval = 500,
    number_of_gradient_checks = 5,
    random_seed = 1234,
)
    checkpointer = InversionCheckpointer(experiment_name, checkpoint_interval)
    observations = woceobservations()
    γ = observations.γ
    y = observations.y
    σ = observations.σ
    b₀ = map(getsurfaceboundary, y)
    σ_b = map(getsurfaceboundary, σ)
    # Keep each tracer's prior, uncertainty, smoothing scale, and bounds together.
    boundary_controls = (
        θ = (prior=b₀.θ, σ=σ_b.θ,
            L=1_000.0, lb=-2.0, ub=35.0),
        S = (prior=b₀.S, σ=σ_b.S,
            L=1_000.0, lb=0.0, ub=45.0),
        δ¹⁸O = (prior=b₀.δ¹⁸O, σ=σ_b.δ¹⁸O,
            L=1_000.0, lb=-10.0, ub=10.0),
        PO₄ = (prior=b₀.PO₄, σ=σ_b.PO₄,
            L=1_000.0, lb=0.0, ub=10.0),
        NO₃ = (prior=b₀.NO₃, σ=σ_b.NO₃,
            L=1_000.0, lb=0.0, ub=45.0),
        O₂ = (prior=b₀.O₂, σ=σ_b.O₂,
            L=1_000.0, lb=0.0, ub=500.0),
    )
    inversion = Inversion(
        observations;
        boundary_conditions=boundary_controls,
        sources=(qPO₄=(prior=log(3.0e-4) * one(onesource(γ, :qPO₄,
                "local source of phosphate", "μmol/kg"; logscale=true)),
            σ=log(50.0),
            lb=log(6.0e-6), ub=log(1.5e-2)),),
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
        limited_memory_max_history = 10,
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
