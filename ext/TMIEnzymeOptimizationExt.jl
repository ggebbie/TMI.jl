module TMIEnzymeOptimizationExt

using Enzyme, LinearAlgebra, NCDatasets, Optimization, Random, SparseArrays, TMI
using OptimizationIpopt: IpoptOptimizer

# TMI controls exceed Enzyme's default type-analysis lattice.
function __init__()
    Enzyme.API.maxtypeoffset!(8192)
    Enzyme.API.maxtypedepth!(12)
end

using Enzyme: Annotation
using TMI: BoundaryCondition, Field, Grid, MassFraction, Source, cartesianindex,
    interpweights, observe, step_cartesian, boundary_smoothness_precision_matrix,
    surfacelaplacianmatrix, watermassmatrix, wet
import TMI: gradient_check, unvec, unvec!
import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth, augmented_primal, needs_primal, needs_shadow, reverse

export Inversion, InversionCheckpointer, runinversion, costfunction,
    data_cost, boundary_prior_cost, boundary_smoothness_cost,
    source_prior_cost, mass_fraction_prior_cost

include("enzyme_optimization/parameters.jl")
include("enzyme_optimization/cost_functions.jl")
include("enzyme_optimization/constraints.jl")
include("enzyme_optimization/checkpoints.jl")
include("enzyme_optimization/optimization.jl")

include("enzyme_optimization/enzyme_rules/gunvec_inplace.jl")
include("enzyme_optimization/enzyme_rules/gwatermassmatrix.jl")
include("enzyme_optimization/enzyme_rules/glu.jl")
include("enzyme_optimization/enzyme_rules/gldiv_field.jl")
# Enzyme 0.13.190 handles TMI accumulator merging without a custom overload.
include("enzyme_optimization/enzyme_rules/gobserve.jl")
end
