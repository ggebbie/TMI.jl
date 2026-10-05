module TMIEnzymeOptimizationExt

using Enzyme, LinearAlgebra, NCDatasets, Optimization, SparseArrays, TMI
using OptimizationIpopt: IpoptOptimizer

# TMI controls exceed Enzyme's default type-analysis lattice.
function __init__()
    Enzyme.API.maxtypeoffset!(8192)
    Enzyme.API.maxtypedepth!(12)
end

using Enzyme: Annotation
using TMI: MassFraction, Source, step_cartesian
import TMI: Inversion, Observations, gradient_check, observe, runinversion, steadyinversion
import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth, augmented_primal, needs_primal, needs_shadow, reverse

include("enzyme_optimization/weighting_matrices.jl")
include("enzyme_optimization/observations.jl")
include("enzyme_optimization/inversion.jl")
include("enzyme_optimization/constraints.jl")
include("enzyme_optimization/output.jl")
include("enzyme_optimization/optimization.jl")

include("enzyme_optimization/enzyme_rules/gunvec_inplace.jl")
include("enzyme_optimization/enzyme_rules/gwatermassmatrix.jl")
include("enzyme_optimization/enzyme_rules/glu.jl")
include("enzyme_optimization/enzyme_rules/gldiv_field.jl")
# Enzyme differentiates `observe` natively, so it has no rule.
end
