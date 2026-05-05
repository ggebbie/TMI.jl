"""
    TMIEnzymeRules

This module collects the custom reverse-mode Enzyme rules for TMI.jl.

Notation in this directory follows the TMI adjoint convention: `ginput` means
`dJ/dinput`, where `J` is a scalar cost function. The included files define
`augmented_primal` and `reverse` methods for operations that Enzyme should not
differentiate by tracing their generic Julia implementation.
"""
module TMIEnzymeRules

using Enzyme
using SparseArrays

using Enzyme: Annotation

using ..TMI: BoundaryCondition, Field, Grid, MassFraction
using ..TMI: cartesianindex, interpweights, observe, step_cartesian, unvec, watermassmatrix

import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth
import Enzyme.EnzymeRules: augmented_primal, needs_primal, needs_shadow, reverse

include("gunvec_boundarycondition.jl")
include("gunvec_field.jl")
include("gwatermassmatrix.jl")
include("gldiv_field.jl")
include("gobserve.jl")

end
