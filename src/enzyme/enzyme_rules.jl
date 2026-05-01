"""
    TMIEnzymeRules

This module collects the custom reverse-mode Enzyme rules for TMI.jl.

The documentation for these rules, including the mathematical derivations,
is in the `src/enzyme/tutorial/` directory.
"""
module TMIEnzymeRules

using Enzyme
using SparseArrays

using Enzyme: Annotation

using ..TMI: BoundaryCondition, Field, Grid, MassFraction
using ..TMI: cartesianindex, step_cartesian, unvec, watermassmatrix

import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth
import Enzyme.EnzymeRules: augmented_primal, needs_primal, needs_shadow, reverse

include("unvec_boundarycondition.jl")
include("unvec_field.jl")
include("watermassmatrix.jl")
include("ldiv_field.jl")

end
