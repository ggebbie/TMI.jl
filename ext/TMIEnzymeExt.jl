# TMIEnzymeExt is loaded automatically
# whenever both `TMI` and `Enzyme` are present in the active environment. 

module TMIEnzymeExt

using Enzyme
using SparseArrays
using TMI

using Enzyme: Annotation
using TMI: BoundaryCondition, Field, Grid, MassFraction, Source
using TMI: cartesianindex, interpweights, observe, step_cartesian, unvec, watermassmatrix, wet

import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth
import Enzyme.EnzymeRules: augmented_primal, needs_primal, needs_shadow, reverse

include("enzyme/gunvec_boundarycondition.jl")
include("enzyme/gunvec_field.jl")
include("enzyme/gwatermassmatrix.jl")
include("enzyme/gldiv_field.jl")
include("enzyme/accumulate_into.jl")
include("enzyme/gobserve.jl")

end
