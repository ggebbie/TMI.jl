# TMIEnzymeExt is loaded automatically
# whenever both `TMI` and `Enzyme` are present in the active environment. 

module TMIEnzymeExt

using Enzyme
using Interpolations
using SparseArrays
using TMI

using Enzyme: Annotation
using TMI: BoundaryCondition, Field, Grid, MassFraction, Source
using TMI: adjustboundarycondition, adjustsource!, cartesianindex, gsetboundarycondition, interpweights, observe, steadyinversion, step_cartesian, unvec, unvec!, watermassmatrix, wet
import TMI: interpweights

import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth
import Enzyme.EnzymeRules: augmented_primal, needs_primal, needs_shadow, reverse

include("enzyme/gunvec_inplace.jl")
include("enzyme/gvec.jl")
include("enzyme/gwatermassmatrix.jl")
include("enzyme/gldiv_field.jl")
include("enzyme/accumulate_into.jl")
include("enzyme/gobserve.jl")
include("enzyme/gadjustsource.jl")

end
