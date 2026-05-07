# TMI.jl Enzyme Rules

This directory contains custom reverse-mode rules for Enzyme.jl.

The files are named with a `g` prefix because they implement adjoint, or gradient,
propagation. In rule code, `ginput` means `dJ/dinput`, where `J` is a scalar cost
function.

- `gunvec_boundarycondition.jl`: vector-to-`BoundaryCondition` scatter/gather.
- `gunvec_field.jl`: vector-to-`Field` scatter/gather.
- `gwatermassmatrix.jl`: mass fractions to sparse water-mass matrix.
- `gldiv_field.jl`: sparse or LU linear solve with a `Field` right-hand side.
- `gobserve.jl`: point observations and their interpolation adjoint.

`enzyme_rules.jl` collects these methods into `TMIEnzymeRules`.
