# TMI.jl Enzyme Rules

This directory contains custom reverse-mode rules used by `TMIEnzymeExt`.

Why these rules exist:
- They make Enzyme run on the current Enzyme tests and gradient scripts in
  this repo.
- Generic Enzyme tracing can fail on some mutating and sparse-operator paths.
- Additional rules may be needed for TMI operators that are not covered by the
  current tests/scripts.

Files:
- `gunvec_inplace.jl`: reverse rules linked to `TMI.unvec!` on tracer-like controls and
  heterogeneous `NamedTuple` control bundles.
- `gvec.jl`: reverse rules linked to `TMI.vec` so packed control vectors map cleanly back
  to TMI structs (`Field`, `BoundaryCondition`, `Source`, `MassFraction`).
- `gwatermassmatrix.jl`: reverse rule linked to `TMI.watermassmatrix(m, γ)`.
- `gldiv_field.jl`: reverse rules linked to `TMI.\\(A, d::Field)` with sparse/LU
  support and `dJ/dA` accumulation.
- `gobserve.jl`: reverse rules linked to `TMI.observe`.
- `gadjustsource.jl`: reverse rule linked to mutating `TMI.adjustsource!`.
- `accumulate_into.jl`: Enzyme accumulation hooks used by TMI overloaded
  arithmetic on `Field`/`BoundaryCondition`/`Source`/`MassFraction`.

If you are new to Enzyme in this repo, start from `gvec.jl` and
`gunvec_inplace.jl`, then read `gldiv_field.jl` for sparse-solve adjoints.
