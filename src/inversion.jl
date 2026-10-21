"""
    struct Observations

Tracer observations for `Inversion`; constructor in the `TMIEnzymeOptimizationExt` extension.

- `values::NamedTuple`: observed values of each tracer at its observation points
- `σ::NamedTuple`: standard deviations of the observed values
- `locs::NamedTuple`: (lon, lat, depth) of each tracer's observations, `nothing` if gridded
- `wis::NamedTuple`: interpolation weights at `locs`, `nothing` if gridded
- `yvec::Vector`: all observed values, also `vec(y)`
- `γ::Grid`: TMI grid
- `Wⁱ::SparseMatrixCSC`: inverse observational error covariance
"""
struct Observations{Y,S,X,I,V,G,M}
    values::Y
    σ::S
    locs::X
    wis::I
    yvec::V
    γ::G
    Wⁱ::M
end

"""
    struct Inversion

Inverse problem for `runinversion`; constructor in the `TMIEnzymeOptimizationExt` extension.

- `y::Observations`: observations
- `b₀`, `q₀`, `m₀`: first-guess boundary conditions, sources, and mass fractions
- `r::NamedTuple`: stoichiometric ratios of each tracer to the sources
- `u::NamedTuple`: preallocated controls
- `ranges::NamedTuple`: blocks `b`, `q`, and `m` of the control vector `uvec`
- `σ`, `lb`, `ub`: standard deviation and limits of each element of `uvec`
- `smoothness::SparseMatrixCSC`: weights of the prior smoothness cost, `∇²ᵀ Diagonal(L⁴/σ²) ∇²`
  for each boundary adjustment with a smoothing length `L`, zero elsewhere in `uvec`
"""
struct Inversion{Y,B,Q,R,M,U,N,V,W}
    y::Y
    b₀::B
    q₀::Q
    r::R
    m₀::M
    u::U
    ranges::N
    σ::V
    lb::V
    ub::V
    smoothness::W
end

"""
    function runinversion(inversion, optimizer; name, iterations, checkpoint_interval,
        checkpoint_directory, data_directory, number_of_gradient_checks)

Solve an `Inversion` with Ipopt. Defined in the `TMIEnzymeOptimizationExt`
extension, loaded with `using Enzyme, Optimization, OptimizationIpopt`.
"""
function runinversion end
