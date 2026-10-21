#=
Inversions of TMI tracers for surface boundary conditions, interior sources, and
mass fractions, with gradients from Enzyme and optimization by Ipopt. Loaded by
`using Enzyme, Optimization, OptimizationIpopt` together with TMI.

A run, as in scripts/full_inversion.jl:
1. `Observations(observed, σ, γ; ...)` in observations.jl: the observed values,
   their σ, and the weighting matrix Wⁱ (built in weighting_matrices.jl).
2. `Inversion(y; b₀, q₀, r, m₀, controls)` in inverse_problem.jl: the first
   guesses, and the σ, limits, and smoothness prior of the adjusted controls.
3. `runinversion(inversion, optimizer; ...)` in this file: Ipopt minimizes
   `costfunction`, which runs the forward model (adjust the first guesses, build
   the water-mass matrix, solve for the steady-state tracers) and adds the data
   cost and prior costs (`costterms`). Gradients come from `costgradient!`.

Enzyme differentiates all of this itself, except four functions whose rules are
in enzyme_rules.jl. The `Observations` and `Inversion` structs and the
`runinversion` stub are defined in src/inversion.jl, since an extension cannot
export names.
=#
module TMIEnzymeOptimizationExt

using Enzyme, LinearAlgebra, NCDatasets, Optimization, SparseArrays, TMI
using OptimizationIpopt: IpoptOptimizer

# TMI controls exceed Enzyme's default type-analysis lattice.
function __init__()
    Enzyme.API.maxtypeoffset!(8192)
    Enzyme.API.maxtypedepth!(12)
end

using Enzyme: Annotation
using TMI: step_cartesian
import TMI: Inversion, Observations, gradient_check, observe, runinversion
import Enzyme.EnzymeRules: AugmentedReturn, RevConfigWidth, augmented_primal, needs_primal, needs_shadow, reverse

"""
    function foreachtracer(f, names)

The result of `f` for each tracer (or source) in `names`, as a `NamedTuple`
keyed by name. Used throughout the extension.

# Arguments
- `f`: function of a name
- `names`: tuple of names
# Output
- `NamedTuple` of `f(name)` keyed by `names`
"""
foreachtracer(f, names) = NamedTuple{names}(map(f, names))

include("weighting_matrices.jl")
include("observations.jl")
include("inverse_problem.jl")
include("enzyme_rules.jl")

"""
    function runinversion(inversion, optimizer; name, iterations, checkpoint_interval,
        checkpoint_directory, data_directory, number_of_gradient_checks=5)

Check the gradient, then minimize `J` with Ipopt subject to the limits on
`uvec` and mass conservation. Ipopt works on `z = uvec ./ σ`, starting at the
first guess `z = 0`. The observations go to
`checkpoint_directory/name_observations.nc` (`writeobservations`), every
`checkpoint_interval`-th iteration to
`checkpoint_directory/name_iteration_NNNNNN.nc`, and the result to
`data_directory/name_final.nc` (`writecheckpoint`).

# Arguments
- `inversion::Inversion`: inverse problem
- `optimizer::IpoptOptimizer`: Ipopt settings; the caller's copy is not modified
- `name`: experiment name for output files
- `iterations`: maximum number of Ipopt iterations
- `checkpoint_interval`: iterations between saved checkpoints
# Optional Arguments
- `checkpoint_directory`, `data_directory`: output directories
- `number_of_gradient_checks`: controls checked in each block before the run; 0 skips the check
# Output
- `file`: NetCDF file with the final iteration
"""
function runinversion(inversion::Inversion, optimizer::IpoptOptimizer; name::AbstractString,
    iterations::Integer, checkpoint_interval::Integer,
    checkpoint_directory::AbstractString=TMI.pkgdir("checkpoints", name),
    data_directory::AbstractString=TMI.pkgdatadir(), number_of_gradient_checks::Integer=5)
    σ = inversion.σ
    workspace = InversionWorkspace(inversion)
    uvec = similar(σ)
    # Ipopt sees the controls divided by their σ, z = uvec ./ σ. It asks for the
    # cost alone at its trial points, most of which it rejects, and for the
    # gradient only at the points it accepts.
    J(z, _=nothing) = (uvec .= σ .* z; costfunction(uvec, inversion, workspace.u, workspace.F))
    gJ!(gz, z, _=nothing) = (uvec .= σ .* z; costgradient!(gz, workspace, uvec); gz .*= σ; nothing)
    z₀ = zero(σ)
    (; cons, cons_j, cons_jac_prototype, lcons, ucons) = massconservation(inversion)
    # Optimization.jl requires an AD type; the gradient comes from gJ! and the Jacobian from cons_j
    f = OptimizationFunction(J, AutoEnzyme(); grad=gJ!, cons, cons_j, cons_jac_prototype)
    problem = OptimizationProblem(f, z₀; lb=inversion.lb ./ σ, ub=inversion.ub ./ σ, lcons, ucons)
    number_of_gradient_checks > 0 && gradient_check(inversion; n=number_of_gradient_checks, workspace)

    # write the observations and the first guess, then every checkpoint_interval-th iteration
    checkpoint(k) = joinpath(checkpoint_directory, "$(name)_iteration_$(lpad(k, 6, '0')).nc")
    writeobservations(joinpath(checkpoint_directory, "$(name)_observations.nc"), inversion)
    writecheckpoint(checkpoint(0), workspace, σ .* z₀; iteration=0)
    callback = (state, _) -> begin
        state.iter > 0 && state.iter % checkpoint_interval == 0 &&
            writecheckpoint(checkpoint(state.iter), workspace, σ .* state.u; iteration=state.iter)
        false
    end
    optimizer = deepcopy(optimizer)
    logfile = joinpath(checkpoint_directory, "ipopt_output.txt")
    optimizer.additional_options["output_file"] = logfile
    println("Starting Ipopt; writing log to $logfile")
    solution = solve(problem, optimizer; maxiters=iterations, verbose=5, callback)
    return writecheckpoint(joinpath(data_directory, "$(name)_final.nc"), workspace, σ .* solution.u;
        iteration=solution.stats.iterations)
end

# METIS ordering: about 40% less fill than UMFPACK's default AMD ordering on TMI
# grids; `lu!` keeps the ordering of the factorization it refactors
const METIS_CONTROL = let control = SparseArrays.UMFPACK.get_umfpack_control(Float64, Int64)
    control[SparseArrays.UMFPACK.JL_UMFPACK_ORDERING] =
        SparseArrays.UMFPACK.LibSuiteSparse.UMFPACK_ORDERING_METIS
    control
end

"""
    struct InversionWorkspace

Storage that the cost function and its Enzyme gradient reuse at every control vector

- `inversion::Inversion`: inverse problem
- `u`: controls, a copy of `inversion.u` that `adjustfirstguess!` overwrites
- `du`: Enzyme shadow of `u`
- `F`: UMFPACK factorization of the water-mass matrix (METIS ordering), refactored
  in place by `lu!` at each control vector
- `dF`: Enzyme shadow of `F`
"""
struct InversionWorkspace{I,U,L}
    inversion::I
    u::U
    du::U
    F::L
    dF::L
end

"""
    function InversionWorkspace(inversion)

Workspace with zero shadows and the factorization of the first-guess water-mass matrix.

# Arguments
- `inversion::Inversion`: inverse problem
# Output
- `InversionWorkspace`
"""
function InversionWorkspace(inversion::Inversion)
    u = deepcopy(inversion.u)
    F = lu(watermassmatrix(inversion.m₀, inversion.y.γ); control=METIS_CONTROL)
    return InversionWorkspace(inversion, u, Enzyme.make_zero(u), F, Enzyme.make_zero(F))
end

"""
    function costfunction(uvec, inversion, u, F)

Cost function `J = Jdata + Jb + Js + Jq + Jm` of the control vector `uvec`: the
data cost plus the prior costs, from `costterms`. Enzyme differentiates this function.

# Arguments
- `uvec`: control vector
- `inversion::Inversion`: observations, first guesses, and weighting matrices
- `u`: preallocated controls, a copy of `inversion.u`; Enzyme differentiates
  through it
- `F`: UMFPACK factorization refactored in place, see `forwardmodel`
# Output
- `J`: cost function
"""
costfunction(uvec::AbstractVector, inversion::Inversion, u, F) =
    sum(costterms(forwardmodel(uvec, inversion, u, F).c, uvec, inversion))

"""
    function costterms(c, uvec, inversion)

The data cost and the prior costs that add up to the cost function:
- `Jdata = nᵀ Wⁱ n`: data cost, with `n = ỹ - y` the modeled minus the observed
  value at each observation point and `Wⁱ` the inverse covariance of the
  observational errors
- `Jb`, `Jq`, `Jm`: prior costs of the adjustments to the boundary conditions,
  sources, and mass fractions, `Σ (u/σ)²` over their blocks of `uvec`
- `Js`: prior smoothness cost of the boundary adjustments, `uvecᵀ S uvec` with
  `S = ∇²ᵀ Diagonal(L⁴/σ²) ∇²` (`inversion.smoothness`) for each boundary
  condition with a smoothing length `L`

# Arguments
- `c::NamedTuple`: steady-state tracers
- `uvec`: control vector
- `inversion::Inversion`: observations, standard deviations `σ`, and `smoothness`
# Output
- `(; Jdata, Jb, Js, Jq, Jm)`
"""
# @inline: as a separate call, Enzyme's gradient through it allocates 30 MiB more on WOCE
@inline function costterms(c::NamedTuple, uvec::AbstractVector, inversion::Inversion)
    (; y, σ, ranges, smoothness) = inversion
    # model minus observations at the observation points
    n = observe(c, y) - vec(y)
    # prior cost of the adjustments in a block of uvec, each with standard deviation σ
    priorcost(block) = sum(i -> (uvec[i] / σ[i])^2, block; init=0.0)
    Jdata = dot(n, y.Wⁱ, n)
    Jb = priorcost(ranges.b)
    Js = dot(uvec, smoothness, uvec)
    Jq = priorcost(ranges.q)
    Jm = priorcost(ranges.m)
    return (; Jdata, Jb, Js, Jq, Jm)
end

"""
    function forwardmodel(uvec, inversion, u, F)

The forward model for a control vector: the first guesses adjusted by `uvec`,
their water-mass matrix, and the steady-state tracers.

# Arguments
- `uvec`: control vector
- `inversion::Inversion`: first guesses
- `u`: preallocated controls, a copy of `inversion.u`
- `F`: UMFPACK factorization of a water-mass matrix on the same grid, refactored
  in place by `lu!`, which keeps its ordering
# Output
- `(; b, q, m, A, c)`: adjusted boundary conditions, sources, and mass
  fractions, water-mass matrix `A`, and steady-state tracers `c` keyed like
  `inversion.b₀`
"""
function forwardmodel(uvec::AbstractVector, inversion::Inversion, u, F)
    γ = inversion.y.γ
    # first guesses adjusted by the controls
    b, q, m = adjustfirstguess!(u, inversion, uvec)
    # water-mass matrix of the adjusted mass fractions, and the steady-state tracers
    A = watermassmatrix(m, γ)
    c = steadytracers(lu!(F, A), b, q, inversion.r, γ)
    return (; b, q, m, A, c)
end

"""
    function steadytracers(Alu, b, q, r, γ)

Steady-state tracers for every boundary condition in `b`, one tracer per
boundary condition (`steadyinversion(Alu, b::NamedTuple, γ; q, r)` instead
combines several boundaries of one tracer). The interior source of each tracer
is the sum of the sources in `q` scaled by its stoichiometric ratios in `r`.

# Arguments
- `Alu`: LU factorization of `A = watermassmatrix(m, γ)`
- `b::NamedTuple`: surface boundary conditions
- `q::NamedTuple`: interior sources
- `r::NamedTuple`: stoichiometric ratios, signed as `r` in `steadyinversion`,
  e.g. `(NO₃ = (qPO₄ = 15.5,),)`; a tracer missing from `r` has no interior source
- `γ::Grid`: TMI grid
# Output
- `c::NamedTuple`: steady-state `Field`s keyed like `b`
"""
function steadytracers(Alu, b::NamedTuple, q::NamedTuple, r::NamedTuple, γ::Grid)
    return foreachtracer(keys(b)) do k
        rₖ = get(r, k, nothing)
        qₖ = isnothing(rₖ) ? nothing : reduce(+, map(j -> rₖ[j] * q[j], keys(rₖ)))
        # r = -1: watermassmatrix(m, γ) has +1 on the interior diagonal, TMI files -1
        steadyinversion(Alu, b[k], γ; q=qₖ, r=-1.0)
    end
end

"""
    function costgradient!(g, workspace, uvec)

Cost function at `uvec`, with its gradient written to `g`, from one Enzyme
reverse pass through `costfunction`.

# Arguments
- `g`: gradient, overwritten
- `workspace::InversionWorkspace`
- `uvec`: control vector
# Output
- `J`: cost function
"""
function costgradient!(g, workspace::InversionWorkspace, uvec)
    fill!(g, 0)
    Enzyme.remake_zero!(workspace.du)
    return Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        Enzyme.Const(costfunction), Enzyme.Active, Enzyme.Duplicated(uvec, g),
        Enzyme.Const(workspace.inversion), Enzyme.Duplicated(workspace.u, workspace.du),
        Enzyme.Duplicated(workspace.F, workspace.dF))[2]
end

"""
    function gradient_check(inversion; n=5, rtol=1e-4, atol=1e-7, ε=1e-4, workspace)

Method of TMI's `gradient_check` for an `Inversion`. Compare the Enzyme gradient
at `uvec = 0` with centered finite differences at up to `n` controls in each of
the `b`, `q`, and `m` blocks, at the midpoints of `n` equal segments of the
block. Each block is compared on its own, so large derivatives in one block
cannot hide errors in another. The finite difference is computed from the
sampled tracers `ỹ`, `J(δ) - J(-δ) = (ỹ₊ - ỹ₋)ᵀ Wⁱ (ỹ₊ + ỹ₋ - 2y)`, which avoids
subtracting two large values of `J`; the prior costs are even in `uvec` and cancel.

# Arguments
- `inversion::Inversion`: inverse problem
# Optional Arguments
- `n`: number of controls checked in each block
- `rtol`, `atol`: tolerances of the comparison within each block
- `ε`: finite-difference step
- `workspace`: `InversionWorkspace` to reuse
# Output
- `NamedTuple` keyed by block of `(indices, ∇J, ∇J_finite)`; errors if a block disagrees
"""
function gradient_check(inversion::Inversion; n::Integer=5, rtol=1e-4, atol=1e-7,
    ε=1e-4, workspace=InversionWorkspace(inversion))
    y = inversion.y
    uvec = zero(inversion.σ)
    g = similar(uvec)
    costgradient!(g, workspace, uvec)
    checks = map(inversion.ranges) do r
        k = min(n, length(r))
        indices = isempty(r) ? Int[] : r[ceil.(Int, (length(r) / k) .* ((1:k) .- 0.5))]
        ∇J_finite = map(indices) do i
            δ = zero(uvec)
            δ[i] = ε
            ỹ₊ = observe(forwardmodel(uvec + δ, inversion, workspace.u, workspace.F).c, y)
            ỹ₋ = observe(forwardmodel(uvec - δ, inversion, workspace.u, workspace.F).c, y)
            dot(ỹ₊ - ỹ₋, y.Wⁱ, ỹ₊ + ỹ₋ - 2vec(y)) / 2ε
        end
        (; indices, ∇J=g[indices], ∇J_finite)
    end
    for (block, check) in pairs(checks)
        isempty(check.indices) && continue
        for (i, g, g_finite) in zip(check.indices, check.∇J, check.∇J_finite)
            println("  $block control $i: Enzyme $g, finite difference $g_finite")
        end
        difference = norm(check.∇J - check.∇J_finite)
        scale = max(norm(check.∇J), norm(check.∇J_finite))
        println("  $block relative difference ", iszero(scale) ? difference : difference / scale)
        isapprox(check.∇J, check.∇J_finite; rtol, atol) ||
            error("gradient check failed for $block controls $(check.indices)")
    end
    return checks
end

"""
    function massconservation(inversion)

Equality constraints `Σ m = 1` in every interior cell, written for the scaled
controls `z = uvec ./ σ` used by Ipopt. The constraints are linear in `z`,
`Σ m = P m₀ + jacobian z`, so the Jacobian is constant.

# Arguments
- `inversion::Inversion`: inverse problem
# Output
- the Optimization.jl arguments `cons`, `cons_j`, `cons_jac_prototype`, `lcons`,
  `ucons`; all `nothing`, Optimization.jl's defaults, when mass fractions are
  not adjusted
"""
function massconservation(inversion::Inversion)
    ranges = inversion.ranges
    isempty(ranges.m) && return (cons=nothing, cons_j=nothing, cons_jac_prototype=nothing,
        lcons=nothing, ucons=nothing)
    γ = inversion.y.γ
    R = linearindex(γ.interior)
    rows = reduce(vcat, [R[wet(m)] for m in inversion.m₀])
    P = sparse(rows, 1:length(rows), 1.0, sum(γ.interior), length(rows))
    Pm₀ = P * vec(inversion.m₀)
    jacobian = sparse(rows, ranges.m, inversion.σ[ranges.m], sum(γ.interior),
        length(inversion.σ))
    cons = (Σm, z, _) -> (mul!(Σm, jacobian, z); Σm .+= Pm₀; nothing)
    cons_j = (out, _, _) -> (copyto!(nonzeros(out), nonzeros(jacobian)); nothing)
    return (cons=cons, cons_j=cons_j, cons_jac_prototype=jacobian,
        lcons=ones(size(jacobian, 1)), ucons=ones(size(jacobian, 1)))
end

"""
    function writeobservations(file, inversion)

Write the observations an inversion fits to the NetCDF `file`: for each tracer
`k`, its observed values `k` and their uncertainty `σk`, at the interior cells
(in `γ.interior` order) for gridded tracers, and with the coordinates
`k_longitude`, `k_latitude`, and `k_depth` for tracers observed at locations.

# Arguments
- `file`: NetCDF file name, replaced if it exists
- `inversion::Inversion`: holds the observations
# Output
- `file`
"""
function writeobservations(file, inversion::Inversion)
    y = inversion.y
    mkpath(dirname(file))
    rm(file; force=true)
    NCDataset(file, "c") do ds
        for k in keys(y.values)
            locs = y.locs[k]
            longname, units = inversion.b₀[k].longname, inversion.b₀[k].units
            variables = [(String(k), y.values[k], longname, units),
                ("σ$k", y.σ[k], "1σ uncertainty in $longname", units)]
            isnothing(locs) || append!(variables, [
                ("$(k)_longitude", first.(locs), "longitude", "°E"),
                ("$(k)_latitude", getindex.(locs, 2), "latitude", "°N"),
                ("$(k)_depth", last.(locs), "depth", "m")])
            observation = "$(k)_observation"
            defDim(ds, observation, length(y.values[k]))
            for (vname, data, vlongname, vunits) in variables
                defVar(ds, vname, Float64.(data), (observation,);
                    attrib=Dict("longname" => vlongname, "units" => vunits))
            end
        end
    end
    println("Saved observations to $file")
    return file
end

"""
    function writecheckpoint(file, workspace, uvec; iteration)

Write the state of the inversion at the control vector `uvec` to the NetCDF
`file`: the modeled tracers, the sources, the water-mass matrix, and the
attributes `iteration`, `J`, its data cost `Jdata` and prior costs `Jb`, `Js`,
`Jq`, and `Jm` (see `costterms`), and their sum `Jcontrol = Jb + Js + Jq + Jm`.

# Arguments
- `file`: NetCDF file name, replaced if it exists
- `workspace::InversionWorkspace`: the inverse problem, and the controls `u` and
  factorization `F` that are overwritten here
- `uvec`: control vector
# Optional Arguments
- `iteration`: optimizer iteration
# Output
- `file`
"""
function writecheckpoint(file, workspace::InversionWorkspace, uvec; iteration::Integer)
    (; inversion, u, F) = workspace
    γ = inversion.y.γ
    (; q, A, c) = forwardmodel(uvec, inversion, u, F)
    terms = costterms(c, uvec, inversion)
    mkpath(dirname(file))
    rm(file; force=true)
    foreach(k -> writefield(file, Field(c[k].tracer, γ, k, c[k].longname, c[k].units)), keys(c))
    foreach(qₖ -> TMI.writesource(file, qₖ), q)
    TMI.watermassmatrix2nc("", A; filenetcdf=file)
    NCDataset(file, "a") do ds
        ds.attrib["iteration"] = iteration
        ds.attrib["J"] = sum(terms)
        foreach(name -> ds.attrib[String(name)] = terms[name], keys(terms))
        ds.attrib["Jcontrol"] = terms.Jb + terms.Js + terms.Jq + terms.Jm
    end
    println("Saved iteration $iteration to $file")
    return file
end

end
