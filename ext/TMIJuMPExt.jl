module TMIJuMPExt

using TMI 
using JuMP

#= HiGHS vs. COSMO
HiGHS is more accurate when fitting both tracers and mass.
Computational time is the same.
COSMO is more robust when fitting mass alone (imperfect fit of tracers).
COSMO should be improved by warm starting and making tolerances more strict.
=#
#using HiGHS
using COSMO


"""
function massfractions(c::NamedTuple, w::NamedTuple; alg = :local)

Create NamedTuple of mass fractions from observations `c`

Doesn't produce a `MassFraction` struct and thus
is named in lower case.

# Arguments
- `c::NamedTuple`: input observations
- `w::NamedTuple`: scale size of observations (used if obs not fit exactly)
- `alg=:local`: default algorithm is `:local`
# Output
- `m::NamedTuple`: collection of mass fractions
"""
function massfractions(c::NamedTuple, w::NamedTuple; alg = :local) # assume: `Field`s inside
    if alg == :local
        return local_solve(c::NamedTuple, w::NamedTuple)
    else
        println("TMI.massfractions: global solution for TMI water-mass matrix not implemented yet")
    end
end

"""
function local_solve(c::NamedTuple, w::NamedTuple)

Given tracers, solve for mass fractions using
a repeated local algorithm.

`local_solve!` is the workhorse for this algorithm and
specifies a default inverse tapering parameter
of `α=100_000`.

# Arguments
- `c::NamedTuple`: `Field` tracers (observations or model output)
- `w::NamedTuple`: scale size of tracers (for weighting observational fits)
# Output
- `m̃::NamedTuple`: `MassFraction` in a NamedTuple collection
"""
function local_solve(c::NamedTuple, w::NamedTuple) # assume: `Field`s inside
    γ = first(c).γ # assumption: all grids match up
    m̃ = TMI.massfractions_isotropic(γ) # good first guess
    TMI.local_solve!(m̃,c,w)
    return m̃
end

function local_quadprog(m0local::Vector, Alocal::Matrix, b::Vector)
    #model =
    #     Model(optimizer_with_attributes(COSMO.Optimizer, "verbose" => false,
    #         "eps_abs" => 1e-4, "eps_rel" = 1e-4))

    model = Model(COSMO.Optimizer);
    #model = Model(HiGHS.Optimizer);
    set_time_limit_sec(model, 0.1)
    set_silent(model);
    ncol = length(m0local)
    @variable(model, 0.0 <= x[i = 1:ncol] <= 1.0 ) 
    @constraint(model, Alocal*x == b)
    #@constraint(model, sum(x) == 1.0)
    @objective(model, Min,
        local_objective_mixing(x,m0local))
    optimize!(model)
    return model, x 
end

local_objective_mixing(x,m0local) =
    sum((x.-m0local).^2)
    
function local_quadprog(Alocal::Matrix,m0local::Vector,w::Vector)
    #settings = COSMO.Settings(verbose = false, eps_abs = 1e-2, eps_rel = 1e-2)
    model =
         Model(optimizer_with_attributes(COSMO.Optimizer, "verbose" => false,
             "eps_abs" => 1e-4, "eps_rel" => 1e-4))
    set_time_limit_sec(model, 0.1)

    #model = Model(COSMO.Optimizer);
    #set_silent(model);
    ncol = length(m0local) 
    @variable(model, 0.0 <= x[i = 1:ncol] <= 1.0 ) 
    @constraint(model, sum(x) == 1.0)
    #@objective(model, Min,
    #    sum(((Alocal[1:end-1,:]*x)./w).^2))
    #    println(local_objective_obs(Alocal,x,w))

    # cut off mass conservation from this objective function
    @objective(model, Min,
        local_objective_obs(Alocal[1:end-1,:],x,w[1:end-1]))
    optimize!(model)
    return model, x
end

local_objective_obs(Alocal,x,w) =
    sum(((Alocal*x)./w).^2)
