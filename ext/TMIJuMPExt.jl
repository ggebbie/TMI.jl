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

import TMI: local_solve, local_quadprog

export local_solve, local_quadprog 

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
    local_solve!(m̃,c,w)
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

local_objective_obs(Alocal,x,w) = sum(((Alocal*x)./w).^2)


function local_solve!(m::Union{NamedTuple,Vector}, c::NamedTuple, w::NamedTuple ; alg = :quadprog)

    γ = first(c).γ
    #Rfull = CartesianIndices(γ.wet)
    Iint = cartesianindex(γ.interior)

    n   = TMI.neighbors(m,γ)
    nrow   = length(c) + 1 # add mass conservation

    # allocate maximum needed
    nmax = maximum(n)
    b = vcat(zeros(nrow-1),1.0)
    mlocal = zeros(nmax)
    nlocal = zeros(nrow)
    ϵ = 1e-6 # tolerance of mass conservation
    Jlimit = 1e-4 # scaled cost function target
    wlocal = vcat([w1 for w1 in w],ϵ)

    for I in Iint

        # well-mixed first guess
        ncol = n.tracer[I]
        m0local = ones(ncol) ./ ncol
        Alocal, single_connection = TMI.local_watermass_matrix(c,m,I,n)
        n0 = b - Alocal*m0local
        J0 = sum((n0./wlocal).^2)/nrow

        if J0 < Jlimit # hit target already?
            mlocal[1:ncol] = m0local
            #println("first guess good enough @ ",I)
        else
            mlocal[1:ncol] = m0local + Alocal\n0
            nlocal[1:nrow] = b - Alocal*mlocal[1:ncol]
            Jlocal = sum((nlocal./wlocal).^2)/nrow

            # check connection, fit, and non-negativity
            if single_connection ||
                ((Jlocal > Jlimit) ||
                    !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ))

                # quadratic programming
                # println("run local quadprog @ ",I)
                model, x = local_quadprog(m0local, Alocal, b)

                if termination_status(model) != OPTIMAL
                    model2, x2 = local_quadprog(Alocal,
                        m0local, wlocal)
                    if termination_status(model2) != OPTIMAL
                        println("TMI.local solve: WARNING: NEVER FOUND A BETTER SOLUTION! @ ",I)
                        mlocal[1:ncol] = m0local
                    else
                        # println("satisfied only mass cons @ ",I)
                        mlocal[1:ncol] = value.(x2)
                    end
                else
                    #println("solved for both perfect data and mass cons @ ",I)
                    mlocal[1:ncol] = value.(x)
                end
            end
        end
        
        # Save into proper mass fractions
        i = 0
        for m1 in m
            if m1.γ.wet[I]
                i += 1
                # should not need to check bounds 
                m1.fraction[I] = mlocal[i]
            end
        end
    end
end

end # module 
