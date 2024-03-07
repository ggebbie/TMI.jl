import Pkg; Pkg.activate(".")
using Revise
using TMI
#using LinearAlgebra

fname = TMI.pkgdatadir("simple_MITgcm.nc")
println(fname)

γ =  Grid(fname, "maskC", "XC", "YC", "Z")
y = (θ = readfield(fname,"THETA",γ),
    S = readfield(fname,"SALT",γ))

# m̃ = massfractions(y)

# pull this function apart
c = y 
γ = first(c).γ
Rfull = CartesianIndices(γ.wet)
Ifirst, Ilast = first(Rfull), last(Rfull)

# loop over interior points (dirichlet bc at surface)
R = copy(γ.R)
Iint = cartesianindex(γ.interior)

@time m0 = TMI.massfractions_isotropic(γ)
@time n   = TMI.neighbors(m0,γ)
@time n2 = TMI.neighbors(γ)  # faster than the previous two put together

# allocate maximum needed
nmax = maximum(n)
nrow = length(c) + 1
b = vcat(zeros(nrow-1),1.0)
mlocal = zeros(nmax)
nlocal = zeros(nrow)
ϵ = 1e-8 # for checking tolerances

I = CartesianIndex(12,8,2)
I = CartesianIndex(13,8,2)
I = CartesianIndex(11,10,2)

ncol = n.tracer[I]
m0local = ones(ncol) ./ ncol
Alocal, single_connection = TMI.local_watermass_matrix_old(c,m0,I,n)
n0 = b - Alocal*m0local


 if sum(abs.(n0[1:nrow])) < ϵ # something small
     mlocal[1:ncol] = m0local
     #println("first guess good enough @ ",I)
 else
     mlocal[1:ncol] = m0local + Alocal\n0
     nlocal[1:nrow] = b - Alocal*mlocal[1:ncol]

     # check fit and check non-negativity
     if single_connection ||
         ((sum(abs.(nlocal)) > ϵ) || !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ))

         # quadratic programming
         println("run local quadprog @ ",I)
         model, x = TMI.local_quadprog(m0local, Alocal, b)

         if termination_status(model) != OPTIMAL
             #     model2, x2 = local_quadprog(Alocal,
             #         m0local, wlocal)
             #     if termination_status(model2) != OPTIMAL
             println("WARNING: NEVER FOUND A BETTER SOLUTION! @ ",I)
             mlocal[1:ncol] = m0local
             #     else
             #         println("satisfied only mass cons @ ",I)
             #         mlocal[1:ncol] = value.(x2)
             #     end
         else
             #println("solved for both perfect data and mass cons @ ",I)
             mlocal[1:ncol] = value.(x)
         end
     end
 end













if sum(abs.(n0)) < ϵ # something small
   mlocal[1:ncol] = m0
else

    # Invert! by maximizing mixing and fitting tracers/mass perfectly
    # attempts to fit tracers and mass conservation perfectly
    #m_local[1:nlocal] = x0 + Alocal'*((Alocal*Alocal')\noise)
    Alocal2 = vcat(Alocal,ones(1,size(Alocal,2)))
    mlocal[1:ncol] = m0local + Alocal2\vcat(n0,0.0)
    nlocal[1:nrow] = vcat(Alocal*mlocal[1:ncol],1-sum(mlocal[1:ncol]))



        # well-mixed first guess
        ncol = n.tracer[I]
        m0local = ones(ncol) ./ ncol
        Alocal, single_connection = local_watermass_matrix_old(c,m,I,n)
        n0 = b - Alocal*m0local

        if sum(abs.(n0[1:nrow])) < ϵ # something small
            mlocal[1:ncol] = m0local
            println("first guess good enough @ ",I)
        else
            mlocal[1:ncol] = m0local + Alocal\n0
            nlocal[1:nrow] = b - Alocal*mlocal[1:ncol]

            # check fit and check non-negativity
            if single_connection ||
                ((sum(abs.(nlocal)) > ϵ) || !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ))

                # quadratic programming
                println("run local quadprog @ ",I)
                model, x = local_quadprog(m0, Alocal, b)

                if termination_status(model) != OPTIMAL
                #     model2, x2 = local_quadprog(Alocal,
                #         m0local, wlocal)
                #     if termination_status(model2) != OPTIMAL
                        println("WARNING: NEVER FOUND A BETTER SOLUTION! @ ",I)
                        mlocal[1:ncol] = m0local
                #     else
                #         println("satisfied only mass cons @ ",I)
                #         mlocal[1:ncol] = value.(x2)
                #     end
                else
                    println("solved for both perfect data and mass cons @ ",I)
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

    #########end
    if single_connection ||
        ((sum(abs.(nlocal)) > ϵ) ||
            !(1.0 - ϵ < sum(abs.(mlocal[1:ncol])) < 1 + ϵ ))

        @time out = TMI.local_quadprog(Alocal,m0local)

        if isnothing(out)
            # relax tracer assumption
            w = [1e-2,1e-3]
            @time out = TMI.local_quadprog(Alocal,m0local,w)
            #sum(out)
            #TMI.local_objective_obs(Alocal,m0local,w)
            #TMI.local_objective_obs(Alocal,out,w)
        end

    end
