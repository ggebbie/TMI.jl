#using TMI: massfractions_isotropic
import Pkg; Pkg.activate(".")
using Revise
using TMI
using LinearAlgebra

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

#nrow   = length(c) + 1 # add mass conservation
# allocate maximum needed
nmax = maximum(n)
#l = zeros(nmax)
#u = ones(nmax)
#b = vcat(zeros(nrow-1),1.0)
mlocal = zeros(nmax)
nlocal = zeros(nrow)
ϵ = 1e-8 # for checking tolerances

I = CartesianIndex(12,8,2)
I = CartesianIndex(13,8,2)

# well-mixed first guess
ncol = n.tracer[I]
m0local = ones(ncol) ./ ncol

# does it already fit the data?
        #Alocal = local_watermass_matrix(c,m,I,ncol.tracer[I])
Alocal, single_connection = TMI.local_watermass_matrix(c,m0, I, n)
n0 = Alocal*m0local

if sum(abs.(n0)) < ϵ # something small
   mlocal[1:ncol] = m0
else

    # Invert! by maximizing mixing and fitting tracers/mass perfectly
    # attempts to fit tracers and mass conservation perfectly
    #m_local[1:nlocal] = x0 + Alocal'*((Alocal*Alocal')\noise)
    Alocal2 = vcat(Alocal,ones(1,size(Alocal,2)))
    mlocal[1:ncol] = m0local + Alocal2\vcat(n0,0.0)
    nlocal[1:nrow] = vcat(Alocal*mlocal[1:ncol],1-sum(mlocal[1:ncol]))

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
