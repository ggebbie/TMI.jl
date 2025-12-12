module TMI

using LinearAlgebra
using SparseArrays
using NetCDF
using Downloads
using GoogleDrive
using Distances
using GibbsSeaWater  
using Distributions
using Optim
using Interpolations
using LineSearches
using NCDatasets
using UnicodePlots
using OrderedCollections
using Downloads
using LinearSolve
using IncompleteLU

export config, download_file,
    surfaceindex, lonindex, latindex, depthindex,
    surfacepatch, 
    layerthickness, cellarea, cellvolume,
    planview, 
    section,
    tracerinit,
    watermassmatrix, watermassdistribution,
    circulationmatrix, boundarymatrix,
    dirichletmatrix, mixedlayermatrix, 
    linearindex, nearestneighbor, updatelinearindex,
    nearestneighbormask, horizontaldistance,
    readtracer, readfield, writefield,
    readsource,
    cartesianindex,
    costfunction_gridded_obs, costfunction_gridded_obs!,
    costfunction_point_obs, costfunction_point_obs!,
    costfunction_gridded_model, costfunction_gridded_model!,
    trackpathways, meanage,
    regeneratedphosphate, preformedphosphate,
    regeneratednitrate, preformednitrate,
    respiredoxygen, preformedoxygen,
    volumefilled, surfaceorigin, synthetic_observations,
    observe,
    steadyclimatology, steadyinversion,
    interpweights, interpindex,
    wetlocation, iswet,
    control2state, control2state!,
    surfacecontrol2field, surfacecontrol2field!,
    sparsedatamap, config2nc, axislabels,
    matrix_zyx2xyz,
    location_obs,
    zerosurfaceboundary,
    onesurfaceboundary,
    getsurfaceboundary,getnorthboundary,geteastboundary,
    getsouthboundary,getwestboundary,
    setboundarycondition!,
    wetmask, interiormask, mixedlayermask, boundarymask,
    adjustboundarycondition, adjustboundarycondition!,
    gsetboundarycondition, setsource!,
    zeros, one, oneunit, ones,
    #maximum, minimum,
    +, -, *, dot,
    zerosource, onesource,
    adjustsource, adjustsource!,
    Grid, Field, BoundaryCondition, vec, unvec!, unvec, wet,
    zerowestboundary, zeronorthboundary,
    zeroeastboundary, zerosouthboundary,
    onewestboundary, onenorthboundary, oneeastboundary, onesouthboundary,
    distancematrix, gaussiandistancematrix, versionlist,
    massfractions, massfractions_isotropic, neighbors, 
    Observations, SourceControls, BoundaryControls, 
    MassFracControls, Controls

import Base: zeros, one, oneunit, ones, \
import Base: maximum, minimum
import Base: +, -, *, /, vec
import LinearAlgebra: dot

# Credit to DrWatson.jl for these functions
# Didn't want to add dependency for these small functions
#projectdir() = dirname(Base.active_project())

# find packagedir even if TMI is not the active project
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

pkgsrcdir() = joinpath(pkgdir(),"src")
pkgsrcdir(args...) = joinpath(pkgsrcdir(), args...)

pkgplotsdir() = joinpath(pkgdir(),"plots")
pkgplotsdir(args...) = joinpath(pkgplotsdir(), args...)

pkgutilsdir() = joinpath(pkgdir(),"utils")
pkgutilsdir(args...) = joinpath(pkgutilsdir(), args...)

include(pkgsrcdir("grid.jl"))
include(pkgsrcdir("field.jl"))
include(pkgsrcdir("source.jl"))
include(pkgsrcdir("config.jl"))
include(pkgsrcdir("boundary_condition.jl"))
include(pkgsrcdir("regions.jl"))
include(pkgsrcdir("mass_fractions.jl"))
include(pkgsrcdir("deprecated.jl"))

include(pkgsrcdir("observations.jl"))
include(pkgsrcdir("controls", "utils.jl"))
include(pkgsrcdir("controls", "boundary.jl"))
include(pkgsrcdir("controls", "source.jl"))
include(pkgsrcdir("controls", "massfrac.jl"))
include(pkgsrcdir("controls", "control_parameters.jl"))

""" 
    function trackpathways(TMIversion,latbox,lonbox)
    Track the pathways of a user-defined water mass.
     Steps: (a) define the water mass by a rectangular surface patch dyed with passive tracer concentration of         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.
     See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `latbox`: min and max latitude of box
- `lonbox`: min and max longitude of box
- `γ`: TMI grid
# Output
- `c`: fraction of water from surface source
"""
function trackpathways(Alu,latbox,lonbox,γ)

    b = surfacepatch(lonbox,latbox,γ)
    c = steadyinversion(Alu,b,γ)

    return c
end

""" 
    function watermassdistribution(TMIversion,latbox,lonbox)
    Track the pathways of a user-defined water mass.
     Steps: (a) define the water mass by an oceanographically-relevant surface patch dyed with passive tracer concentration of one
         (b) propagate the dye with the matrix A, with the result being the fraction of water originating from the surface region.
     See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `region`: name of pre-defined surface region
- `γ`: TMI grid
# Output
- `g`: water-mass fraction
"""
function watermassdistribution(TMIversion,Alu,region,γ)
    
    b = surfaceregion(TMIversion,region) # version 2 of this routine
    g = steadyinversion(Alu,b,γ)
    return g
end

""" 
    function regeneratednutrient(TMIversion,Alu,γ)

    Regenerated (i.e., accumulated, remineralized) nutrient

# Arguments
- `tracer::Union{String,Symbol}`: tracer name
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
- `r`: optional stoichiometric ratio relative to PO₄
# Output
- `PO₄ᴿ`: regenerated phosphate
"""
function regeneratednutrient(tracer,
    TMIversion,Alu,γ;r=1)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # configure meta data
    c = readfield(TMIfile,tracer,γ) 
    name = Symbol(c.name,"ᴿ")
    if c.longname == "dissolved oxygen"
        longname = "respired "*c.longname
    else
        longname = "regenerated "*c.longname
    end    
    
    ## read phosphate source
    qPO₄ = readsource(TMIfile,"qPO₄",γ)

    # zero boundary condition
    b = zerosurfaceboundary(γ,name,longname,c.units)
    cᴿ = steadyinversion(Alu,b,γ,q=qPO₄,r=r)
    return cᴿ
end
regeneratedphosphate(TMIversion,Alu,γ) = regeneratednutrient("PO₄",TMIversion,Alu,γ,r=1)
regeneratednitrate(TMIversion,Alu,γ) = regeneratednutrient("NO₃",TMIversion,Alu,γ,r=15.5)
respiredoxygen(TMIversion,Alu,γ) = regeneratednutrient("O₂",TMIversion,Alu,γ,r=-170.0)

""" 
    function preformednutrient(tracer::Union{String,Symbol},TMIversion,Alu,γ)

    Preformed (i.e., NO accumulation or remineralization) nutrient
# Arguments
- `tracer::Union{Symbol,String}`: tracer name
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `c★`: preformed tracer
"""
function preformednutrient(tracer::Union{String,Symbol},TMIversion,Alu,γ)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # get meta-data
    c = readfield(TMIfile,tracer,γ)
    name = Symbol(tracer,"★")
    longname = "preformed "*c.longname

    b = getsurfaceboundary(Field(c.tracer,c.γ,name,longname,c.units))
    return steadyinversion(Alu,b,γ) 
end

"""
preformedphosphate(TMIversion,Alu,γ) = preformednutrient("PO₄",TMIversion,Alu,γ)
"""
preformedphosphate(TMIversion,Alu,γ) = preformednutrient("PO₄",TMIversion,Alu,γ)

"""
preformednitrate(TMIversion,Alu,γ) = preformednutrient("NO₃",TMIversion,Alu,γ)
"""
preformednitrate(TMIversion,Alu,γ) = preformednutrient("NO₃",TMIversion,Alu,γ)

"""
preformedoxygen(TMIversion,Alu,γ) = preformednutrient("O₂",TMIversion,Alu,γ)
"""
preformedoxygen(TMIversion,Alu,γ) = preformednutrient("O₂",TMIversion,Alu,γ)

"""
preformedcarbon13(TMIversion,Alu,γ) = preformednutrient("δ¹³C",TMIversion,Alu,γ)
"""
preformedcarbon13(TMIversion,Alu,γ) = preformednutrient("δ¹³C",TMIversion,Alu,γ)

""" 
    function meanage(TMIversion,Alu,γ)
    Mean or ideal age
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `a`: mean age [yr]
"""
function meanage(TMIversion,Alu,γ)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    if TMIfile[end-1:end] == "nc"

        #F = ncread(file,"F")
        ## read age source
        F₀ = readfield(TMIfile,"F₀",γ)
        qPO₄ = readsource(TMIfile,"qPO₄",γ) # use this to define mixed-layer

        # better to define a Source type
        Iq = findall(x -> x > 0,qPO₄.tracer)

        qa = zeros(γ)
        qa.tracer[Iq] = 1 ./ F₀.tracer[Iq]
        # zero boundary condition
        b₀ = zerosurfaceboundary(γ)
        a = steadyinversion(Alu,b₀,γ,q=qa)

    else
        
        error("not implemented for mat input file")
    end
        
    return a
end

function local_residence_time(TMIfile, A, γ)

    ## read age source
    F₀ = readfield(TMIfile,"F₀",γ)

    # qa = age source
    qa = zeros(γ)

    # no source in mixed layer
    mixedlayer = mixedlayermask( A, γ)

    # no source on boundary
    boundary = (γ.wet .&&  .!γ.interior)

#    qPO₄ = readsource(TMIfile,"qPO₄",γ) # use this to define mixed-layer
    # better to define a Source type
 #   Iq = findall(x -> x > 0,qPO₄.tracer)

    for i in eachindex(F₀.tracer)
        if (!mixedlayer[i] && !boundary[i]) 
            qa.tracer[i] = 1 / F₀.tracer[i]
        end
    end
    return qa
end

function meanage(TMIfile, A, b, γ)

    # ## read age source
    # F₀ = readfield(TMIfile,"F₀",γ)
    # qPO₄ = readsource(TMIfile,"qPO₄",γ) # use this to define mixed-layer

    # # better to define a Source type
    # Iq = findall(x -> x > 0,qPO₄.tracer)
    
    # # qa = age source
    # qa = zeros(γ)
    # qa.tracer[Iq] = 1 ./ F₀.tracer[Iq]

    # bmask = zeros(γ)
    # setboundarycondition!(bmask, b1)
    # # setboundarycondition!(bmask, b_south)
    # # setboundarycondition!(bmask, b_up)
    # # setboundarycondition!(bmask, b_lo)

    # vbmask = vec(bmask)
    # sumA = sum(A, dims=2)
    # for i in eachindex(vbmask)
    #     if vbmask[i] > 0.0
    #         println("i=",i)
    #         println("sum(A)=",sumA[i])
    #         (sumA[i] != 1.0) && error(" row should sum to one ")
    #         qa.tracer[i] = 0.0
    #         #qa.tracer[γ.wet][i] = 0.0
    #     end
    # end

    # for i in eachindex(bmask)
    #     if bmask[i] > 0.0
    #         println("i=",i)
    #         println("sum(A)=",sumA[i])
    #         (sumA[i] != 1.0) && error(" row should sum to one ")
    #         qa.tracer[i] = 0.0
    #     end
    # end

    qa = local_residence_time(TMIfile, A, γ)

    println("max age source ", maximum(qa))
    println("min age source ", minimum(qa))
    # zero boundary condition
    a = steadyinversion(lu(A), b, γ, q=qa)
#    a = steadyinversion(A,b0,γ,q=qa)

    return a#, qa
end

""" 
    function volumefilled(TMIversion)
    Find the ocean volume that has originated from each surface box.
     This is equivalent to solving a sensitivity problem:
     The total volume is V = vᵀ c , where v is the volume of each box 
     and c is the fraction of volume from a given source which
     satisfies the equation A c = d.                     
     Next, dV/d(d) = A⁻ᵀ v, and dV/d(d) is exactly the volume originating from each source.

     See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011. 
# Arguments
- `TMIversion`: version of TMI water-mass/circulation model
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI.grid
# Output
- `volume`: log10 of global ocean volume filled by a surface region, exists at surface, therefore given BoundaryCondition type
"""
function volumefilled(TMIversion,Alu,γ)::BoundaryCondition

    v = cellvolume(γ)
    area = cellarea(γ)

    vfilled = 0*area # answer will go here

    # effectively take inverse of transpose A matrix.
    #dVdd = zeros(γ) # pre-allocate array
    dVdd = Alu'\v #[γ.wet]
    dimval = surfaceindex(γ)
    # scale the sensitivity value by surface area so that converging meridians are taken into account.
    I = γ.I
    #volume = zeros(Float64,length(γ.lon),length(γ.lat))
    volume = zeros(γ.wet[:,:,dimval])
    # this step could use a function with γ.I argument

    for ii ∈ eachindex(I)
        if I[ii][3] == dimval
            volume[I[ii][1],I[ii][2]] = dVdd.tracer[I[ii][1],I[ii][2],dimval]./area.tracer[I[ii][1],I[ii][2]]
        end
    end
             
    volume = log10.(volume)
    ∂V∂b  = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[dimval],3,dimval,γ.wet[:,:,dimval],:V,"volume filled by surface gridcell","log₁₀(m³/m²)")
    
    return  ∂V∂b 
end

"""
     function effective_endmember(TMIversion,Alu,field,region,γ)

Effective (i.e., importance-weighted) endmember tracer value
calculated according to Gebbie and Huybers 2011.
"""
function effective_endmember(TMIversion,Alu,field::Field,region,γ::Grid)
    b = surfaceregion(TMIversion,region) # version 2 of this routine
    return effective_endmember(Alu,field,b,γ)
end
function effective_endmember(Alu,field::Field,b::BoundaryCondition,γ::Grid)
    volweight, bc = effective_endmember_sums(Alu,field,b,γ)
    
    # caution: did not expicitly check that the wet masks matched up
    return sum(volweight.tracer[wet(volweight)].*
               bc.tracer[wet(bc)].*
               b.tracer[wet(b)])/
           sum(volweight.tracer[wet(volweight)].*
               b.tracer[wet(b)])
end
function effective_endmember(Alu,field::Field,b::NamedTuple,γ::Grid; control_volume=γ.wet)

    numer = 0.0
    denom = 0.0 
    for bi in b 
        volweight, bc = effective_endmember_sums(Alu,field,bi,γ; control_volume=control_volume)
        numer += sum(volweight.tracer[wet(volweight)].*
                     bc.tracer[wet(bc)].*
                     bi.tracer[wet(bi)])

        denom += sum(volweight.tracer[wet(volweight)].*
                     bi.tracer[wet(bi)])
    end
    
    # caution: did not expicitly check that the wet masks matched up
    return numer/denom 
end

"""
     function effective_endmember_sums(Alu,field::Field,b::BoundaryCondition,γ::Grid)

Intermediate quantities for computing effective endmembers
"""
function effective_endmember_sums(Alu,field::Field,b::BoundaryCondition,γ::Grid; control_volume = γ.wet)
    v = cellvolume(γ)

    for i in eachindex(control_volume)
        if !control_volume[i] 
            v.tracer[i] = 0.0
        end
    end
    
    # effectively take inverse of transpose A matrix.
    dVdd = Alu'\v 

    # pick out the relevant part along the boundary
    volweight = getboundarycondition(dVdd,b.dim,b.dimval) 
    bc = getboundarycondition(field,b.dim,b.dimval)
    return volweight, bc
end

""" 
    function surfaceorigin(TMIversion,loc)
     Find the surface origin of water for some interior box 
     This is equivalent to solving a sensitivity problem:
     The mass fraction at a location `loc` of interest is 
    `c[loc] = δᵀ c`, where `δ` samples the location of the global mass-fraction variable, c.
    Then the sensitivity of `c[loc]` is: d(c[loc])/d(d) = A⁻ᵀ δ.
    The derivative is solved using the constraint: Ac = d.
    The sensitivity is exactly the mass fraction originating from each source.      
    This problem is mathematically similar to determining how the ocean is filled.
# Arguments
- `loc`: location (lon,lat,depth) of location of interest
- `Alu`: LU decomposition of water-mass matrix A
- `γ`: TMI grid
# Output
- `origin`: surface map of fraction of source water for a given location, log10 of effective depth, in terms of a BoundaryCondition
"""
function surfaceorigin(loc,Alu,γ::Grid)::BoundaryCondition

    δ = interpweights(loc,γ)
    dvlocdd = zeros(γ.wet); # pre-allocate c
    dvlocdd[γ.wet] = Alu'\δ[γ.wet]
    dimval = surfaceindex(γ)
    # origin is defined at sea surface
    #origin = view(dvlocdd,:,:,1)

    println(sum(dvlocdd[:,:,dimval][γ.wet[:,:,dimval]]))
    dvlocdd = log10.(dvlocdd[:,:,dimval])
    small_cutoff = -10 # 1e-10
    replace!(x -> x < small_cutoff ? small_cutoff : x,dvlocdd) 
    origin = BoundaryCondition(dvlocdd,(γ.lon,γ.lat),γ.depth[dimval],3,dimval,γ.wet[:,:,dimval],:origin,"surface origin","log₁₀(m³)")
    
    return origin
end

"""
function steadyclimatology_optim(u₀,fg!,iterations)
     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u₀
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)
    subject to Aỹ = d₀ + Γ u₀.                 
    W⁻ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u₀`:
- `fg!`: compute cost function and gradient in place
- `iterations`: number of optimization iterations
"""
function steadyclimatology_optim(u₀,fg!,iterations)
#function steadyclimatology(u₀,Alu,d₀,y,W⁻,fg!,γ)

    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    #out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(),Optim.Options(show_trace=true, iterations = iterations))

    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(linesearch = LineSearches.BackTracking()),Optim.Options(show_trace=true, iterations = iterations))

    return out    
end

function steadyclimatology(Alu,b,u,y,W⁻,γ)
    uvec = vec(u)
    F,G = costfunction_gridded_obs(uvec,Alu,b,u,y,W⁻,γ)
    fg!(F,G,x) = costfunction_gridded_obs!(F,G,x,Alu,b,u,y,W⁻,γ)
    fg(x) = costfunction_gridded_obs(x,Alu,b,u,y,W⁻,γ)
    f(x) = fg(x)[1]

    J₀,gJ₀ = fg(uvec)
    iterations = 10
    out = steadyclimatology_optim(uvec,fg!,iterations)
    return out, f, fg, fg!
end

"""
    function sparsedatamap(u₀::Vector{T},Alu,b::BoundaryCondition{T},u::BoundaryCondition{T},y::Vector{T},W⁻,wis,locs,Q⁻,γ::Grid;iterations=10) where T <: Real

     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u₀
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - y)ᵀ W⁻¹ (ỹ - y)
    subject to Aỹ = d₀ + Γ u₀.                 
    W⁻ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u₀`:
- `Alu`:
- `b`: first guess of boundary conditions and interior sources
- `y`: observations on 3D grid
- `W⁻`: weighting matrix best chosen as inverse error covariance matrix
- `fg!`: compute cost function and gradient in place
- `γ`: grid
"""
function sparsedatamap_optim(u₀::Vector,Alu,b::Union{BoundaryCondition,NamedTuple},u::Union{BoundaryCondition,NamedTuple},y::Vector,W⁻,wis::Vector,locs,Q⁻,γ::Grid;q = nothing, r = 1.0,iterations=10) #where {N1, N2, T <: Real}
#function sparsedatamap(u₀::Vector{T},Alu,b::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Vector{T},W⁻,wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}}},locs,Q⁻,γ::Grid,iterations=10) where {N1, N2, T <: Real}

     fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q₀=q,r=r)
    
    # a first guess: observed surface boundary conditions are perfect.
    # set surface boundary condition to the observations.
    out = optimize(Optim.only_fg!(fg!), u₀, LBFGS(linesearch = LineSearches.BackTracking()),Optim.Options(show_trace=true, iterations = iterations))

    return out    
end

function sparsedatamap(Alu,b,u,y,W⁻,wis,locs,Q⁻,γ;q = nothing, r = 1.0,iterations=10)
    fg(x) = costfunction_point_obs(x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=q,r=r)
    f(x) = fg(x)[1]
    #J0 = f(uvec)
    #J₀,∂J₀∂u = fg(uvec)
    uvec = vec(u)
    fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q₀=q,r=r)
    out = sparsedatamap_optim(uvec,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=q,r=r,iterations=iterations)
    return out, f, fg, fg!
end

function gradient_check(uvec,f,fg,fg!)
    # check with forward differences
    ϵ = 1e-3
    ii = rand(1:length(uvec))
    println("Location for test =",ii)
    δu = copy(uvec); δu[ii] += ϵ
    ∇f_finite = (f(δu) - f(uvec))/ϵ

    J₀,∂J₀∂u = fg(uvec)
    fg!(J₀,∂J₀∂u,(uvec+δu)./2) # J̃₀ is not overwritten
    ∇f = ∂J₀∂u[ii]
    println("∇f=",∇f)

    # error less than 10 percent?
    println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
    return ∇f, ∇f_finite
end


"""
    function writefield(file,field)

    Write a Field to NetCDF.
 
    Use NCDatasets so that Unicode is correct

# Arguments
- `file`: TMI NetCDF file name
- `field::Field`: a TMI.Field struct
# Output
- none
# Side-effect
- write to `file`
"""
function write(file,field::Union{Source{T},Field{T}}) where T <: Real

    if !isfile(file)
        # create new NetCDF file
        ds = Dataset(file,"c")

        TMIgrids, TMIgridsatts = griddicts(field.γ)

        # Define the dimension "lon" and "lat" with the size 100 and 110 resp.
        defDim(ds,"lon",length(field.γ.lon))
        defDim(ds,"lat",length(field.γ.lat))
        defDim(ds,"depth",length(field.γ.depth))

        # Define a global attribute
        ds.attrib["title"] = "TMI output"

        vlon = defVar(ds,"lon",Float64,["lon"],
               attrib = OrderedDict(TMIgridsatts["lon"]))
        vlon[:] = field.γ.lon

        vlat = defVar(ds,"lat",Float64,["lat"],
               attrib = OrderedDict(TMIgridsatts["lat"]))
        vlat[:] = field.γ.lat

        vdepth = defVar(ds,"depth",Float64,["depth"],
               attrib = OrderedDict(TMIgridsatts["depth"]))
        vdepth[:] = field.γ.depth
    

        v = defVar(ds,String(field.name),Float64,("lon","lat","depth"),
                  attrib = OrderedDict("longname" => field.longname,
                                "units" => field.units))
        v[:,:,:] = field.tracer

        close(ds)

    else
        # assumption: on the same grid
        ds = Dataset(file,"a")

        println(field.name)
        v = defVar(ds,String(field.name),Float64,("lon","lat","depth"),
                  attrib = OrderedDict("longname" => field.longname,
                                "units" => field.units))
        v[:,:,:] = field.tracer
        close(ds)
    end
        
    return nothing
end

"""
    function depthindex(I) 
    Get the k-index (depth level) from the Cartesian index
"""
function depthindex(I)
    T = eltype(I[1])
    k = Vector{T}(undef,length(I))
    [k[n]=I[n][3] for n ∈ eachindex(I)]
    return k
end

"""
    function lonindex(I) 
    Get the i-index (lon index) from the Cartesian index
"""
function lonindex(I)
    T = eltype(I[1])
    i = Vector{T}(undef,length(I))
    [i[n]=I[n][1] for n ∈ eachindex(I)]
    return i
end

"""
    function latindex(I) 
    Get the j-index (latitude index) from the Cartesian index
"""
function latindex(I)
    T = eltype(I[1])
    j = Vector{T}(undef,length(I))
    [j[n]=I[n][2] for n ∈ eachindex(I)]
    return j
end

function findindex(I,indexfunc,indexnumber)
    Ifound = findall(indexfunc(I) .== indexnumber)
    return Ifound
end
    
"""
    function surfaceindex(I) 
    Get the vector-index where depth level == 1 and it is ocean.
"""
surfaceindex(I) = findindex(I,depthindex,1)

"""
    function southindex(I) 
    Get the vector-index on the southern open boundary
"""
southindex(I) = findindex(I,latindex,1)

"""
    function northindex(I) 
    Get the vector index on the northern open boundary
"""
northindex(I) = findindex(I,latindex,maximum(latindex(I)))

"""
    function westindex(I) 
    Get the vector index on the western open boundary
"""
westindex(I) = findindex(I,lonindex,1)

"""
    function eastindex(I) 
    Get the vector index on the northern open boundary
"""
eastindex(I) = findindex(I,lonindex,maximum(lonindex(I)))

"""
    Horizontal area of grid cell
"""
function cellarea(γ)
    dx = zonalgriddist(γ)
    dy = haversine((γ.lon[1],γ.lat[1])
                  ,(γ.lon[1],γ.lat[2]))

    area = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
    fill!(area,0.0)

    # to calculate area everywhere
    #[area[i,j] = dx[j] * dy for i ∈ eachindex(γ.lon) for j ∈ eachindex(γ.lat)]

    dimval = surfaceindex(γ) # surface

    # to calculate sea surface area
    I = γ.I
    [area[I[ii][1],I[ii][2]] = dx[I[ii][2]] * dy for ii ∈ eachindex(I) if I[ii][3] == dimval]

    dim = 3 # 3rd dimension is fixed 

    # is it really a Boundary Condition? (sorta, but more of a 2D Field)
    return BoundaryCondition(area,(γ.lon,γ.lat),γ.depth[dimval],dim,dimval,γ.wet[:,:,dimval],
             :area,"cell area","m²")
    #return area
end

"""
    function distancematrix(γ;surface = true)

    Matrix with size of surface points squared

    Each entry gives distance in km between surface points
    Gives only horizontal distance.
"""
function distancematrix(γ;surface = true)
    if surface
        dimval = surfaceindex(γ)
        nsfc = sum(γ.wet[:,:,dimval])
        Dh = zeros(nsfc,nsfc)
        lonlist = [γ.lon[γ.I[ii][1]] for ii in 1:nsfc]
        latlist = [γ.lat[γ.I[ii][2]] for ii in 1:nsfc]

        for ii in 1:nsfc
            loc1 = (lonlist[ii],latlist[ii])
            for jj in ii:nsfc
                loc2 = (lonlist[jj],latlist[jj])
                Dh[ii,jj] = TMI.haversine(loc1,loc2)./1000.0 # m -> km
                Dh[jj,ii] = Dh[ii,jj] # make symmetric, save time
            end
        end
    else
        error("not implemented for whole ocean")
    end
    return Dh
end

"""
    function gaussiandistancematrix(γ,σ,L)

    uses distance matrix plus a lengthscale `L` (km)
    and a size of the diagonal `σ`

    returns values with inverse gaussian weighting
"""
function gaussiandistancematrix(γ,σ,L)

    adhoc_factor = 0.01 # to make non-negative matrix
    Dh = distancematrix(γ,surface=true)
    factor1 = (1-adhoc_factor)*σ^2
    factor2 = adhoc_factor*σ^2
    N= size(Dh,1)
    Dg = factor1 .* exp.(-(Dh./L).^2) + factor2*Diagonal(ones(N))
end

"""
    function cellvolume(γ)::Field

    Volume of each grid cell.
"""
function cellvolume(γ)
    dz = layerthickness(γ)
    area = cellarea(γ)
    volume = Array{Float64,3}(undef,length(γ.lon),length(γ.lat),length(γ.depth))
    fill!(volume,0.0)

    # for volume everywhere
    # [volume[i,j,k] = area[i,j] * dz[k] for i ∈ eachindex(γ.lon) for j ∈ eachindex(γ.lat) for k ∈ eachindex(γ.depth)]

    # for ocean volume only
    I = γ.I
    [volume[I[ii]] = area.tracer[I[ii][1],I[ii][2]] * dz[I[ii][3]] for ii ∈ eachindex(I)]

    # turn it into a Field
    return Field(volume,γ,:vol,"cell volume","m³")
end

function layerthickness(γ::Grid)
    zface= (γ.depth[1:end-1].+γ.depth[2:end])./2;
    dz = ([zface[1] ; diff(zface); 500]);
    return dz
end

function zonalgriddist(γ::Grid)
    dx = similar(γ.lat)
    for j in eachindex(γ.lat)
        dx[j] = haversine((γ.lon[1],γ.lat[j])
                         ,(γ.lon[2],γ.lat[j]))
    end
    return dx
end

"""
    function nearestneighbormask
    Make a 3D tracer field that is 1 at location 
    of nearest neighbor, 0 elsewhere
# Arguments
- `loc`: location in a 3-tuple (lon,lat,depth)
- `γ`: TMI.grid
# Output
- `δ`: nearest neighbor mask 3D field
"""
function nearestneighbormask(loc,γ::Grid,N=1)

    Inn, Rnn = nearestneighbor(loc,γ,N)

    # preallocate
    δ = falses(size(γ.wet))
    #Array{BitArray,3}(undef,size(γ.wet))
    #fill!(δ,zero(Bool))
    δ[Inn] = 1 
    return δ
end

"""
    function nearestneighbor(loc,γ)
    return the Cartesian index and linear index 
    of the nearest N neighbors
# Arguments
- `loc`: 3-tuple of lon,lat,depth location
- `γ`: TMI.grid
# Output
- `Inn`: Cartesian indices of nearest neighbor
#- `Rnn`: linear indices of nearest neighbor, Removed from code
"""
function nearestneighbor(loc,γ,N=1)

    xydist = horizontaldistance(loc[1:2],γ)

    if N==1
        ijdist,ijmin = findmin(xydist[γ.wet[:,:,1]])
        kdist,kmin = findmin(abs.(loc[3] .- γ.depth))
    elseif N > 1
        ijmin = sortperm(xydist[γ.wet[:,:,1]])
        kmin = sortperm(abs.(loc[3] .- γ.depth))
    end

    if N == 1
        Inn = CartesianIndex.(γ.I[ijmin][1],γ.I[ijmin][2],kmin)
    elseif N > 1
        Inn = Vector{CartesianIndex}(undef,N)
        cN2 = ceil(Integer,N/2)
        for ii in 1:cN2
            # translate ijmin into imin, jmin
            Inn[ii] = CartesianIndex.(γ.I[ijmin[ii]][1],γ.I[ijmin[ii]][2],kmin[1])
        end
        for ii in 1:floor(Integer,N/2)
            Inn[cN2+ii] = CartesianIndex.(γ.I[ijmin[ii]][1],γ.I[ijmin[ii]][2],kmin[2])
        end        
    end
    
    return Inn
end

"""
    function horizontaldistance(loc,γ)

# Arguments
- `loc`: 3-tuple of lon,lat,depth location
- `γ`: TMI.grid
# Output
- `hordist`: horizontal distance to nearest tracer grid points
"""
function horizontaldistance(loc,γ::Grid)

    # hordist will have same type as lon,lat,depth
    T = eltype(γ.lon)
    
    # pre-allocate horizontal distance
    hordist = Matrix{T}(undef,length(γ.lon),length(γ.lat))
    # will give NaN with appropriate precision
    fill!(hordist,zero(T)/zero(T))
    
    # calculate haversine horizontal distance on sphere
    [hordist[γ.I[ii]] = haversine((loc[1],loc[2]),(γ.lon[γ.I[ii][1]],γ.lat[γ.I[ii][2]]))
       for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]
    return hordist
end


"""
    function interpindex(loc::T, γ::Grid{R,1}) where {T,R}

    1D grid specialization of `interpindex`.
    See the 3D version for more details.
# Arguments
- `loc`: location
- `γ`: TMI grid
# Output
- `wis`: weighted indexes for interpolation
"""
function interpindex(loc::T, γ::Grid{R,1}) where {T,R}
    loc_on_grid = shiftloc(loc, γ)
    nodes = (γ.axes[1],)
    return Interpolations.weightedindexes(
        (Interpolations.value_weights,),
        (Gridded(Linear()),),
        nodes,
        (loc_on_grid,),
    )
end

"""
    function interpindex(loc::AbstractVector{T}, γ::Grid{R,1}) where {T,R}

    1D grid specialization of `interpindex` for a vector location.
    See the 3D version for more details.
# Arguments
- `loc`: a 1-element vector location
- `γ`: TMI grid
# Output
- `wis`: weighted indexes for interpolation
"""
function interpindex(loc::AbstractVector{T}, γ::Grid{R,1}) where {T,R}
    if length(loc) != 1
        throw(ArgumentError("interpindex 1D expects a length-1 vector, got length $(length(loc))"))
    end
    return interpindex(loc[1], γ)
end

"""
function interpindex(loc,γ)
    Weights for linear interpolation.
    The derivative of linear interpolation is needed in sensitivity studies.
    ReverseDiff.jl could find this quantity automatically.
    Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.
# Arguments
- `c`: a temporary tracer field, would be nice to make it unnecessary
- `loc`: (lon,lat,depth) tuple of a location of interest
- `γ`: TMI grid
# Output
- `δ`: weights on a 3D tracer field grid
"""
function interpindex(loc::Tuple{T, T, T},γ) where T

    # Handle a grid mismatch.
    loc_on_grid = shiftloc(loc,γ)
    
    # Handle longitudinal periodic condition (i.e., wraparound)
    lon = vcat(copy(γ.lon),γ.lon[1]+360.)
    list = vcat(1:length(γ.lon),1)
    nodes = (lon,γ.lat,γ.depth)

    # eliminate need to pass tracer value
    wis = Interpolations.weightedindexes((Interpolations.value_weights,),((Gridded(Linear()), Gridded(Linear()), Gridded(Linear()))),nodes,loc_on_grid)

    # issue, some of weighted points may be NaNs in tracer field
    # handle this in the Interpolations.jl routines
    # may involve chaging Gridded(Linear()) above
    return wis
end

"""
    function shiftloc(loc::T, γ::Grid{R,1}) where {T,R}
    1D grid specialization of `shiftloc`.
    Check if location is on the 1D grid.
# Arguments
- `loc`: location
- `γ`: TMI grid
# Output
- `loc`: original location if on grid
"""
function shiftloc(loc::T, γ::Grid{R,1}) where {T,R}
    x = γ.axes[1]

    if length(x) == 1
        if !isapprox(loc, x[1]; atol = zero(T))
            error("location not on grid")
        end
        return loc
    end

    xmin, xmax = extrema(x)
    if loc < xmin || loc > xmax
        error("location not on grid")
    end

    return loc
end

"""
    function shiftloc(loc::AbstractVector{T}, γ::Grid{R,1}) where {T,R}
    1D grid specialization of `shiftloc` for a vector location.
# Arguments
- `loc`: a 1-element vector location
- `γ`: TMI grid
# Output
- `loc`: original location if on grid
"""
function shiftloc(loc::AbstractVector{T}, γ::Grid{R,1}) where {T,R}
    if length(loc) != 1
        throw(ArgumentError("shiftloc 1D expects a length-1 vector, got length $(length(loc))"))
    end
    return shiftloc(loc[1], γ)
end

"""
    shiftloc(loc::Tuple, γ)

Wrap longitudes into the grid domain (allowing a half-cell overlap) while
leaving latitude/depth untouched. Errors if the adjusted longitude is still
off-grid.
"""
function shiftloc(loc::Tuple{T, T, T},γ) where T
    # accounts for a half grid cell of overlap space
    newlon = loc[1]
    westlon = (3/2)*γ.lon[1] - (1/2)*γ.lon[2]
    eastlon =   (3/2)*γ.lon[end] - (1/2)*γ.lon[end-1]
    while newlon <  westlon
        newlon += 360.0
    end
    
    while newlon > eastlon
        newlon -= 360.0
    end

    if newlon > eastlon || newlon < westlon
        error("location not on grid")
    end

    return (newlon, loc[2], loc[3])    
end

"""
    function interpweights(loc::T, γ::Grid{R,1}; wis = nothing) where {T,R}
    1D grid specialization of `interpweights`.
# Arguments
- `loc`: location
- `γ`: TMI grid
- `wis`: pre-computed weighted indexes
# Output
- `δ`: weights on a 1D tracer field grid
"""
function interpweights(loc::T, γ::Grid{R,1}; wis = nothing) where {T,R}
    if isnothing(wis)
        wis = interpindex(loc, γ)
    end

    δ = zeros(γ.wet)
    for ii in 1:length(wis[1].weights)
        δ[wis[1].istart + ii - 1] += wis[1].weights[ii]
    end

    if iszero(sum(filter(!isnan,δ)))
        δ = nothing
    elseif sum(filter(!isnan,δ)) < 1.0
        δ ./= sum(filter(!isnan,δ))
    end
    
    return δ
end

"""
    function interpweights(loc::AbstractVector{T}, γ::Grid{R,1}; wis = nothing) where {T,R}
    1D grid specialization of `interpweights` for a vector location.
# Arguments
- `loc`: 1-element vector location
- `γ`: TMI grid
- `wis`: pre-computed weighted indexes
# Output
- `δ`: weights on a 1D tracer field grid
"""
function interpweights(loc::AbstractVector{T}, γ::Grid{R,1}; wis = nothing) where {T,R}
    if length(loc) != 1
        throw(ArgumentError("interpweights 1D expects a length-1 vector, got length $(length(loc))"))
    end
    return interpweights(loc[1], γ; wis)
end

"""
function interpweights(loc,γ)
    Weights for linear interpolation.
    The derivative of linear interpolation is needed in sensitivity studies.
    ReverseDiff.jl could find this quantity automatically.
    Instead we dig into the Interpolations.jl package to find the weights that are effectively the partial derivatives of the function.
# Arguments
- `loc`: (lon,lat,depth) tuple of a location of interest
- `γ`: TMI grid
# Output
- `δ`: weights on a 3D tracer field grid
"""
function interpweights(loc::Tuple{T, T, T},γ; wis = nothing) where T

    # handle wraparound
    # repeated (unnecessarily?) in interpindex
    lon = vcat(copy(γ.lon),γ.lon[1]+360)
    list = vcat(1:length(γ.lon),1)
    if isnothing(wis)
        wis = interpindex(loc,γ)
    end

    # translate to weights via
    #http://juliamath.github.io/Interpolations.jl/latest/devdocs/
    δ = zeros(γ.wet)

    # changes in δwrap i=91 are translated back to δ i=1
    δwrap = view(δ,list,:,:)
    for ii = 1:2
        for jj = 1:2
            for kk = 1:2
                δwrap[wis[1].istart+ii-1,wis[2].istart+jj-1,wis[3].istart+kk-1] +=
                wis[1].weights[ii]*wis[2].weights[jj]*wis[3].weights[kk]
            end
        end
    end

    # if some adjacent points are dry, then re-normalize to keep this interpolation as an average.
    # The hope is that the interpolation is stable with this approach, but other side effects are likely.
    # Note that this should be handled earlier, like in the interpindex section. For this reason, there could be an inconsistency in the global map function.
    if iszero(sum(filter(!isnan,δ)))
        δ = nothing
    elseif sum(filter(!isnan,δ)) < 1.0
        δ ./= sum(filter(!isnan,δ))
    end
    
    return δ
end

function zerosource(γ::Grid,name=:none,longname="unknown",units="unknown";logscale=false)::Source
    T = eltype(γ.depth)
    tracer = Array{T}(undef,size(γ.interior))
    tracer[γ.interior] .= zero(T)
    tracer[.!γ.interior] .= zero(T)/zero(T)
    return Source(tracer,γ,name,longname,units,logscale)
end

function onesource(γ::Grid,name=:none,longname="unknown",units="unknown";logscale=false)::Source
    T = eltype(γ.depth)
    tracer = Array{T}(undef,size(γ.interior))
    tracer[γ.interior] .= one(T)
    tracer[.!γ.interior] .= zero(T)/zero(T)
    return Source(tracer,γ,name,longname,units,logscale)
end

""" 
    function ones(γ::Grid,name=:none,longname="unknown",units="unknown")::Field

      initialize tracer field of ones on TMI grid
      using a Field struct and constructor
# Arguments
- `γ`::TMI.Grid
# Output
- `d`::Field,  3d tracer field with NaN on dry points
"""
function ones(γ::Grid,name=:none,longname="unknown",units="unknown")::Field

    # use depth (could have been lon, lat)
    # to get element type
    T = eltype(γ.depth)
    
    # preallocate
    tracer = Array{T}(undef,size(γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    tracer[γ.wet] .= Base.one(T) # add Base: error "should import Base"
    tracer[.!γ.wet] .= zero(T)/zero(T) # NaNs with right type
    d = Field(tracer,γ,name,longname,units)

    return d
end

""" 
   function oneunit, help for gridded Interpolations
"""
function one(field::Field{T})::Field{T} where T <: Real

    # use depth (could have been lon, lat)
    # to get element type
    #T = eltype(field.γ.depth)
    #println(T)
    
    # preallocate
    tracer = Array{T}(undef,size(field.γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    #println("calling one with ",T)
    tracer[field.γ.wet] .= Base.one(T) # add Base: error "should import Base"
    tracer[.!field.γ.wet] .= zero(T)/zero(T) # NaNs with right type

    d = Field(tracer,field.γ,field.name,field.longname,field.units)

    return d
end

""" 
   function oneunit, help for gridded Interpolations
"""
function one(T::Type{Field})
    TMIversion = "modern_90x45x33_GH10_GH12"
    TMIfile = download_file(TMIversion)
    γ = Grid(TMIfile)

    return TMI.ones(γ)
end

function one(field::Field{Float64})

    # use depth (could have been lon, lat)
    # to get element type
    #T = eltype(field.γ.depth)
    #println(T)
    
    # preallocate
    T = Float64
    tracer = Array{T}(undef,size(field.γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    #println("calling one with ",T)
    tracer[field.γ.wet] .= one(T) # add Base: error "should import Base"
    tracer[.!field.γ.wet] .= zero(T)/zero(T) # NaNs with right type

    d = Field(tracer,field.γ,field.name,field.longname,field.units)

    return d
end

function Field(field::Field{Float64})

    # use depth (could have been lon, lat)
    # to get element type
    #T = eltype(field.γ.depth)
    #println(T)
    
    # preallocate
    T = Float64
    tracer = Array{T}(undef,size(field.γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    # println("calling one with ",T)
    tracer[field.γ.wet] .= one(T) # add Base: error "should import Base"
    tracer[.!field.γ.wet] .= zero(T)/zero(T) # NaNs with right type

    d = Field(tracer,field.γ,field.name,field.longname,field.units)

    return d
end

function /(field::Field{Float64},scalar::Float64)

    # preallocate
    T = Float64
    tracer = Array{T}(undef,size(field.γ.wet))

    # set ocean to zero, land to NaN
    # consider whether land should be nothing or missing
    # println("calling one with ",T)
    tracer[field.γ.wet] ./= scalar # add Base: error "should import Base"
    tracer[.!field.γ.wet] .= zero(T)/zero(T) # NaNs with right type

    d = Field(tracer,field.γ,field.name,field.longname,field.units)

    return d
end

""" 
    function zeros(wet,ltype=Float64)
    initialize tracer field on TMI grid
    This version will give an array
# Arguments
- `wet`::BitArray mask of ocean points
- `ltype`:: optional type argument, default=Float64
# Output
- `d`:: 3d tracer field with NaN on dry points
"""
function zeros(wet,ltype=Float64)
    # preallocate
    d = Array{ltype}(undef,size(wet))

    # set ocean to zero, land to NaN
    d[wet] .= zero(ltype)
    d[.!wet] .= zero(ltype)/zero(ltype) # NaNs with right type
    return d
end


"""
    zero(c::Field) = zeros(c.γ)
"""
Base.zero(c::Field) = zeros(c.γ)


# Define maximum for Field to not include NaNs
Base.maximum(c::Union{Field,Source,BoundaryCondition}) = maximum(c.tracer[wet(c)])
Base.minimum(c::Union{Field,Source,BoundaryCondition}) = minimum(c.tracer[wet(c)])

"""
    Base.length(c::Union{Field,Source,BoundaryCondition}) = length(c.tracer[wet(c)])

    Extend `length` to give the number of wet (i.e., ocean) gridcells.
"""
Base.length(c::Union{Field,Source,BoundaryCondition}) = length(c.tracer[wet(c)])

"""
    Iterate over Field
"""
#Base.iterate(c::Field) =  (c.tracer[c.γ.I[1]],1)
#Base.iterate(c::Field,state) = state < length(c) ? (c.tracer[c.γ.I[state+1]],state+1) : nothing
#Base.iterate(c::Field) =  (c.tracer[wet(c)][1],1)
#Base.iterate(c::Field,state) = state < length(c) ? (c.tracer[wet(c)][state+1],state+1) : nothing

Base.iterate(c::Field) =  (c[1],1)
Base.iterate(c::Field,state) = state < length(c) ? (c[state+1],state+1) : nothing

Base.getindex(c::Field,i::Int) = c.tracer[wet(c)][i]

"""
    Specialize Base.sum(c::Field)

    so that it doesn't use the slow iteration method
"""
Base.sum(c::Field) = sum(c.tracer[wet(c)])

# """
#     Specialize Base.sum(c::BoundaryCondition)

#     so that it doesn't use the slow iteration method
# """
# Base.sum(c::BoundaryCondition) = sum(c.tracer[wet(c)])

"""
    `function \\(A,d::Field)::Field`
    Define left division for Fields
    Need two slashes to prevent invalid escape
"""
function \(A,d::Field{T})::Field{T} where T <: Real
    # initialize output
    c = zeros(d.γ,d.name,d.longname,d.units)
    c.tracer[c.γ.wet] = A\d.tracer[d.γ.wet]
    return c
end

"""
    `function +(c::BoundaryCondition,d::BoundaryCondition)::BoundaryCondition`
    Define addition for Fields
"""
function add!(c::T,d::T) where T <: Union{Source,Field,BoundaryCondition}
    if wet(c) != wet(d) # check conformability
        error("TMI type not conformable for addition")
    end
    # a strange formulation to do in-place addition
    c.tracer[wet(c)] += d.tracer[wet(d)]
end

function Base.:+(c::T,d::T) where T <: Union{Source,Field,BoundaryCondition}
    e = deepcopy(c)
    add!(e,d)
    return e
end

function subtract!(c::T,d::T) where T <: Union{Source,Field,BoundaryCondition}
    if wet(c) != wet(d) # check conformability
        error("TMI type not conformable for addition")
    end
    # a strange formulation to do in-place addition
    c.tracer[wet(c)] -= d.tracer[wet(d)]
end

function Base.:-(c::T,d::T) where T <: Union{Source,Field,BoundaryCondition}
    e = deepcopy(c)
    subtract!(e,d)
    return e
end

"""
    `function *(C,d::Field)::Field`
    Define scalar or matrix multiplication for fields

    one argument can be a number, but not both (type piracy?)
"""
# the order doesn't matter when multiplying by a scalar

Base.:*(c::Number,d::Union{Field,BoundaryCondition,Source}) = d*c
# right matrix multiply not handled
function Base.:*(c::AbstractArray,d::Union{Field,BoundaryCondition,Source})
    e = deepcopy(d)
    mul!(c,e)
    return e
end
function Base.:*(d::T,c::Union{Number,T}) where T <: Union{Field,BoundaryCondition,Source}
    e = deepcopy(d)
    mul!(e,c)
    return e
end

function mul!(d::Union{Field,BoundaryCondition,Source},C::Number)
    d.tracer[wet(d)] *= C #*d.tracer[wet(d)]
end
function mul!(c::T,d::T) where T <: Union{Field,BoundaryCondition,Source}
    # initialize output
    if wet(c) != wet(d) # check conformability
        error("Fields not conformable for multiplication")
    end

    if !isequal(d.units,c.units)
        error("Units not consistent:",d.units," vs ",c.units)
    end
    c.tracer[wet(c)] .*= d.tracer[wet(d)]
end

# matrix multiplication with arrays/matrices: order matters
# d is mutated
function mul!(c::AbstractArray,d::Union{Field,BoundaryCondition,Source})
    # initialize output
    if size(c,2) != sum(wet(d)) # check conformability
        error("Number/Matrix-Field/BC/Source multiplication: not right size")
    end
    d.tracer[wet(d)] = c*d.tracer[wet(d)]
end
function mul!(d::Union{Field,BoundaryCondition,Source},c::AbstractArray)
    # initialize output
    if size(c,1) != sum(wet(d)) # check conformability
        error("Field/BC/Source-Number/Matrix multiplication: not right size")
    end
    d.tracer[wet(d)] = d.tracer[wet(d)]*c
end

"""
    `function *(c::Field,d::Field)::Field`
    Field by field multiplication is element-by-element.
"""
function dot(c::Field{T},d::Field{T})::T where T <: Real
    e =  c.tracer[c.γ.wet]' * d.tracer[d.γ.wet]
    return e
end

"""
    function vec(u)

    Turn a collection of controls into a vector
    for use with Optim.jl. 
    An in-place version of this function would be handy.
"""
vec(u::Field) = u.tracer[u.γ.wet]
vec(u::Source) = u.tracer[u.γ.interior]
"""
    vec(u::NamedTuple)
Turn a NamedTuple of `Field`, `BoundaryCondition`, or `Source` structs into a single vector.
Skips `nothing` values. Recursively unwraps nested NamedTuples.
"""
function vec(u::NamedTuple) 
    # grab all non-nothing values
    nn = collect(filter(!isnothing, values(u)))

    # if everything is `nothing`, return an empty vector
    if isempty(nn)
        return Float64[]   # or Float64[] or whatever element type you prefer
    end

    # use the first non-nothing entry to figure out T
    ufirst = first(nn)

    if ufirst isa Union{Field, BoundaryCondition, Source}
        T = eltype(ufirst.tracer)
    elseif ufirst isa AbstractVector # Added this clause to handle NamedTuple of Vectors
        T = eltype(ufirst)
    elseif ufirst isa NamedTuple 
        #handling the case of tuples of tuples
        uufirst = first(filter(!isnothing, values(ufirst)))
        T = eltype(uufirst.tracer)
    end
    #println(typeof(ufirst)) # Commented out debug print
    # T = eltype(u)
    uvec = Vector{T}(undef,0)
    for v in u
        if !isnothing(v)
        #append!(uvec,v.tracer[v.wet])
            append!(uvec,vec(v))
        end
    end
    return uvec
end

"""
    unvec(u₀, uvec) -> u
Create a new object `u` of the same type as `u₀` (`Field`, `BoundaryCondition`, `Source`, or `NamedTuple`) and fill it with the values from `uvec`.
This is the inverse of `vec`.
# Arguments
- `u₀`: A `Field`, `BoundaryCondition`, `Source`, or `NamedTuple` to use as a template.
- `uvec`: A vector of values.
# Output
- `u`: A new object with values from `uvec`.
"""
function unvec(u₀::Union{NamedTuple,Field,BoundaryCondition, Source},uvec::Vector) #where T <: Real
    u = deepcopy(u₀)
    unvec!(u,uvec)
    return u
end

"""
    unvec(u₀::Source, uvec::Vector)
Specialization of `unvec` for `Source` types. Zeros out the tracer field before filling.
"""
function unvec(u₀::Source,uvec::Vector) #where T <: Real
    u = deepcopy(u₀)
    u.tracer .= 0.0
    unvec!(u,uvec)
    return u
end


"""
    wet(u::NamedTuple) -> Vector{Bool}

Return a concatenated boolean vector indicating wet points for all `Field`, `BoundaryCondition`, or `Source` structs within a `NamedTuple`.
Recursively handles nested `NamedTuple`s.
"""
function wet(u::NamedTuple)
    result = Bool[]
    for v in u
        if !isnothing(v)
            append!(result, wet(v))
        end
    end
    return result
end

"""
    unvec!(u, uvec)
Fill an object `u` (`Field`, `BoundaryCondition`, `Source`, or `NamedTuple`) in-place with values from `uvec`.
This is the inverse of `vec`.
# Arguments
- `u`: A `Field`, `BoundaryCondition`, `Source`, or `NamedTuple` to be modified.
- `uvec`: A vector of values.
"""
function unvec!(u::Union{BoundaryCondition{T},Field{T},Source{T}}, uvec::Vector{T}; idx::Int = 1, return_idx::Bool = false) where T <: Real
    mask = wet(u)
    data = u.tracer
    @inbounds for I in eachindex(mask)
        if mask[I]
            data[I] = uvec[idx]
            idx += 1
        end
    end
    return return_idx ? idx : nothing
end

"""
    unvec!(u::NamedTuple, uvec::Vector; ...)
In-place `unvec!` for a `NamedTuple` of `Field`, `BoundaryCondition`, or `Source` structs. Recursively fills each non-nothing element.
"""
function unvec!(u::NamedTuple, uvec::Vector; idx::Int = 1, return_idx::Bool = false) #where {N, T <: Real}
    for v in u
        if !isnothing(v)
            idx = unvec!(v, uvec; idx = idx, return_idx = true)
        end
    end
    return return_idx ? idx : nothing
end

"""
    one!(u::NamedTuple)
Set all wet points of all `Field`, `BoundaryCondition`, or `Source` structs within a `NamedTuple` to one, in-place.
"""
function one!(u::NamedTuple) #where {N, T <: Real}
    for v in u
        if !isnothing(v)
            one!(v)
        end          
    end
end

"""
    one!(u::Union{BoundaryCondition, Field, Source})
Set all wet points of a `BoundaryCondition`, `Field`, or `Source` struct to one, in-place.
"""
function one!(u::Union{BoundaryCondition, Field, Source}) #where {N, T <: Real}
    if length(u) > 1
       u.tracer[u.γ.wet] .= 1.0
    else 
        u.tracer .= 1.0
    end
end

"""
    one(u::NamedTuple) -> NamedTuple
Create a new `NamedTuple` of the same structure as `u` where all wet points of all `Field`, `BoundaryCondition`, or `Source` structs are set to one.
"""
function Base.one(u::NamedTuple) #where {N, T <: Real}
    uz = deepcopy(u)
    for v in uz
        if !isnothing(v)
            one!(v)  
        end          
    end
    return uz
end

"""
    zero!(u::NamedTuple)
Set all wet points of all `Field`, `BoundaryCondition`, or `Source` structs within a `NamedTuple` to zero, in-place. Land points are set to `NaN`.
"""
function zero!(u::NamedTuple) #where {N, T <: Real}
    for v in u
        if !isnothing(v)
            zero!(v)
        end          
    end
end

"""
    zero!(u::Union{BoundaryCondition, Field, Source})
Set all wet points of a `BoundaryCondition`, `Field`, or `Source` struct to zero, in-place. Land points are set to `NaN`.
"""
function zero!(u::Union{BoundaryCondition, Field, Source}) #where {N, T <: Real}
    if length(u) > 1
       u.tracer .= NaN
       u.tracer[wet(u)] .= 0.0
    else 
        u.tracer .= 0.0
    end
end

"""
    zero(u::NamedTuple) -> NamedTuple
Create a new `NamedTuple` of the same structure as `u` where all wet points of all `Field`, `BoundaryCondition`, or `Source` structs are set to zero.
"""
function Base.zero(u::NamedTuple) #where {N, T <: Real}
    uz = deepcopy(u)
    for v in uz
        if !isnothing(v)
            zero!(v)  
        end          
    end
    return uz
end


""" 
    function synthetic_observations(TMIversion,variable)
    Synthetic observations that are a contaminated version of real observations
    This version: gridded observations
# Arguments
- `TMIversion::String`: version of TMI water-mass/circulation model
- `variable::String`: variable name to use as template
# Output
- `y`: contaminated observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `θtrue`: real observations, 3D field
"""
function synthetic_observations(TMIversion,variable,γ)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    θtrue = readfield(TMIfile,variable,γ)
    σθ = readfield(TMIfile,"σ"*variable,γ)

    #ntrue = zeros(γ)
    #ntrue = zeros(γ.wet)
    #ntrue += rand(Normal(),length(σθ[γ.wet])) .* σθ[γ.wet]

    ntrue = Field(rand(Normal(),size(γ.wet)),γ,θtrue.name,θtrue.longname,θtrue.units)
    ntrue *= σθ
    y = θtrue + ntrue
    #y = θtrue .+ ntrue

    # get cost function (J) based on model misfit
    # here the data-model misfit is weighted by the expected error

    # weighting matrix
    #W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    W⁻ = (1/sum(γ.wet)) .* Diagonal(1 ./σθ.tracer[γ.wet].^2)
    return y, W⁻, θtrue
end
 
""" 
    function synthetic_observations(TMIversion,variable,locs)
    Synthetic observations that are a contaminated version of real observations
    This version: observations with random (uniform) spatial sampling
# Arguments
- `TMIversion::String`: version of TMI water-mass/circulation model
- `variable::String`: variable name to use as template
- `N`: number of observations
# Output
- `y`: contaminated observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `ytrue`: uncontaminated observations, 3D field
- `locs`: 3-tuples of locations for observations
- `wis`: weighted indices for interpolation to locs sites
"""
function synthetic_observations(TMIversion,variable,γ,N,σ=nothing)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    
    θtrue = readfield(TMIfile,variable,γ)
    replace!(θtrue.tracer,NaN=>0.0)

    if isnothing(σ)
        σθ = readfield(TMIfile,"σ"*variable,γ)
        replace!(σθ.tracer,NaN=>0.0)
    end

    # get random locations that are wet (ocean)
    locs = Vector{Tuple{Float64,Float64,Float64}}(undef,N)
    [locs[i] = wetlocation(γ) for i in eachindex(locs)]

    # get weighted interpolation indices
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    ytrue = observe(θtrue,wis,γ)
    if isnothing(σ)
        σtrue = observe(σθ,wis,γ)
    else
        σtrue = σ * ones(N)
    end

    #ntrue = rand(Normal.(zeros(N),σtrue),N)# .* σtrue
    ntrue = rand(Normal(),N).*σtrue
    y = ytrue .+ ntrue

    # weighting matrix
    #W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
    W⁻ = (1/N) .* Diagonal(1 ./σtrue.^2)
    return y, W⁻, θtrue, ytrue, locs, wis
end

"""
    function observe
    Take a observation at location given by weights wis
"""
function observe(c::Field{T},wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}}},γ::Grid)::Vector{T} where T <: Real

        # look at total weight, < 1 if there are land points
    # later make sure total weight = 1 for proper average
    sumwis = Vector{Float64}(undef,length(wis))
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)
    [sumwis[i] = Interpolations.InterpGetindex(wetwrap)[wis[i]...] for i in eachindex(wis)]

    # sample the true field at these random locations
    y = Vector{Float64}(undef,length(wis))
    replace!(c.tracer,NaN=>0.0)
    ywrap = view(c.tracer,list,:,:)
    [y[i] = Interpolations.InterpGetindex(ywrap)[wis[i]...]/sumwis[i] for i in eachindex(wis)]

    return y
end

"""
    function observe(c,loc,γ)

    Extend the TMI.observe method to use locations rather than weighted interpolations.
"""
function observe(c::Field{T},loc::Vector{Tuple{T,T,T}},γ::Grid) where T <: Real
    # observe at locs.
    N = length(loc)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(loc[i],γ) for i in 1:N]

    y = observe(c,wis,γ)

    return y
end

"""
    observe(c::Field{T}, wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}}}, γ::Grid{R,1}) where {T<:Real,R}

1D grid specialization of `observe` when interpolation weights have already
been computed. Normalizes by the wet-mask weight to avoid bias near land.
# Arguments
- `c`: tracer `Field`.
- `wis`: weighted interpolation indices (1D).
- `γ`: grid.
# Output
- `y`: interpolated tracer values.
"""
function observe(c::Field{T}, wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}}}, γ::Grid{R,1})::Vector{T} where {T <: Real,R}

    sumwis = Vector{Float64}(undef,length(wis))
    [sumwis[i] = Interpolations.InterpGetindex(γ.wet)[wis[i]...] for i in eachindex(wis)]

    y = Vector{Float64}(undef,length(wis))
    replace!(c.tracer,NaN=>0.0)
    [y[i] = Interpolations.InterpGetindex(c.tracer)[wis[i]...]/sumwis[i] for i in eachindex(wis)]

    return y
end

"""
    observe(c::Field{T}, loc::Vector{T}, γ::Grid{R,1}) where {T <: Real,R}
    1D grid specialization of `observe` to use locations rather than weighted interpolations.
# Arguments
- `c`: a `Field` of tracer values.
- `loc`: a vector of 1D locations.
- `γ`: TMI grid.
# Output
- `y`: observed tracer values at specified locations.
"""
function observe(c::Field{T}, loc::Vector{T}, γ::Grid{R,1}) where {T <: Real,R}
    N = length(loc)
    wis = Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}}}(undef,N)
    [wis[i] = interpindex(loc[i],γ) for i in 1:N]
    return observe(c,wis,γ)
end

"""
    function gobserve(gy::Vector{T}, c::Field{T}, locs; wis = nothing) where T <: Real

Adjoint of `observe`: back-propagate observation sensitivities `gy` to a field
using provided interpolation weights `wis` or by computing them from `locs`.
# Arguments
- `gy`: observation gradients.
- `c`: tracer `Field`.
- `locs`: observation locations (used if `wis` is not provided).
- `wis`: optional precomputed interpolation weights.
# Output
- `gc`: gradient field.
"""
function gobserve(gy::Vector{T},c::Field{T},locs; wis =nothing) where T <: Real

    #initialize gc this sneaky way
    gc = 0.0 * c
    for ii in eachindex(gy)
        # interpweights repeats some calculations
        gc.tracer[c.γ.wet] .+= gy[ii] * interpweights(locs[ii],c.γ; wis = nothing)[c.γ.wet]
    end

    return gc
end



"""
    function location_obs(field, locs, γ)
"""
function location_obs(field, locs, γ)

    tlength = length(size(field)) > 3 ? size(field)[1] : 1
    #determine weights for locations 
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    #
    sumwis = Vector{Float64}(undef,length(wis))
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)
    [sumwis[i] = Interpolations.InterpGetindex(wetwrap)[wis[i]...] for i in eachindex(wis)]
    
    field_sample = Matrix{Float64}(undef,(length(wis),tlength))
    replace!(field,NaN=>0.0)
    if tlength > 1
        for t in 1:tlength
            field_wrap = view(field[t,:,:,:], list, :, :)
            [field_sample[i,t] = field_wrap[wis[i]...]/sumwis[i] for i in eachindex(wis)]
        end

        #below code is an attempt to replace for loop 
        
        #field_wrap = view(field, :, list, :, :) #handle wraparound
        #[field_sample[i,:] = field_wrap[:,wis[i][1],wis[i][2],wis[i][3]]/sumwis[i] for i in eachindex(wis)]

    else
    
    field_wrap = view(field,list,:,:)
    [field_sample[i] = field_wrap[wis[i]...]/sumwis[i] for i in eachindex(wis)]
    end
    
    return field_sample
end

""" 
    function costfunction_gridded_obs(uvec::Vector{T},Alu,b::BoundaryCondition{T},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real

    squared model-data misfit for gridded data
    controls are a vector input for Optim.jl
# Arguments
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
- `u`: controls, field format
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary conditions
- `y`: observations on grid
- `Wⁱ`: inverse of W weighting matrix for observations
- `γ`: grid
"""
function costfunction_gridded_obs(uvec,Alu,b₀::Union{BoundaryCondition,NamedTuple},u₀::Union{BoundaryCondition,NamedTuple},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid{T}) where {T <: Real}

    # turn uvec into a boundary condition
    u = unvec(u₀,uvec)

    b = adjustboundarycondition(b₀,u) #b += u # easy case where u and b are on the same boundary
    n = steadyinversion(Alu,b,γ) - y  # gives the misfit
    J = n ⋅ (Wⁱ * n) # dot product

    # adjoint equations
    gy = -2Wⁱ * n
    gb = gsteadyinversion( gy, Alu, b, γ)
    gu = gadjustboundarycondition(gb,u)
    guvec = vec(gu)

    return J, guvec
end

"""
    function costfunction_gridded_obs!(J,guvec,uvec::Vector{T},Alu,b₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T}}}},u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}
"""
function costfunction_gridded_obs!(J,guvec,uvec::Vector{T},Alu,b₀::Union{BoundaryCondition{T},NamedTuple},u₀::Union{BoundaryCondition{T},NamedTuple},y::Field{T},Wⁱ::Diagonal{T, Vector{T}},γ::Grid{T}) where {T <: Real}

    # turn uvec into a boundary condition
    u = unvec(u₀,uvec)

    b = adjustboundarycondition(b₀,u) #b += u # easy c
    #adjustboundarycondition!(b,u) # easy case where u and b are on the same boundary
    y -= steadyinversion(Alu,b,γ)  # gives the misfit

    if guvec != nothing
        # adjoint equations
        gy = -2Wⁱ * y
        gb = gsteadyinversion( gy, Alu, b, γ)
        gu = gadjustboundarycondition(gb,u)
        #guvec = vec(gu)

        # next block just to modify the contents
        tmp = vec(gu)
        for (ii,vv) in enumerate(tmp)
            guvec[ii] = vv
        end
    end
    
    if J !=nothing
        return  y ⋅ (Wⁱ * y) # dot product
    end
end

""" 
    function costfunction_point_obs(uvec::Vector{T},Alu,b₀::BoundaryCondition{T},u₀::BoundaryCondition{T},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid;q=nothing,r=1.0) where T <: Real

    Squared model-data misfit for pointwise data.
    Controls are a vector input for Optim.jl.
    Core numerics handled by `costfunction_point_obs`.
    
# Arguments
- `uvec`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary condition
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `Q⁻`: weights for control vector
- `γ`: grid
# Optional
- `q::Field`: interior source
- `r::Number`: scalar factor for source
# Output
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
"""
function costfunction_point_obs(uvec::Vector,Alu,b::Union{BoundaryCondition,NamedTuple},u::Union{BoundaryCondition,NamedTuple},y::Vector,Wⁱ::Diagonal,wis,locs,Q⁻,γ::Grid;q=nothing,r=1.0) 

    J = 0.0
    guvec = 0.0.*uvec # same size
    #gu = unvec(u,0 .* uvec)
    J = costfunction_point_obs!(J,guvec,uvec,Alu,b,u,y,Wⁱ,wis,locs,Q⁻,γ,q₀=q,r=r)
        
    return J, guvec
end

""" 
    function costfunction_point_obs!(J,guvec,uvec::Vector{T},Alu,b₀::BoundaryCondition{T},u₀::BoundaryCondition{T},y::Vector{T},Wⁱ::Diagonal{T, Vector{T}},wis,locs,Q⁻,γ::Grid) where T <: Real

    squared model-data misfit for pointwise data
    controls are a vector input for Optim.jl
    Issue #1: couldn't figure out how to nest with costfunction_obs!
    
# Arguments
- `J`: cost function of sum of squared misfits
- `guvec`: derivative of cost function wrt to controls
- `uvec`: controls, vector format
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary condition
- `y`: pointwise observations
- `Wⁱ`: inverse of W weighting matrix for observations
- `wis`: weights for interpolation (data sampling, E)
- `locs`: data locations (lon,lat,depth)
- `Q⁻`: weights for control vector
- `γ`: grid
"""
function costfunction_point_obs!(J,guvec::Union{Nothing,Vector},uvec::Vector,Alu,b₀::Union{BoundaryCondition,NamedTuple},u₀::Union{BoundaryCondition,NamedTuple},y::Vector,Wⁱ::Diagonal,wis,locs,Q⁻,γ::Grid;q₀=nothing,r=1.0)

    u = unvec(u₀,uvec)
    b = adjustboundarycondition(b₀,u) # combine b₀, u

    if !isnothing(q₀)
        # careful with scope of c
        q = adjustsource(q₀,u)
        c = steadyinversion(Alu,b,γ,q=q,r=r)
    else
        c = steadyinversion(Alu,b,γ)
    end

    # observe at right spots
    ỹ = observe(c,wis,γ)
    n = ỹ - y

    if guvec != nothing
        ## start adjoint model

        # initialize to zero
        gu = unvec(u,0 .* uvec)
        
        gu_ = 2*(Q⁻*uvec) # control penalty gradient
        gn = 2Wⁱ * n
        gỹ = gn
        gc = gobserve(gỹ,c,locs)

        if !isnothing(q₀)
            gb,gq = gsteadyinversion(gc,Alu,b,γ,q=q,r=r)
            gadjustsource!(gu,gq,q) # pass q to linearize logscale version
        else
            gb = gsteadyinversion(gc, Alu, b, γ)
        end
        
        gadjustboundarycondition!(gu,gb)
        gu_ += vec(gu)
        for (ii,vv) in enumerate(gu_)
            guvec[ii] = vv
        end
    end

    if J !=nothing
        # control penalty and gradient
        Jcontrol = uvec'*(Q⁻*uvec)
        Jdata = n ⋅ (Wⁱ * n) # dot product
        #println("Jcontrol:",Jcontrol)
        #println("Jdata:",Jdata)
        return Jdata + Jcontrol
    end
end

""" 
    function costfunction_gridded_model(convec::Vector{T},non_zero_indices,y::Field{T},u,A0,c,q,Wⁱ::Diagonal{T, Vector{T}},Qⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real

    squared model-data misfit for gridded data
    controls are a vector input for Optim.jl
# Arguments
- `convec`: concatenated control vecotr incuding u and f
- `J`: cost function of sum of squared misfits
- `gJ`: derivative of cost function wrt to controls
- `u`: tracer controls, field format
- `non_zero_indices`: Non-zero indices for reconstruction of water-mass matrix A
- `c`: tracer concentrations from GCM
- `Wⁱ`: inverse of W weighting matrix for observations
- `Qⁱ`: inverse of Q weighting matrix for tracer conservation
- `γ`: grid
"""
function costfunction_gridded_model(convec,non_zero_indices,u₀::Field{T},A0,y::Vector{T},c,q,Wⁱ::Diagonal{T, Vector{T}},Qⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real
    ulength = sum(γ.wet)
    
    #control vectors
    uvec=convec[begin:ulength]
    ufvec = convec[ulength+1:end]

    Actl = sparse(non_zero_indices[:, 1], non_zero_indices[:, 2], ufvec)
    A=A0 + Actl
    dummy,dummy,fguess  = findnz(A0)
    dummy,dummy,fnow  = findnz(A) 
    onesvec = ones(size(q))
    csum = Wⁱ * uvec+c


    # find lagrange multipliers
    muk = transpose(A) * Qⁱ * (A * c - q)
    dAcdf = spzeros(length(ufvec),length(c))
    dA1df = spzeros(length(ufvec),length(c))
    for ii in eachindex(ufvec)
          dAcdf[ii,non_zero_indices[ii, 1]] = c[non_zero_indices[ii, 2]]
          dA1df[ii,non_zero_indices[ii, 2]] = 1
    end
    

    dAdf_terms = dAcdf * Qⁱ * (A * c - q) + dA1df * Qⁱ * (A * onesvec - onesvec)

    J =  uvec ⋅ uvec + transpose(A * c - q) * Qⁱ * (A*c - q) - 
          2 * transpose(muk)*( Wⁱ * uvec+c-y) +
          transpose(A * onesvec - onesvec) * Qⁱ * (A* onesvec - onesvec)

    # adjoint equations
    guvec = zeros(length(convec))

    for (ii,vv) in enumerate(convec)
        if ii <= ulength 
          #this is the derivative of the cost function wrt the part of the control vector
          # associated with the tracer concentration
          guvec[ii] =  2 * uvec[ii] - (2 * transpose(muk) * Wⁱ)[ii]
        else
          #this is the derivative of the cost function wrt the part of the control vector
          # associated with the transport vector
          guvec[ii]=2 * dAdf_terms[ii-ulength]#2 * convec[ii]+
        end
    end

    return J , guvec
end

"""
    function costfunction_gridded_model!(J,guvec,convec::Vector{T},non_zero_indices,u₀::Union{BoundaryCondition{T},NamedTuple{<:Any, NTuple{N2,BoundaryCondition{T}}}},c,y::Field{T},Wⁱ::Diagonal{T, Vector{T}},Qⁱ::Diagonal{T, Vector{T}},γ::Grid) where {N1, N2, T <: Real}
"""
function costfunction_gridded_model!(J,guvec,convec::Vector{T},non_zero_indices,u₀::Field{T},A0,y::Vector{T},c,q,Wⁱ::Diagonal{T, Vector{T}},Qⁱ::Diagonal{T, Vector{T}},γ::Grid) where T <: Real

    ulength = sum(γ.wet)
    uvec = convec[begin:ulength]
    ufvec = convec[ulength+1:end]
    
    Actl = sparse(non_zero_indices[:, 1], non_zero_indices[:, 2], ufvec)
    A=A0 + Actl
    dummy,dummy,fguess  = findnz(A0)
    dummy,dummy,fnow  = findnz(A)
    onesvec = ones(size(q))
    csum = Wⁱ * uvec+c

    # find lagrange multipliers
    muk = transpose(A) * Qⁱ * (A * c - q)
    dAcdf = spzeros(length(ufvec),length(c))
    dA1df = spzeros(length(ufvec),length(c))
    for ii in eachindex(ufvec)
          dAcdf[ii,non_zero_indices[ii, 1]] = c[non_zero_indices[ii, 2]]
          dA1df[ii,non_zero_indices[ii, 1]] = 1
    end
    dAdf_terms = dAcdf * Qⁱ * (A * c - q) + dA1df * Qⁱ * (A * onesvec - onesvec)

    if guvec != nothing
        tmp = guvec
        for (ii,vv) in enumerate(tmp)
            if ii <= ulength
               guvec[ii] = 2 * uvec[ii]-(2 * transpose(muk) * Wⁱ)[ii]
            else
               guvec[ii]=2 * dAdf_terms[ii-ulength]#2 * convec[ii] +
            end
        end
    end
    
    if J !=nothing
        return uvec ⋅ uvec + transpose(A * c - q) * Qⁱ * (A*c - q)-
                  2 * transpose(muk)*( Wⁱ * uvec+c-y)+
               transpose(A * onesvec - onesvec) * Qⁱ * (A* onesvec - onesvec)
    end
end

""" 
    function steadyinversion(Alu,b;q=nothing,r=1.0)
    invert for a steady-state tracer distribution
# Arguments
- `Alu`: LU decomposition of water-mass matrix
- `b`: boundary condition, assumed to be surface boundary condition
- `γ`::Grid
# Optional Arguments
- `q`: interior sources/sinks of phosphate
- `r`: stochiometric ratio of tracer:phosphate
# Output
- `c`::Field, steady-state tracer distribution
"""
function steadyinversion(Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU, SparseMatrixCSC},b::BoundaryCondition,γ::Grid{T};q=nothing,r=1.0)::Field{T} where T <: Real

    # preallocate Field for equation constraints
    d = zeros(γ,b.name,b.longname,b.units)
    
    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end

    c = Alu \ d
    return c
end


""" 
    function gsteadyinversion(gc,Alu,b;q=nothing,r=1.0)

    ADJOINT invert for a steady-state tracer distribution

# Arguments
- `Alu`: LU decomposition of water-mass matrix
- `b`: BoundaryCondition
- `γ`::Grid
# Optional Arguments
- `q`: interior sources/sinks of phosphate
- `r`: stochiometric ratio of tracer:phosphate
# Output
- `c`::Field, steady-state tracer distribution
"""
function gsteadyinversion(gc::Field,Alu,b::Union{BoundaryCondition,NamedTuple},γ::Grid;q=nothing,r=1.0) #where T <: Real
    #println("running adjoint steady inversion")
    gd = Alu' \ gc
    gb = gsetboundarycondition(gd,b)

    if !isnothing(q)
        gq = gsetsource(gd,q,r)
        return gb,gq
    else
        return gb
    end

end


# function steadyinversion(Alu,b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T,R,N2,B}}},γ::Grid{R};q=nothing,r=1.0)::Field{R} where {N1, N2 <: Integer, T <: Real, R <: Real, B <: AbstractMatrix{T}}
function steadyinversion(Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU, SparseMatrixCSC},b::NamedTuple,γ::Grid{T};q=nothing,r=1.0)::Field{T} where {T <: Real}

    # preallocate Field for equation constraints
    d = zeros(γ,first(b).name,first(b).longname,first(b).units)

    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end

    # define ldiv with fields
    #c = zeros(d.γ,b.name,b.longname,b.units)
    c = Alu \ d

    return c
end

"""
    steadyinversion(Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU, SparseMatrixCSC},
                    b::NamedTuple{tracer_names}, q::NamedTuple, γ; r=1.0)

Multi-tracer steady inversion using an LU/UMFPACK factorization of the
watermassmatrix. Applies per-tracer boundary conditions and sources, returning
a `NamedTuple` of steady tracer fields.
# Arguments
- `Alu`: LU/UMFPACK factorization of the watermassmatrix.
- `b`: `NamedTuple` of boundary conditions keyed by tracer.
- `q`: `NamedTuple` of interior sources keyed by tracer.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- `c_nt`: `NamedTuple` of steady tracer `Field`s.
"""
function steadyinversion(Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU, SparseMatrixCSC},
                        b::NamedTuple{tracer_names, S}, q::NamedTuple,
                        γ::Grid{T}; r=1.0) where {T <: Real, tracer_names, S}
                        
    c_nt = map(b, q) do b_i, q_i
            steadyinversion(Alu, b_i, γ; q = q_i, r = r)
            end
    return c_nt

end


"""
    lsolve(cache::LinearSolve.LinearCache, d::Field, γ) -> Field

Solve `A * c = d` with LinearSolve (using the supplied `LinearCache`, which may
wrap a preconditioner such as IncompleteLU). Inserts the solution into a new
`Field`, placing values only on wet points.
# Arguments
- `cache`: `LinearSolve.LinearCache` containing operator/preconditioner.
- `d`: right-hand-side `Field`.
- `γ`: grid.
# Output
- `c`: solution `Field`.
"""
function lsolve(cache::LinearSolve.LinearCache, d::Field{T}, γ) where T
    c = zeros(γ,d.name,d.longname,d.units)
    cache.b .= vec(d)
    sol = solve!(cache)
    c.tracer[wet(c)] .= sol.u
    return c
end

"""
    steadyinversion(cache::LinearSolve.LinearCache, b::BoundaryCondition, γ; q=nothing, r=1.0)

Steady-state tracer inversion solved with LinearSolve/Iterative solvers (via
the provided `LinearCache`, often with an ILU preconditioner) instead of a
direct `\\` on an LU factorization. Applies boundary conditions `b`, optional
interior sources `q` scaled by `r`, then solves for the steady tracer field.
# Arguments
- `cache`: `LinearSolve.LinearCache` (operator and preconditioner).
- `b`: surface `BoundaryCondition`.
- `γ`: grid.
# Optional Arguments
- `q`: interior sources/sinks.
- `r`: stoichiometric ratio applied to `q`.
# Output
- `c`: steady-state tracer `Field`.
"""
function steadyinversion(cache::LinearSolve.LinearCache,b::BoundaryCondition,γ::Grid{T};q=nothing,r=1.0)::Field{T} where T <: Real

    # preallocate Field for equation constraints
    d = zeros(γ,b.name,b.longname,b.units)
    
    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end

    return lsolve(cache, d, γ)
end



# function steadyinversion(Alu,b::NamedTuple{<:Any, NTuple{N1,BoundaryCondition{T,R,N2,B}}},γ::Grid{R};q=nothing,r=1.0)::Field{R} where {N1, N2 <: Integer, T <: Real, R <: Real, B <: AbstractMatrix{T}}
"""
    steadyinversion(cache::LinearSolve.LinearCache, b::NamedTuple, γ; q=nothing, r=1.0)

NamedTuple overload of the LinearSolve-based steady inversion; applies each
boundary condition and optional source before solving with the provided
`LinearCache` (e.g., iterative solve with ILU preconditioner).
# Arguments
- `cache`: `LinearSolve.LinearCache`.
- `b`: `NamedTuple` of boundary conditions.
- `γ`: grid.
# Optional Arguments
- `q`: `NamedTuple` of interior sources/sinks.
- `r`: stoichiometric ratio applied to `q`.
# Output
- `c`: steady-state tracer `Field`.
"""
function steadyinversion(cache::LinearSolve.LinearCache,b::NamedTuple,γ::Grid{T};q=nothing,r=1.0)::Field{T} where {T <: Real}

    # preallocate Field for equation constraints
    d = zeros(γ,first(b).name,first(b).longname,first(b).units)

    # update d with the boundary condition b
    setboundarycondition!(d,b)

    if !isnothing(q)
        # apply interior sources
        # negative because of equation arrangement
        setsource!(d,q,r)
    end

    return lsolve(cache, d, γ)

end

"""
    steadyinversion(cache::LinearSolve.LinearCache, b::NamedTuple{tracer_names}, q::NamedTuple,
                    γ; r=1.0, c_obs=nothing)

Multi-tracer steady inversion using LinearSolve/Iterative methods. Optionally
modifies the operator per tracer (e.g., subtracting a decay matrix) and builds
a new `LinearCache` with ILU preconditioning when needed. Returns a
`NamedTuple` of steady tracer fields.
# Arguments
- `cache`: `LinearSolve.LinearCache` (base operator/preconditioner).
- `b`: `NamedTuple` of boundary conditions keyed by tracer.
- `q`: `NamedTuple` of interior sources keyed by tracer.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
- `c_obs`: optional observations that may supply a per-tracer `decay_rate_matrix`.
# Output
- `c_nt`: `NamedTuple` of steady tracer `Field`s.
"""
function steadyinversion(cache::LinearSolve.LinearCache,b::NamedTuple{tracer_names, S}, q::NamedTuple,
                        γ::Grid{T}; r=1.0, c_obs = nothing) where {T <: Real, tracer_names, S} # c_obs is now OPTIONAL
    n_tracers = length(tracer_names)
    c_results = Vector{Field}(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        local_cache = cache # Start with the original cache

        # Check for decay rate for this specific tracer, with !isnothing(c_obs)
        if !isnothing(c_obs) && haskey(c_obs, name) && !isnothing(c_obs[name].decay_rate)
            @error "Decay rates in steady inversion not implemented yet"
            # A_modified = cache.A - c_obs[name].decay_rate_matrix # Modify A from the original cache
            # prob = LinearProblem(A_modified, Vector{Float64}(undef, size(A_modified,1)))
            # P = ilu(A_modified; τ=0.01)
            # local_cache = init(prob, KrylovJL_GMRES(); Pl=P) # Create a new local_cache
        end

        b_i = get(b, name, nothing)
        q_i = get(q, name, nothing)  # q is a NamedTuple
        c_results[i] = steadyinversion(local_cache, b_i, γ; q = q_i, r = r) # Use the potentially modified local_cache
    end
    c_nt = NamedTuple{tracer_names}(Tuple(c_results))
    return c_nt

end

"""
    gsteadyinversion(gc::Field, c::Field, A, Alu, b, q, γ; r=1.0)

Adjoint of `steadyinversion` for a single tracer using a direct factorization
(`lu`/UMFPACK). Returns gradients with respect to boundary conditions `b`,
optional sources `q`, and the watermassmatrix `A`.
# Arguments
- `gc`: adjoint `Field` of the solution.
- `c`: forward solution `Field`.
- `A`: watermassmatrix.
- `Alu`: LU/UMFPACK factorization of `A`.
- `b`: boundary condition(s) used in the forward solve.
- `q`: optional interior sources/sinks.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- `gb`, `gq`, `gA`: gradients w.r.t. boundary, sources, and `A` (with `gq` possibly `nothing`).
"""
function gsteadyinversion(gc::Field,c::Field,A, Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU}, 
                        b::Union{BoundaryCondition,NamedTuple}, q::Union{Nothing,Source},γ::Grid; r=1.0) #where T <: Real
    #println("running adjoint steady inversion")
    gd = Alu' \ gc

    #ga = gd * c', but some A entries are fixed (1 or 0)
    #we just need to consider the entries that propogate into the adjoint
    rows, cols, _ = findnz(A) 
    offdiag_mask = rows .!= cols 
    #may need to call the other steady inversion
    # only consider the off-diagonal entries
    rows_offdiag = rows[offdiag_mask]
    cols_offdiag = cols[offdiag_mask]

    vals_offdiag = -vec(gd)[rows_offdiag] .* vec(c)[cols_offdiag]
    gA = sparse(rows_offdiag, cols_offdiag, vals_offdiag, size(A,1), size(A,2))

    gb = gsetboundarycondition(gd,b)  #update bc based on preexisting object

    if !isnothing(q)
        gq = gsetsource(gd,q,r) #update source based on preexisting object
        return gb,gq, gA
    else
        return gb, nothing, gA
    end
end

"""
    gsteadyinversion(gc::NamedTuple, c::NamedTuple, A, Alu, b, q, γ; r=1.0)

Multi-tracer adjoint of `steadyinversion` using a direct factorization
(`lu`/UMFPACK). Aggregates `gA` across tracers and returns per-tracer gradients
for boundary conditions and sources.
# Arguments
- `gc`: `NamedTuple` of adjoint fields.
- `c`: `NamedTuple` of forward solutions.
- `A`: watermassmatrix.
- `Alu`: LU/UMFPACK factorization of `A`.
- `b`: `NamedTuple` of boundary conditions.
- `q`: `NamedTuple` of interior sources/sinks.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- `gb`, `gq`, `gA`: per-tracer boundary/source gradients and aggregated `gA`.
"""
function gsteadyinversion(gc::NamedTuple, c::NamedTuple, A, Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU}, 
                          b::NamedTuple{b_names, B}, q::NamedTuple{q_names, Q}, 
                          γ::Grid; r=1.0) where {b_names, B, q_names, Q}

    tracer_names = keys(gc)
    gA_total = spzeros(size(A,1), size(A,2))

    #should i preallocate here?
    gA_rows = Int[]
    gA_cols = Int[]
    gA_vals = Float64[]

    n_tracers = length(tracer_names)
    gb_results = Vector{Union{B.parameters...}}(undef, n_tracers)
    gq_results = Vector{Union{Q.parameters...}}(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)  # Added enumerate
        b_i = get(b, name, nothing)
        q_i = get(q, name, nothing)  # q is a NamedTuple

        gb_i, gq_i, gA_i = gsteadyinversion(gc[name], c[name], A, Alu, b_i, q_i, γ; r=r)
        gb_results[i] = gb_i
        gq_results[i] = gq_i

        # appending will be equivalent to gA_total += gA_i
        rows, cols, vals = findnz(gA_i)
        append!(gA_rows, rows)
        append!(gA_cols, cols)
        append!(gA_vals, vals)
    end

    gA_total = sparse(gA_rows, gA_cols, gA_vals, size(A,1), size(A,2))

    # gb = NamedTuple{tracer_names}(Tuple(gb_results))
    # gq = NamedTuple{tracer_names}(Tuple(gq_results))
    gb = (;zip(tracer_names, gb_results)...)
    gq = (;zip(tracer_names, gq_results)...)
    return gb, gq, gA_total
end


"""
    gsteadyinversion!(gA_vals, gc, c, Alu, b, q, rows_offdiag, cols_offdiag; r=1.0)

In-place adjoint of `steadyinversion` using a direct factorization. Accumulates
`gA_vals` on the off-diagonal structure of the watermassmatrix and returns
boundary/source gradients.
# Arguments
- `gA_vals`: buffer for gradients of the watermassmatrix (off-diagonal entries).
- `gc`: adjoint `Field` of the solution.
- `c`: forward solution `Field`.
- `Alu`: LU/UMFPACK factorization of the watermassmatrix.
- `b`: boundary condition(s).
- `q`: optional interior sources/sinks.
- `rows_offdiag`, `cols_offdiag`: off-diagonal structure of the watermassmatrix.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- `gb`, `gq`: gradients for boundary conditions and sources.
"""
function gsteadyinversion!(
    gA_vals::AbstractVector,     # preallocated, length == length(rows_offdiag)
    gc::Field, c::Field, Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU},
    b::Union{BoundaryCondition,NamedTuple},
    q::Union{Nothing,Source},
    rows_offdiag::AbstractVector{<:Integer},
    cols_offdiag::AbstractVector{<:Integer};
    r::Real = 1.0,
)
    wet_inds = findall(wet(gc))

    # Solve adjoint system
    gd = Alu' \ gc    # later you could make this ldiv! if you have a workspace

    @inbounds @simd for k in eachindex(rows_offdiag)
        ip = rows_offdiag[k]     # row index of wet point 
        jp = cols_offdiag[k]     # col index of wet point 

        il = wet_inds[ip]        # full-grid index for gd
        jl = wet_inds[jp]        # full-grid index for c

        gA_vals[k] += -gd.tracer[il] * c.tracer[jl]
    end

    gb = gsetboundarycondition(gd, b)
    gq = isnothing(q) ? nothing : gsetsource(gd, q, r)

    return gb, gq
end

function gsteadyinversion(gc::NamedTuple, c::NamedTuple, A, Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU},
    b::NamedTuple{b_names, B}, q::NamedTuple{q_names, Q}; r = 1.0) where {b_names, B, q_names, Q}

    tracer_names = keys(gc)

    # Get off-diagonal structure of A once
    rows, cols, _ = findnz(A)
    offdiag_mask = rows .!= cols
    gA_rows = rows[offdiag_mask]
    gA_cols = cols[offdiag_mask]

    # Allocate gA_vals, to be accumulated by each tracer
    Tval = eltype(gc[first(tracer_names)].tracer)
    gA_vals = zeros(Tval, length(gA_rows))

    n_tracers = length(tracer_names)
    gb_results = Vector{Union{B.parameters...}}(undef, n_tracers)
    gq_results = Vector{Union{Q.parameters...}}(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        b_i = get(b, name, nothing)
        q_i = get(q, name, nothing)

        gb_i, gq_i = gsteadyinversion!(
            gA_vals, gc[name], c[name], Alu,
            b_i, q_i, gA_rows, gA_cols; r = r)

        gb_results[i] = gb_i
        gq_results[i] = gq_i
    end

    gA_total = sparse(gA_rows, gA_cols, gA_vals, size(A,1), size(A,2))

    gb = (; zip(tracer_names, gb_results)...)
    gq = (; zip(tracer_names, gq_results)...)

    return gb, gq, gA_total
end


"""
    gsteadyinversion!(gA_vals, gc, c, cache, b, q, rows_offdiag, cols_offdiag, γ; r=1.0)

In-place adjoint of `steadyinversion` solved with LinearSolve (e.g., iterative +
ILU). Uses precomputed off-diagonal structure of the watermassmatrix to
accumulate `gA_vals` and returns boundary/source gradients.
# Arguments
- `gA_vals`: buffer for gradients of the watermassmatrix (off-diagonal entries).
- `gc`: adjoint `Field` of the solution.
- `c`: forward solution `Field`.
- `cache`: `LinearSolve.LinearCache` for the watermassmatrix (or variant).
- `b`: boundary condition(s).
- `q`: optional interior sources/sinks.
- `rows_offdiag`, `cols_offdiag`: off-diagonal structure indices of the watermassmatrix.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- `gb`, `gq`: gradients for boundary conditions and sources.
"""
function gsteadyinversion!(
    gA_vals::AbstractVector,     # preallocated, length == length(rows_offdiag)
    gc::Field, c::Field, 
    cache::LinearSolve.LinearCache,
    b::Union{BoundaryCondition,NamedTuple},
    q::Union{Nothing,Source},
    rows_offdiag::AbstractVector{<:Integer},
    cols_offdiag::AbstractVector{<:Integer}, 
    γ::Grid;
    r::Real = 1.0,
)
    wet_inds = findall(wet(gc))

    # Solve adjoint system
    gd = lsolve(cache, gc, γ) #Alu' \ gc    # later you could make this ldiv! if you have a workspace

    @inbounds @simd for k in eachindex(rows_offdiag)
        ip = rows_offdiag[k]     # row index of wet point 
        jp = cols_offdiag[k]     # col index of wet point 

        il = wet_inds[ip]        # full-grid index for gd
        jl = wet_inds[jp]        # full-grid index for c

        gA_vals[k] += -gd.tracer[il] * c.tracer[jl]
    end

    gb = gsetboundarycondition(gd, b)
    gq = isnothing(q) ? nothing : gsetsource(gd, q, r)

    return gb, gq
end

"""
    gsteadyinversion(gc::NamedTuple, c::NamedTuple, A, cache::LinearSolve.LinearCache,
                     b::NamedTuple, q::NamedTuple, γ; r=1.0, c_obs)

Multi-tracer adjoint for a LinearSolve-based steady inversion. Supports
per-tracer decay via `c_obs[name].decay_rate_matrix` by building tracer-specific
adjoint caches. Aggregates `gA` across tracers.
# Arguments
- `gc`: `NamedTuple` of adjoint fields.
- `c`: `NamedTuple` of forward solutions.
- `A`: watermassmatrix.
- `cache`: `LinearSolve.LinearCache` (base adjoint operator/preconditioner).
- `b`: `NamedTuple` of boundary conditions.
- `q`: `NamedTuple` of interior sources/sinks.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
- `c_obs`: observation info containing optional `decay_rate_matrix` per tracer.
# Output
- `gb`, `gq`, `gA`: per-tracer boundary/source gradients and aggregated `gA`.
"""
function gsteadyinversion(gc::NamedTuple, c::NamedTuple, A, cache::LinearSolve.LinearCache,
    b::NamedTuple{b_names, B}, q::NamedTuple{q_names, Q}, γ::Grid; r = 1.0, c_obs) where {b_names, B, q_names, Q} # c_obs is now REQUIRED

    tracer_names = keys(gc)

    # Get off-diagonal structure of A once
    rows, cols, _ = findnz(A)
    offdiag_mask = rows .!= cols
    gA_rows = rows[offdiag_mask]
    gA_cols = cols[offdiag_mask]

    # Allocate gA_vals, to be accumulated by each tracer
    Tval = eltype(gc[first(tracer_names)].tracer)
    gA_vals = zeros(Tval, length(gA_rows))

    n_tracers = length(tracer_names)
    gb_results = Vector{Union{B.parameters...}}(undef, n_tracers)
    gq_results = Vector{Union{Q.parameters...}}(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        local_cache_adjoint = cache # Start with the original adjoint cache
        if haskey(c_obs, name) && !isnothing(c_obs[name].decay_rate)
            @error "Decay rates in gsteadyinversion not yet implemented" 
            # A_modified = A - c_obs[name].decay_rate_matrix
            # A_modified_t = sparse(transpose(A_modified))
            # prob_t = LinearProblem(A_modified_t, Vector{Float64}(undef, size(A_modified_t,1)))
            # P_t = ilu(A_modified_t; τ=0.01)
            # local_cache_adjoint = init(prob_t, KrylovJL_GMRES(); Pl=P_t)
        end

        b_i = get(b, name, nothing)
        q_i = get(q, name, nothing)

        gb_i, gq_i = gsteadyinversion!(
            gA_vals, gc[name], c[name], local_cache_adjoint, # Pass local_cache_adjoint here
            b_i, q_i, gA_rows, gA_cols, γ; r = r)

        gb_results[i] = gb_i
        gq_results[i] = gq_i
    end

    gA_total = sparse(gA_rows, gA_cols, gA_vals, size(A,1), size(A,2))

    gb = (; zip(tracer_names, gb_results)...)
    gq = (; zip(tracer_names, gq_results)...)

    return gb, gq, gA_total
end


"""
    gsteadyinversion!(gdub, gduq, gA_vals, gc, c, Alu, b, q, rows_offdiag, cols_offdiag, γ; r=1.0)

In-place adjoint of `steadyinversion` using a direct LU/UMFPACK factorization of
the watermassmatrix. Accumulates `gA_vals` on the off-diagonal structure and
writes boundary/source adjoints into provided buffers.
# Arguments
- `gdub`: destination buffer for boundary-condition adjoints (`BoundaryCondition` or `NamedTuple`, or `nothing` to skip).
- `gduq`: destination buffer for source adjoints (`Source`/`NamedTuple`, or `nothing`).
- `gA_vals`: buffer for gradients of the watermassmatrix (off-diagonal entries).
- `gc`: adjoint `Field` of the solution.
- `c`: forward solution `Field`.
- `Alu`: LU/UMFPACK factorization of the watermassmatrix.
- `b`: boundary condition(s).
- `q`: optional interior sources/sinks.
- `rows_offdiag`, `cols_offdiag`: off-diagonal structure of the watermassmatrix.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- None (buffers mutated in-place).
"""
function gsteadyinversion!(
    gdub::Union{Nothing,BoundaryCondition,NamedTuple},
    gduq::Union{Nothing,Source,NamedTuple},
    gA_vals::AbstractVector,
    gc::Field, c::Field, 
    cache::LinearSolve.LinearCache,
    b::Union{BoundaryCondition,NamedTuple},
    q::Union{Nothing,Source},
    rows_offdiag::AbstractVector{<:Integer},
    cols_offdiag::AbstractVector{<:Integer}, 
    γ::Grid;
    r::Real = 1.0,
)
    wet_inds = findall(wet(gc))

    gd = lsolve(cache, gc, γ)

    @inbounds @simd for k in eachindex(rows_offdiag)
        ip = rows_offdiag[k]
        jp = cols_offdiag[k]

        il = wet_inds[ip]
        jl = wet_inds[jp]

        gA_vals[k] += -gd.tracer[il] * c.tracer[jl]
    end

    if !isnothing(gdub)
        adjustboundarycondition!(gdub, gsetboundarycondition(gd, b); r = 1.0)
    end
    if !isnothing(gduq) && !isnothing(q)
        adjustsource!(gduq, gsetsource(gd, q, r); r = 1.0)
    end

    return nothing
end


"""
    gsteadyinversion!(gb_dest, gq_dest, gA_vals, gc, c, A, cache, b, q, γ; r=1.0)

In-place variant that writes boundary/source adjoints directly into preallocated
`gb_dest`/`gq_dest` buffers (e.g., `controls.gdub`/`controls.gduq`). Accumulates
`gA_vals` in the same way as the allocating version. Does not return `gb_dest`/`gq_dest`.
"""
function gsteadyinversion!(
    gdub::NamedTuple, gduq::NamedTuple,
    gc::NamedTuple, c::NamedTuple, A, cache::LinearSolve.LinearCache,
    b::NamedTuple{b_names, B}, q::NamedTuple{q_names, Q}, γ::Grid; r = 1.0, c_obs
) where {b_names, B, q_names, Q} 

    tracer_names = keys(gc)

    rows, cols, _ = findnz(A)
    offdiag_mask = rows .!= cols
    gA_rows = rows[offdiag_mask]
    gA_cols = cols[offdiag_mask]

    Tval = eltype(gc[first(tracer_names)].tracer)
    gA_vals = zeros(Tval, length(gA_rows))

    for name in tracer_names
        local_cache_adjoint = cache # Start with the original adjoint cache
        if haskey(c_obs, name) && !isnothing(c_obs[name].decay_rate)
            @error "Decay rates in gsteadyinversion not yet impleme
            # A_modified = A - c_obs[name].decay_rate_matrix
            # A_modified_t = sparse(transpose(A_modified))
            # prob_t = LinearProblem(A_modified_t, Vector{Float64}(undef, size(A_modified_t,1)))
            # P_t = ilu(A_modified_t; τ=0.01)
            # local_cache_adjoint = init(prob_t, KrylovJL_GMRES(); Pl=P_t)
        end

        b_i = get(b, name, nothing)
        q_i = get(q, name, nothing)

        gsteadyinversion!(
            gdub[name], gduq[name],
            gA_vals, gc[name], c[name], local_cache_adjoint, # Pass local_cache_adjoint here
            b_i, q_i, gA_rows, gA_cols, γ; r = r)
    end

    gA_total = sparse(gA_rows, gA_cols, gA_vals, size(A,1), size(A,2))
    return gA_total
end

"""
    gsteadyinversion!(gdub, gduq, gA_vals, gc, c, Alu, b, q, rows_offdiag, cols_offdiag, γ; r=1.0)

In-place adjoint of `steadyinversion` using a direct LU/UMFPACK factorization of
the watermassmatrix. Accumulates `gA_vals` on the off-diagonal structure and
writes boundary/source adjoints into provided buffers.
# Arguments
- `gdub`: destination buffer for boundary-condition adjoints (`BoundaryCondition`/`NamedTuple` or `nothing`).
- `gduq`: destination buffer for source adjoints (`Source`/`NamedTuple` or `nothing`).
- `gA_vals`: buffer for gradients of the watermassmatrix (off-diagonal entries).
- `gc`: adjoint `Field` of the solution.
- `c`: forward solution `Field`.
- `Alu`: LU/UMFPACK factorization of the watermassmatrix.
- `b`: boundary condition(s).
- `q`: optional interior sources/sinks.
- `rows_offdiag`, `cols_offdiag`: off-diagonal structure of the watermassmatrix.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
# Output
- None (buffers mutated in-place).
"""
function gsteadyinversion!(
    gdub::Union{Nothing,BoundaryCondition,NamedTuple},
    gduq::Union{Nothing,Source,NamedTuple},
    gA_vals::AbstractVector,
    gc::Field, c::Field, 
    Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU},
    b::Union{BoundaryCondition,NamedTuple},
    q::Union{Nothing,Source},
    rows_offdiag::AbstractVector{<:Integer},
    cols_offdiag::AbstractVector{<:Integer}, 
    γ::Grid;
    r::Real = 1.0,
)
    wet_inds = findall(wet(gc))

    gd = Alu' \  gc

    @inbounds @simd for k in eachindex(rows_offdiag)
        ip = rows_offdiag[k]
        jp = cols_offdiag[k]

        il = wet_inds[ip]
        jl = wet_inds[jp]

        gA_vals[k] += -gd.tracer[il] * c.tracer[jl]
    end

    if !isnothing(gdub)
        adjustboundarycondition!(gdub, gsetboundarycondition(gd, b); r = 1.0)
    end
    if !isnothing(gduq) && !isnothing(q)
        adjustsource!(gduq, gsetsource(gd, q, r); r = 1.0)
    end

    return nothing
end


"""
    gsteadyinversion!(gdub, gduq, gA_vals, gc, c, A, Alu, b, q, γ; r=1.0, c_obs)

Multi-tracer in-place adjoint using a direct LU/UMFPACK factorization of the
watermassmatrix. Accumulates `gA_vals` over off-diagonal entries and writes
per-tracer boundary/source adjoints into provided buffers.
# Arguments
- `gdub`: `NamedTuple` of boundary-condition adjoint buffers.
- `gduq`: `NamedTuple` of source adjoint buffers.
- `gA_vals`: buffer for gradients of the watermassmatrix (off-diagonal entries).
- `gc`: `NamedTuple` of adjoint fields.
- `c`: `NamedTuple` of forward solutions.
- `A`: watermassmatrix.
- `Alu`: LU/UMFPACK factorization of `A`.
- `b`: `NamedTuple` of boundary conditions.
- `q`: `NamedTuple` of interior sources/sinks.
- `γ`: grid.
# Optional Arguments
- `r`: stoichiometric ratio applied to `q`.
- `c_obs`: unused (kept for API symmetry with LinearSolve variant).
# Output
- `gA_total`: sparse matrix of accumulated watermassmatrix gradients.
"""
function gsteadyinversion!(
    gdub::NamedTuple, gduq::NamedTuple,
    gc::NamedTuple, c::NamedTuple, A, Alu::Union{LU, SparseArrays.UMFPACK.UmfpackLU},
    b::NamedTuple{b_names, B}, q::NamedTuple{q_names, Q}, γ::Grid; r = 1.0, c_obs
) where {b_names, B, q_names, Q} 

    tracer_names = keys(gc)

    rows, cols, _ = findnz(A)
    offdiag_mask = rows .!= cols
    gA_rows = rows[offdiag_mask]
    gA_cols = cols[offdiag_mask]

    Tval = eltype(gc[first(tracer_names)].tracer)
    gA_vals = zeros(Tval, length(gA_rows))

    for name in tracer_names
        b_i = get(b, name, nothing)
        q_i = get(q, name, nothing)

        gsteadyinversion!(
            gdub[name], gduq[name],
            gA_vals, gc[name], c[name], Alu, # Pass local_cache_adjoint here
            b_i, q_i, gA_rows, gA_cols, γ; r = r)
    end

    gA_total = sparse(gA_rows, gA_cols, gA_vals, size(A,1), size(A,2))
    return gA_total
end


# """
#     function gsteadyinversion(gc::Field{T},Alu,b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}},γ::Grid;q=nothing,r=1.0)::Field{T} where {N, T <: Real}

#     ADDJOINT steady inversion for b::NamedTuple
# """
# function gsteadyinversion(gc::Field{T},Alu,b::NamedTuple{<:Any, NTuple{N,BoundaryCondition{T}}},γ::Grid;q=nothing,r=1.0) where {N, T <: Real}

#     #println("running adjoint steady inversion")

#     # sensitivity of Field of equation constraints
#     gd = Alu' \ gc

#     # still need to develop this
#     # and to find a way to return the value
#     if !isnothing(q)
#         gq = gsetsource!(gd,d,q,r)
#     end

#     # update d with the boundary condition b
#     gb = gsetboundarycondition(gd,b)

#     return gb
# end

"""
    function wetlocation(γ)
    Get (lon,lat,depth) tuples of wet locations.
    Allow a location to be wet if at least one out of 8 nearby gridpoints is wet.
    Certainly "wet" gridpoints could be defined more strictly.
# Arguments
- `γ`: TMI.grid
# Output
- `loc`: lon,lat,depth """
function wetlocation(γ)

    confirmwet = false
    neighbors  = 8
    while !confirmwet
        loc = (rand(minimum(γ.lon):0.1:maximum(γ.lon)),
               rand(minimum(γ.lat):0.1:maximum(γ.lat)),
               rand(minimum(γ.depth):1.0:maximum(γ.depth)))

        iswet(loc,γ) && return loc
        println("dry point, try again")
    end # if not, then start over.
end

"""
    function iswet(loc, γ, neighbors)
    Get (lon,lat,depth) tuples of wet locations.
    Allow a location to be wet if at least one out of 8 nearby gridpoints is wet.
    Certainly "wet" gridpoints could be defined more strictly.
# Arguments
- `loc`: lon,lat,depth
- `γ`
- `neighbors`
# Output
- Boolean (true for ocean point)
"""
function iswet(loc,γ,neighbors)
    # two approaches
    # approach 2
    # find 8 nearest neighbors
    Inn = nearestneighbor(loc,γ,neighbors)

    # are any of them wet?
    for ii = 1:neighbors
        if γ.wet[Inn[ii]]
            return true
        end
    end
    return false
end

function iswet(loc,γ)

    # wetness bounded by 0 and 1
    # should be an argument
    # 1 = very strict
    # 0 = all points
    wetness = 0.2
    
    wis = interpindex(loc,γ)

    # handle wraparound
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)

    # are any of them wet?
    # interpolate ones and zeros on to this loc.
    # if there is land nearby, the interpolated value
    # will be greater than 0.
    # this criterion only requires on land point nearby,
    # where nearby is one of the 8 corners of the cube that contains loc
    return Interpolations.InterpGetindex(wetwrap)[wis...] > wetness
end

wet(a::BoundaryCondition) = a.wet
wet(a::Field) = a.γ.wet
wet(a::Source) = a.γ.interior

function _read3d(file,tracername)
    ds = Dataset(file,"r")
    v = ds[tracername]

    # eliminate Union{Missing} types
    T = eltype(v[1,1,1])
    c = convert(Array{T,3},v[:,:,:])
    
    # load an attribute
    if "units" in keys(v.attrib)
        units = v.attrib["units"]
    else
        error("TMI._read3d: units not found")
    end

    if "longname" in keys(v.attrib)
        longname = v.attrib["longname"]
    elseif "long_name" in keys(v.attrib)
        longname = v.attrib["long_name"]
    else
        error("TMI._read3d: longname not found")
    end
        
    close(ds)
    return c, units, longname
end

end
