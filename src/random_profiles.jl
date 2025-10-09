@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)

"""
    function wetlocation(γ)
    Get (lon,lat,depth) tuples of wet locations.
    Allow a location to be wet if at least one out of 8 nearby gridpoints is wet.
    Certainly "wet" gridpoints could be defined more strictly.
# Arguments
- `γ`: TMI.grid
# Output
- `loc`: lon,lat,depth """
function wetsurfacelocation(γ; sampling=:uniform)
    #Paban Modified
        confirmwet = false
        neighbors  = 8
        while !confirmwet
            if sampling == :uniform
                loc = (rand(minimum(γ.lon):0.1:maximum(γ.lon)),
                    rand(minimum(γ.lat):0.1:maximum(γ.lat)),
                    γ.depth[1])
                iswet(loc,γ) && return loc[1:2]
                println("dry point, try again")
            elseif sampling == :spherical
                #from the following source: 
                #https://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
                randlat = acos(2rand()-1) - 90 #random number -90, 90
                randlon = 2π * rand()
                loc = (randlon, randlat, γ.depth[1])
                iswet(loc,γ) && return loc[1:2]
                println("dry point, try again")
            else
                throw("specified method not implemented")

            end
        end # if not, then start over.
end

function wetsurfacelocation(γ1, γ2; sampling_method=:uniform)
    #Paban Modified
        confirmwet = false
        neighbors  = 8
        min_lon = minimum([minimum(γ1.lon), minimum(γ2.lon)])
        max_lon = maximum([maximum(γ1.lon), maximum(γ2.lon)])
        min_lat= minimum([minimum(γ1.lat), minimum(γ2.lat)])
        max_lat = maximum([maximum(γ1.lat), maximum(γ2.lat)])

        while !confirmwet
            if sampling_method == :uniform
                loc = (rand(min_lon:0.1:max_lon),
                    rand(min_lat:0.1:max_lat),
                    γ1.depth[1]) #use the first grid's depth grid, probably not the best but a quick fixx
                (iswet(loc,γ1) * iswet(loc,γ2)) && return loc[1:2]
                println("dry point, try again")
            elseif sampling_method == :spherical
                #from the following source: 
                #https://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
                randlat = rad2deg(acos(2rand()-1)) - 90 #random number -90, 90
                randlon = rad2deg(2π * rand())
                # randlat = acos((2 * rand(-1:0.0005:1))-1) - 90 #random number -90, 90
                # randlon = 2π * rand(-1:0.0005:1)

                loc = (randlon, randlat, γ1.depth[1])
                (iswet(loc,γ1) * iswet(loc,γ2)) && return loc[1:2]
                println("dry point, try again")
            else
                throw("specified method not implemented")
            end
        end # if not, then start over.
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
- `y`: contaminated profile of observations on 3D grid
- `W⁻`: appropriate weighting (inverse covariance) matrix for these observations,
- `ytrue`: uncontaminated observations, 3D field
- `locs`: 3-tuples of locations for observations
- `wis`: weighted indices for interpolation to locs sites
"""
function random_profiles(TMIversion,variable,γ,N; σ=0.0, locs = nothing)

    TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")
    
    θtrue = readfield(TMIfile,variable,γ)
    replace!(θtrue.tracer,NaN=>0.0)

    # get random locations that are wet (ocean)
    if isnothing(locs)
        locs = Vector{Tuple{Float64,Float64}}(undef,N)
        [locs[i] = wetsurfacelocation(γ) for i in eachindex(locs)]
    end

    # get weighted interpolation indices
    Nl = length(locs)
    Nz = length(γ.depth)

    yp = Array{Float32}(undef,(Nz, Nl))
    yptrue = Array{Float32}(undef,(Nz, Nl))

    for (il, loc) in enumerate(locs)
        wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,Nz)
        
        [wis[id] = interpindex(tuplejoin(loc, γ.depth[id]),γ) for id in 1:Nz]

        ytrue = observe(θtrue,wis,γ)
        σtrue = σ * ones(Nz)

        ntrue = rand(Normal(),Nz).*σtrue

        yp[:, il] .= ytrue .+ ntrue
        yptrue[:, il] .= ytrue

    end

    return yp, θtrue, yptrue, locs
end

function random_profiles(θtruetracer::AbstractArray{T, 3},γ,N; σ=0.0, locs = nothing) where T <: Real
    
    θtrue = copy(θtruetracer)
    replace!(θtrue,NaN=>0.0)

    # get random locations that are wet (ocean)
    if isnothing(locs)
        locs = Vector{Tuple{Float64,Float64}}(undef,N)
        [locs[i] = wetsurfacelocation(γ) for i in eachindex(locs)]
    else
        println("using predfined locations")
    end

    # get weighted interpolation indices
    Nl = length(locs)
    Nz = length(γ.depth)

    yp = Array{Float32}(undef,(Nz, Nl))
    yptrue = Array{Float32}(undef,(Nz, Nl))

    for (il, loc) in enumerate(locs)
        wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,Nz)
        
        [wis[id] = interpindex(tuplejoin(loc, γ.depth[id]),γ) for id in 1:Nz]

        ytrue = observe(θtrue,wis,γ)
        σtrue = σ * ones(Nz)

        ntrue = rand(Normal(),Nz).*σtrue

        yp[:, il] .= ytrue .+ ntrue
        yptrue[:, il] .= ytrue

    end

    return yp, θtrue, yptrue, locs
end


"""
    function observe
    Take a observation at location given by weights wis
"""
function observe(ctracer::AbstractArray{T, 3},wis::Vector{Tuple{Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}, Interpolations.WeightedAdjIndex{2,T}}},γ::Grid)::Vector{T} where T <: Real

        # look at total weight, < 1 if there are land points
    # later make sure total weight = 1 for proper average
    sumwis = Vector{Float64}(undef,length(wis))
    list = vcat(1:length(γ.lon),1)
    wetwrap = view(γ.wet,list,:,:)
    [sumwis[i] = Interpolations.InterpGetindex(wetwrap)[wis[i]...] for i in eachindex(wis)]

    # sample the true field at these random locations
    y = Vector{Float64}(undef,length(wis))
    replace!(ctracer,NaN=>0.0)
    ywrap = view(ctracer,list,:,:)
    [y[i] = Interpolations.InterpGetindex(ywrap)[wis[i]...]/sumwis[i] for i in eachindex(wis)]

    return y
end