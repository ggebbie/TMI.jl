using Revise, TMI, Test

@testset "TMI.jl" begin

    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH10_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    
    ############################
    ## trackpathways
    @testset "trackpathways" begin

        @test isapprox(maximum(sum(A,dims=2)),1.0)
        @test minimum(sum(A,dims=2))> -1e-14
        
        #- define the surface patch by the bounding latitude and longitude.
        latbox = [50. , 60.]; # 50 N -> 60 N, for example.

        # mutable due to wraparound: don't use an immutable tuple
        lonbox = [-50. , 0.]; # 50 W -> prime meridian

        c = trackpathways(Alu,latbox,lonbox,γ)

        @test maximum(c[γ.wet]) ≤ 1.0
        @test minimum(c[γ.wet]) ≥ 0.0
    end
    
    #######################################
    ## example: regeneration
    @testset "regeneration" begin

        PO₄ᴿ = regeneratedphosphate(TMIversion,Alu,γ)
        @test maximum(PO₄ᴿ[γ.wet]) < 10
        @test minimum(PO₄ᴿ[γ.wet]) ≥ 0
    end
    
    ########################################
    ## example: howoceansfilled
    @testset "howoceansfilled" begin

        # are the areas and volumes consistent?
        v = cellvolume(γ)
        area = cellarea(γ)
        @test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
        
        # effectively take inverse of transpose A matrix.
        dVdd = tracerinit(γ.wet); # pre-allocate c
        dVdd[γ.wet] = Alu'\v[γ.wet]

        volume = volumefilled(TMIversion,Alu,γ)

        # volumefill positive at surface?
        @test sum(volume .< 0) == 0
        @test minimum(volume) ≤ 0.0
    end
    
    ####################################
    ## example: surfaceorigin
    @testset "surfaceorigin" begin

        # choose linear or nearest neighbor interpolation
        linearinterp = true #true
        
        # randomized location.
        loc = wetlocation(γ)

        # choose linear or nearest neighbor interpolation

        if linearinterp
            # temporary tracer to initialize linear interp
            ctmp = tracerinit(γ.wet)
            δ = interpweights(loc,γ)
        else
            # alternate choice
            δ = nearestneighbormask(loc,γ)
        end

        @test isapprox(sum(filter(!isnan,δ)),1.0) 

        origin = surfaceorigin(loc, Alu, γ)
        
        @test isapprox(sum(filter(!isnan,origin)),1)
        @test isapprox(sum(filter(!isnan,origin)),1)
        @test minimum(filter(!isnan,origin)) ≥ 0
        @test maximum(filter(!isnan,origin)) ≤ 1
    end
    
    #########################################
    ## formerly filterdata
    @testset "steadyclimatology" begin

        # first guess of change to surface boundary conditions
        # ocean values are 0
        u₀ = zeros(Float64,sum(γ.wet[:,:,1]))

        # take synthetic, noisy observations
        y, W⁻, ctrue = sample_observations(TMIversion,"θ",γ)

        # a first guess: observed surface boundary conditions are perfect.
        # set surface boundary condition to the observations.
        # below surface = 0 % no internal sinks or sources.
        d₀ = tracerinit(γ.wet)
        d₀[:,:,1] = y[:,:,1]

        # check gradients in misfit_gridded_data!
        fg(x) = costfunction_obs(x,Alu,d₀,y,W⁻,γ)
        f(x) = fg(x)[1]
        J̃₀,gJ₀ = fg(u₀)
        fg!(F,G,x) = costfunction_obs!(F,G,x,Alu,d₀,y,W⁻,γ)
        #fg!(J̃₀,gJ₀,u₀)
        # filter the data with an Optim.jl method

        iterations = 5
        out = steadyclimatology(u₀,fg!,iterations)
        # out = steadyclimatology(u₀,Alu,y,d₀,W⁻,fg!,γ)

        # check with forward differences
        ϵ = 1e-3
        ii = rand(1:sum(γ.wet[:,:,1]))
        println("Location for test =",ii)
        δu = copy(u₀); δu[ii] += ϵ
        ∇f_finite = (f(δu) - f(u₀))/ϵ
        println(∇f_finite)

        fg!(J̃₀,gJ₀,(u₀+δu)./2) # J̃₀ is not overwritten
        ∇f = gJ₀[ii]
        println(∇f)
        
        # error less than 10 percent?
        println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
        @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

        # was cost function decreased?
        @test out.minimum < J̃₀

        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,gJ̃ = fg(ũ)
        @test J̃ < J̃₀

    end
    
    #########################################
    ## globalmap
    @testset "sparsedatamap" begin
        # first guess of change to surface boundary conditions
        # how many randomly sampled observations?
        N = 20

        # ocean values are 0
        u₀ = zeros(Float64,sum(γ.wet[:,:,1]))

        # take synthetic, noisy observations
        y, W⁻, ctrue, locs, wis = sample_observations(TMIversion,"θ",γ,N)

        # does this help optimization stay stable?
        #W⁻ *= 1.0/100.0
            
        # make a silly first guess for surface
        d₀ = tracerinit(γ.wet)
        [d₀[γ.I[ii]] = 15.0 for ii ∈ eachindex(γ.I) if γ.I[ii][3] == 1]

        # permit surface deviations on order of 5°C
        Q⁻ = 1.0/(5.0^2)
        #Q⁻ = 10.0
        
        # gradient check
        # check with forward differences
        fg(x) = costfunction(x,Alu,d₀,y,W⁻,wis,locs,Q⁻,γ)
        f(x) = fg(x)[1]
        J̃₀,gJ₀ = fg(u₀)
        fg!(F,G,x) = costfunction!(F,G,x,Alu,d₀,y,W⁻,wis,locs,Q⁻,γ)

        ϵ = 1e-3 # size of finite perturbation
        # Note: ϵ=1e-5 fails tests sometimes due to no finite difference at all
        # Problem with types or rounding or precision?
        
        ii = rand(1:sum(γ.wet[:,:,1]))
        ii = 1804
        println("gradient check location=",ii)
        δu = copy(u₀); δu[ii] += ϵ
        ∇f_finite = (f(δu) - f(u₀))/ϵ
        println("∇f_finite=",∇f_finite)

        fg!(J̃₀,gJ₀,(u₀+δu)./2) # J̃₀ is not overwritten
        ∇f = gJ₀[ii]
        println("∇f=",∇f)

        # error less than 10 percent?
        println("Percent error=",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
        @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

        iterations = 5
        # optimize the sparse data map with an Optim.jl method
        out = sparsedatamap(u₀,Alu,d₀,y,W⁻,wis,locs,Q⁻,γ,iterations)

        # was cost function decreased?
        @test out.minimum < J̃₀

        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,gJ̃ = fg(ũ)
        @test J̃ < J̃₀

    end
    
end
