using Test
using TMI

@testset "TMI.jl" begin

    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH10_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

    @testset "steadyinversion" begin

        yPO₄ = readfield(TMIfile,"PO₄",γ)
        bPO₄ = getsurfaceboundary(yPO₄)
        PO₄pre = steadyinversion(Alu,bPO₄,γ)
        qPO₄ = TMI.readsource(TMIfile,"qPO₄",γ)
        b₀ = zerosurfaceboundary(γ)
        PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)
        PO₄total = PO₄ᴿ + PO₄pre
        PO₄direct = steadyinversion(Alu,bPO₄,γ,q=qPO₄)

        ## how big is the maximum difference?
        # could replace with abs function
        @test maximum(PO₄direct - PO₄total) < 0.1
        @test minimum(PO₄direct - PO₄total) > -0.1

        ## oxygen distribution, just be sure it runs
        yO₂ = readfield(TMIfile,"O₂",γ)
        bO₂ = getsurfaceboundary(yO₂)
        O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)

    end

    @testset "steadyinversion with Source type" begin

        yPO₄ = readfield(TMIfile,"PO₄",γ)
        bPO₄ = getsurfaceboundary(yPO₄)
        PO₄pre = steadyinversion(Alu,bPO₄,γ)
        qPO₄ = readfield(TMIfile,"qPO₄",γ)
        b₀ = zerosurfaceboundary(γ)
        PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)
        PO₄total = PO₄ᴿ + PO₄pre
        PO₄direct = steadyinversion(Alu,bPO₄,γ,q=qPO₄)

        ## how big is the maximum difference?
        # could replace with abs function
        @test maximum(PO₄direct - PO₄total) < 0.1
        @test minimum(PO₄direct - PO₄total) > -0.1

        ## oxygen distribution, just be sure it runs
        yO₂ = readfield(TMIfile,"O₂",γ)
        bO₂ = getsurfaceboundary(yO₂)
        O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)

    end

    ############################
    ## trackpathways
    @testset "trackpathways" begin

        @test isapprox(maximum(sum(A,dims=2)),1.0)
        @test minimum(sum(A,dims=2))> -1e-14
        
        latbox = [50,60]; # 50 N -> 60 N, for example.
        lonbox = [-50, 0]; # 50 W -> prime meridian

        c = trackpathways(Alu,latbox,lonbox,γ)

        @test maximum(c) ≤ 1.0
        @test minimum(c) ≥ 0.0
    end

    @testset "watermassdistribution" begin
        list = ("GLOBAL","ANT","SUBANT",
                "NATL","NPAC","TROP","ARC",
                "MED","ROSS","WED","LAB","GIN",
                "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
                "TROPATL","TROPPAC","TROPIND")
        region = list[2]
        g = watermassdistribution(TMIversion,Alu,region,γ)
        @test maximum(g) ≤ 1.0
        @test minimum(g) ≥ 0.0

    end
    
    #######################################
    ## example: regeneration
    @testset "regeneration" begin

        PO₄ᴿ = regeneratedphosphate(TMIversion,Alu,γ)
        @test maximum(PO₄ᴿ) < 5
        @test minimum(PO₄ᴿ) ≥ 0
    end

    #######################################
    ## example: mean age
    @testset "meanage" begin

        a = meanage(TMIversion,Alu,γ)
        @test maximum(a) < 3000
        @test minimum(a) ≥ 0
    end

    ########################################
    ## example: howoceansfilled
    @testset "howoceansfilled" begin

        # are the areas and volumes consistent?
        v = cellvolume(γ)
        area = cellarea(γ)
        @test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
        volume = volumefilled(TMIversion,Alu,γ)
        # volumefill no smaller than smallest box?
        @test exp10(minimum(volume)) ≥ 5.0
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
            ctmp = zeros(γ.wet)
            δ = interpweights(loc,γ)
        else
            # alternate choice
            δ = nearestneighbormask(loc,γ)
        end
        @test isapprox(sum(filter(!isnan,δ)),1.0) 

        origin = surfaceorigin(loc, Alu, γ)
        
        #@test isapprox(sum(filter(!isnan,origin)),1)
        #@test isapprox(sum(filter(!isnan,origin)),1)
        #@test minimum(origin) ≥ -20
        @test maximum(origin) ≤ 0 # log10(1) = 0
    end
    
    #########################################
    ## formerly filterdata
    @testset "steadyclimatology" begin

        for ii = 1:2
            y, W⁻, ctrue = synthetic_observations(TMIversion,"θ",γ)
            if ii == 1
                println("NamedTuple type")
                u = (;surface = zerosurfaceboundary(γ))
                b = (;surface = getsurfaceboundary(y))
            else
                println("BoundaryCondition type")
                u = zerosurfaceboundary(γ)
                b = getsurfaceboundary(y)
            end
            
            uvec = vec(u)
            F,G = costfunction_gridded_obs(uvec,Alu,b,u,y,W⁻,γ)
            fg!(F,G,x) = costfunction_gridded_obs!(F,G,x,Alu,b,u,y,W⁻,γ)
            fg(x) = costfunction_gridded_obs(x,Alu,b,u,y,W⁻,γ)
            f(x) = fg(x)[1]
            J₀,gJ₀ = fg(uvec)

            iterations = 10
            out = steadyclimatology(uvec,fg!,iterations)

            # check with forward differences
            ϵ = 1e-3
            ii = rand(1:sum(γ.wet[:,:,1]))
            println("Location for test =",ii)
            δu = copy(uvec); δu[ii] += ϵ
            ∇f_finite = (f(δu) - f(uvec))/ϵ
            println(∇f_finite)

            fg!(J₀,gJ₀,(uvec+δu)./2) # J̃₀ is not overwritten
            ∇f = gJ₀[ii]
            println(∇f)
            
            # error less than 10 percent?
            println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
            @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

            # was cost function decreased?
            @test out.minimum < J₀

            # reconstruct by hand to double-check.
            ũ = out.minimizer
            J,gJ = fg(ũ)
            @test J < J₀
        end
    end
    
    #########################################
    ## globalmap
    @testset "sparsedatamap" begin

        using Statistics, Interpolations
        
        N = 20
        y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N)

        for ii = 1:2
            if ii == 1
                println("NamedTuple type")
                u = (;surface = zerosurfaceboundary(γ))
                b = (;surface = mean(y) * onesurfaceboundary(γ))
            else
                println("BoundaryCondition type")
                u = zerosurfaceboundary(γ)
                b = mean(y) * onesurfaceboundary(γ)
            end
            uvec = vec(u)
            σb = 5.0
            Q⁻ = 1.0/(σb^2)
            fg(x) = costfunction_point_obs(x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ)
            f(x) = fg(x)[1]
            J0 = f(uvec)
            J̃₀,∂J₀∂u = fg(uvec)
            fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ)

            ϵ = 1e-3 # size of finite perturbation
            ii = rand(1:sum(γ.wet[:,:,1]))
            println("gradient check location=",ii)
            δu = copy(uvec); δu[ii] += ϵ
            ∇f_finite = (f(δu) - f(uvec))/ϵ
            println("∇f_finite=",∇f_finite)

            ∂J₀∂u = 0.0 .* uvec
            fg!(J̃₀,∂J₀∂u,(uvec+δu)./2) # J̃₀ is not overwritten
            ∇f = ∂J₀∂u[ii]
            println("∇f=",∇f)

            # error less than 10 percent?
            println("Percent error=",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
            @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1
            iters = 5
            out = sparsedatamap(uvec,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,iterations=iters)

            # was cost function decreased?
            @test out.minimum < J̃₀
            # reconstruct by hand to double-check.
            ũ = out.minimizer
            J̃,∂J̃∂ũ = fg(ũ)
            @test J̃ < J̃₀
        end
    end
    @testset "sourcemap" begin

        using Statistics, Interpolations

        yPO₄ = readfield(TMIfile,"PO₄",γ)
        bPO₄ = getsurfaceboundary(yPO₄)

        #PO₄pre = steadyinversion(Alu,bPO₄,γ)
        qPO₄ = readfield(TMIfile,"qPO₄",γ)
        #b₀ = bPO₄ #zerosurfaceboundary(γ)
        #PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)

        N = 20
        y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"PO₄",γ,N)

        #u = (;surface = zerosurfaceboundary(γ))
        u = (; source = zeros(γ))
        #b = (; surface = bPO₄)
        b = (; surface = bPO₄) # assume a surface boundary condition

        PO₄ = steadyinversion(Alu,b,γ,q=qPO₄)
        uvec = vec(u)

        σq = 0.1
        Q⁻ = 1.0/(σq^2) # how well is q (source) known?
        fg(x) = costfunction_point_obs(x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=qPO₄)
        f(x) = fg(x)[1]
        J0 = f(uvec)
        J̃₀,∂J₀∂u = fg(uvec)
        fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q₀=qPO₄)

        ϵ = 1e-3 # size of finite perturbation
        ii = rand(1:sum(γ.wet[:,:,1]))
        println("gradient check location=",ii)
        δu = copy(uvec); δu[ii] += ϵ
        ∇f_finite = (f(δu) - f(uvec))/ϵ
        println("∇f_finite=",∇f_finite)

        fg!(J̃₀,∂J₀∂u,(uvec+δu)./2) # J̃₀ is not overwritten
        ∇f = ∂J₀∂u[ii]
        println("∇f=",∇f)

        # error less than 10 percent?
        println("Percent error=",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
        @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1
        iters = 5
        out = sparsedatamap(uvec,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=qPO₄,r=1.0,iterations=iters)
        # was cost function decreased?
        @test out.minimum < J̃₀
        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,∂J̃∂ũ = fg(ũ)
        @test J̃ < J̃₀
    end
    @testset "source_and_surfacemap" begin

        using Statistics, Interpolations

        yPO₄ = readfield(TMIfile,"PO₄",γ)
        bPO₄ = getsurfaceboundary(yPO₄)

        qPO₄ = readfield(TMIfile,"qPO₄",γ)

        N = 20
        y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"PO₄",γ,N)

        u = (; surface = zerosurfaceboundary(γ), source = zeros(γ))
        #b = (; surface = bPO₄)
        b = (; surface = bPO₄) # assume a surface boundary condition

        PO₄ = steadyinversion(Alu,b,γ,q=qPO₄)
        uvec = vec(u)

        σq = 0.1
        Q⁻ = 1.0/(σq^2) # how well is q (source) known?
        fg(x) = costfunction_point_obs(x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=qPO₄)
        f(x) = fg(x)[1]
        J0 = f(uvec)
        J̃₀,∂J₀∂u = fg(uvec)
        fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q₀=qPO₄)
        
        ϵ = 1e-3 # size of finite perturbation
        ii = rand(1:length(uvec))
        println("gradient check location=",ii)
        δu = copy(uvec); δu[ii] += ϵ
        ∇f_finite = (f(δu) - f(uvec))/ϵ
        println("∇f_finite=",∇f_finite)

        fg!(J̃₀,∂J₀∂u,(uvec+δu)./2) # J̃₀ is not overwritten
        ∇f = ∂J₀∂u[ii]
        println("∇f=",∇f)

        # error less than 10 percent?
        println("Percent error=",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
        @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1
        iters = 5
        out = sparsedatamap(uvec,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=qPO₄,r=1.0,iterations=iters)
        # was cost function decreased?
        @test out.minimum < J̃₀
        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,∂J̃∂ũ = fg(ũ)
        @test J̃ < J̃₀
    end

    @testset "unvec" begin
        u = (;surface = onesurfaceboundary(γ),
             west = 2 * onewestboundary(γ),
             east = 3 *oneeastboundary(γ))

        vecu = vec(u)
        
        u2 = (;surface = zerosurfaceboundary(γ),
             west = 2 * zerowestboundary(γ),
             east = 3 * zeroeastboundary(γ))

        u3 = unvec(u2,vecu)
        @test u.surface.tracer[10,10] == u3.surface.tracer[10,10]
        @test u.east.tracer[10,10] == u3.east.tracer[10,10]
        @test u.west.tracer[10,10] == u3.west.tracer[10,10]
        @test u.surface.tracer[10,10] == 1.0
        @test u.west.tracer[10,10] == 2.0
        @test u.east.tracer[10,10] == 3.0
        
        unvec!(u2,vecu)
        @test u.surface.tracer[10,10] == u2.surface.tracer[10,10]
        @test u.east.tracer[10,10] == u2.east.tracer[10,10]
        @test u.west.tracer[10,10] == u2.west.tracer[10,10]
        @test u.surface.tracer[10,10] == 1.0
        @test u.west.tracer[10,10] == 2.0
        @test u.east.tracer[10,10] == 3.0

    end
end

# @testset "regional" begin

#     shellscript = TMI.pkgsrcdir("read_nc_nordic_lowresolution.sh")
#     run(`sh $shellscript`)

# # NordicSeas_30-Dec-2020_low.mat
#     TMIversion = "NordicSeas_30-Dec-2020_low";
#     TMIfile = TMIversion*".nc"
#     !ispath(TMI.pkgdatadir()) && mkpath(TMI.pkgdatadir())
#     mv(joinpath(pwd(),TMIfile),TMI.pkgdatadir(TMIfile),force=true)

# end
