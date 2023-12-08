using Test
using TMI

function compare_controls(ubc,ubc2,testval)
    irand = rand(1:sum(ubc.wet))
    iloc = findall(ubc.wet)[irand]
    @test ubc.tracer[iloc] == ubc2.tracer[iloc]
    @test ubc.tracer[iloc] == testval
end

@testset "TMI.jl" begin

    @testset "regional" begin
        # error-prone and long = run it first
        TMIversion = "nordic_201x115x46_B23"
        A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

        @testset "sourcemap" begin

            using Statistics, Interpolations

            lscale = true
            surfacetoo = false
            
            yPO₄ = readfield(TMIfile,"PO₄",γ)
            bPO₄ = (; surface = getsurfaceboundary(yPO₄),
                    north = getnorthboundary(yPO₄),
                    east = geteastboundary(yPO₄),
                    south = getsouthboundary(yPO₄),
                    west = getwestboundary(yPO₄))
            qPO₄ = readsource(TMIfile,"qPO₄",γ) # truth
            q₀ = 1e-2*onesource(γ) # first guess
            
            N = 20
            σ = 0.01 # incomplete first guess error, use scalar
            y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"PO₄",γ,N,σ)

            # Adjust surface only
            #u = (; surface = zerosurfaceboundary(γ), source = zerosource(γ))

            # Adjust interior sources only
            u = (; source = zerosource(γ,logscale=lscale))

            # Adjust interior sources and lateral boundary conditions
            u = (; source = zerosource(γ,logscale=lscale),
                 north = zeronorthboundary(γ),
                 east = zeroeastboundary(γ),
                 south = zerosouthboundary(γ),
                 west = zerowestboundary(γ))

            PO₄true = steadyinversion(Alu,bPO₄,γ,q=qPO₄)
            PO₄₀ = steadyinversion(Alu,bPO₄,γ,q=q₀)
            uvec = vec(u)

            σq = 1.0
            Q⁻ = 1.0/(σq^2) # how well is q (source) known?
            fg(x) = costfunction_point_obs(x,Alu,bPO₄,u,y,W⁻,wis,locs,Q⁻,γ,q=q₀)
            f(x) = fg(x)[1]
            J0 = f(uvec)
            J̃₀,∂J₀∂u = fg(uvec)
            fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,bPO₄,u,y,W⁻,wis,locs,Q⁻,γ,q₀=q₀)

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
            out = sparsedatamap(uvec,Alu,bPO₄,u,y,W⁻,wis,locs,Q⁻,γ,q=q₀,r=1.0,iterations=iters)
            # was cost function decreased?
            @test out.minimum < J̃₀
            # reconstruct by hand to double-check.
            ũ = out.minimizer
            J̃,∂J̃∂ũ = fg(ũ)
            @test J̃ < J̃₀
            #@test J̃ < 3N # too strict if first guess is bad

            if lscale
                #b̃ = adjustboundarycondition(b₀,unvec(u,ũ)) # combine b₀, u
                q̃ = TMI.adjustsource(q₀,unvec(u,ũ))

                # next test likes to fail if ≥ 0
                # sources may permit negative values
                #@test minimum(q̃) ≥ -0.1 

                # get new boundary conditions
                b̃ = TMI.adjustboundarycondition(bPO₄,unvec(u,ũ))

                # are lateral boundaries adjusted?
                @test ~iszero(minimum(b̃.south-bPO₄.south))
            end
        end
    end

    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH11_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

    @testset "regions" begin
        # test that all global values sum to 1.
        Nr = length(TMI.regionlist())
        sumb = zeros(Int,Nr)
        for i in 1:Nr
            # read v1 of regions from NetCDF file: used Floating point numbers for mask
            sumb[i] = sum(TMI.surfaceregion(TMIversion,TMI.regionlist()[i]).tracer)
        end

        @test sum(sumb[1]) == sum(sumb[2:8])
        @test sum(sumb[1]) == sum(sumb[vcat(5,7:Nr)])

    end

    @testset "steadyinversion" begin

        yPO₄ = readfield(TMIfile,"PO₄",γ)
        bPO₄ = getsurfaceboundary(yPO₄)
        PO₄pre = steadyinversion(Alu,bPO₄,γ)
        qPO₄ = readsource(TMIfile,"qPO₄",γ)
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

        # this line changes
        qPO₄ = readsource(TMIfile,"qPO₄",γ)
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

        #b = surfacepatch(lonbox,latbox,γ)
        c = trackpathways(Alu,latbox,lonbox,γ)

        @test maximum(c) ≤ 1.0
        @test minimum(c) ≥ 0.0
    end

    @testset "watermassdistribution" begin
        list = TMI.regionlist()
        region = list[2]

        #b = TMI.surfaceregion(TMIversion,region,γ)
        #g = steadyinversion(Alu,b,γ)
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
        #@test sum(0. .< v./area .< 1000.)/length(γ.I) == 1
        volume = volumefilled(TMIversion,Alu,γ)
        # volumefill no smaller than smallest box?
        @test exp10(minimum(volume)) ≥ 4.9
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

        @test 0.99 < sum(exp10.(origin.tracer[origin.wet])) < 1.01
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
            ii = rand(1:length(uvec))
            println("Location for test =",ii)
            δu = copy(uvec); δu[ii] += ϵ
            ∇f_finite = (f(δu) - f(uvec))/ϵ

            fg!(J₀,gJ₀,(uvec+δu)./2) # J̃₀ is not overwritten
            ∇f = gJ₀[ii]
            
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

        using Statistics, Interpolations, LinearAlgebra
        
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

            for jj = 1:2
                if jj == 1
                    Q⁻ = 1.0/(σb^2) # a scalar
                    # spatially uniform first-guess expected size
                elseif jj==2
                    lengthscale = 1000.0
                    Dg = gaussiandistancematrix(γ,σb,lengthscale)
                    Q⁻ = inv(cholesky(Dg))
                end
            
                fg(x) = costfunction_point_obs(x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ)
                f(x) = fg(x)[1]
                J0 = f(uvec)
                J̃₀,∂J₀∂u = fg(uvec)
                fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ)

                ϵ = 1e-3 # size of finite perturbation
                ii = rand(1:length(uvec))
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
    end

    @testset "sourcemap" begin

        using Statistics, Interpolations

        for lscale in (false,true)
            for surfacetoo in (false,true)
                
                yPO₄ = readfield(TMIfile,"PO₄",γ)
                bPO₄ = getsurfaceboundary(yPO₄)
                qPO₄ = readsource(TMIfile,"qPO₄",γ)
                q₀ = 1e-2*onesource(γ)
                
                N = 20
                y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"PO₄",γ,N)

                #u = (; source = zerosource(γ))
                if surfacetoo
                    u = (; surface = zerosurfaceboundary(γ), source = zerosource(γ))
                else
                    u = (; source = zerosource(γ,logscale=lscale))
                end
                
                b = (; surface = bPO₄) # surface boundary condition

                PO₄true = steadyinversion(Alu,b,γ,q=qPO₄)
                PO₄₀ = steadyinversion(Alu,b,γ,q=q₀)
                uvec = vec(u)

                σq = 1.0
                Q⁻ = 1.0/(σq^2) # how well is q (source) known?
                fg(x) = costfunction_point_obs(x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q=q₀)
                f(x) = fg(x)[1]
                J0 = f(uvec)
                J̃₀,∂J₀∂u = fg(uvec)
                fg!(F,G,x) = costfunction_point_obs!(F,G,x,Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,q₀=q₀)

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
                #@test J̃ < 3N # too strict if first guess is bad

                if lscale
                    #b̃ = adjustboundarycondition(b₀,unvec(u,ũ)) # combine b₀, u
                    q̃ = TMI.adjustsource(q₀,unvec(u,ũ))
                    #@test minimum(q̃) ≥ 0 # likes to fail from time to time
                end
            end
        end
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
        compare_controls(u.surface,u3.surface,1.0)
        compare_controls(u.east,u3.east,3.0)
        compare_controls(u.west,u3.west,2.0)
        
        unvec!(u2,vecu)
        compare_controls(u.surface,u2.surface,1.0)
        compare_controls(u.east,u2.east,3.0)
        compare_controls(u.west,u2.west,2.0)

    end

end
