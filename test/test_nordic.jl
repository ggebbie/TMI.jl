    @testset "nordic" begin
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
