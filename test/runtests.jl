using TMI, Test

@testset "TMI.jl" begin

    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH10_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

    @testset "steadyinversion" begin

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

        u = zerosurfaceboundary(γ)
        u₀ = u.tracer[u.wet]
        y, W⁻, ctrue = synthetic_observations(TMIversion,"θ",γ)
        b = getsurfaceboundary(y)

        # check gradients with respect to costfunction!
        fg(x) = costfunction_obs(x,Alu,b,y,W⁻,γ)
        f(x) = fg(x)[1]
        J̃₀,gJ₀ = fg(u₀)
        fg!(F,G,x) = costfunction_obs!(F,G,x,Alu,b,y,W⁻,γ)

        iterations = 10
        out = steadyclimatology(u₀,fg!,iterations)

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

        N = 20
        u = zerosurfaceboundary(γ)
        u₀ = u.tracer[u.wet]
        y, W⁻, ctrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N)
        b = mean(y) * onesurfaceboundary(γ)
        σb = 5.0
        Q⁻ = 1.0/(σb^2)

        # gradient check
        # check with forward differences
        fg(x) = costfunction(x,Alu,b,y,W⁻,wis,locs,Q⁻,γ)
        f(x) = fg(x)[1]
        J0 = f(u₀)
        J̃₀,gJ₀ = fg(u₀)
        fg!(F,G,x) = costfunction!(F,G,x,Alu,b,y,W⁻,wis,locs,Q⁻,γ)

        ϵ = 1e-3 # size of finite perturbation
        # Note: ϵ=1e-5 fails tests sometimes due to no finite difference at all
        # Problem with types or rounding or precision?
        
        ii = rand(1:sum(γ.wet[:,:,1]))
        println("gradient check location=",ii)
        δu = copy(u₀); δu[ii] += ϵ
        ∇f_finite = (f(δu) - f(u₀))/ϵ
        println("∇f_finite=",∇f_finite)

        fg!(J̃₀,gJ₀,(u₀+δu)./4) # J̃₀ is not overwritten
        ∇f = gJ₀[ii]
        println("∇f=",∇f)

        # error less than 10 percent?
        println("Percent error=",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
        @test abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

        iterations = 5
        # optimize the sparse data map with an Optim.jl method
        out = sparsedatamap(u₀,Alu,b,y,W⁻,wis,locs,Q⁻,γ,iterations)

        # was cost function decreased?
        @test out.minimum < J̃₀

        # reconstruct by hand to double-check.
        ũ = out.minimizer
        J̃,gJ̃ = fg(ũ)
        @test J̃ < J̃₀

    end

    @testset "transientsimulation" begin
        using Interpolations, NaNMath, DifferentialEquations, LinearAlgebra

        latbox = [50,60]
        lonbox = [-50,0]
        d = surfacepatch(lonbox, latbox, γ) 
        dsfc =  d[:,:,1][γ.wet[:,:,1]]

        #following make_initial_conditions.m
        c0 = B * dsfc 

        #Fixed euler timestep approximation
        c = c0
        Δt = 1e-3 #this becomes unstable if you go any lower
        T  = 0.1
        Nt = T/Δt
        for tt = 1:Nt
            # forward Euler timestep
            c += L*c*Δt
            println("Σ c = ",sum(c))
        end
        gain_euler = sum(c .- c0)
        
        #Solving differential equation for fixed case 
        u0 = c0
        du = similar(u0)
        f(du,u,p,t) = mul!(du, L, u) 
        tspan = (0.0,T)
        func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
        prob = ODEProblem(func, u0, tspan)
        println("Solving fixed ODE")
        @time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,calck=false)
        println("ODE solved")

        #put sol into time x lon x lat x depth 
        sol_array = zeros((length(sol.t), 90,45,33))
        [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

        #stability check
        stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
        @test stable
        println("fixed bc stable: ", stable)

        #gain check - tracer concentration should increase 
        gain_ode = NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])
        println("Gain = ", gain_ode)
        @test gain_ode ≥ 0.0

        #compare forward euler timestep approx and solved ODE results 
        gain_error = abs(gain_ode - gain_euler)/(abs(gain_ode) + abs(gain_euler))
        @test gain_error < 0.1
        
        println("Gain percent error ",200gain_error,"%")

        #varying case stability check
        tsfc = [0, T]
        Csfc = zeros((2, length(dsfc)))
        Csfc[1, :] .= 1
        τ = 1/12
        li = LinearInterpolation(tsfc, 1:length(tsfc))
        LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
        BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
        Cb = similar(Csfc[1,:])
        surface_ind = findall(x->x[3] ==1, γ.I)

        p = (Csfc,surface_ind,τ,L,B,li,LC,BF,Cb) #parameters
        f(du, u, p, t) = varying!(du, u, p, t)
        func = ODEFunction(f, jac_prototype=L)
        prob = ODEProblem(func, u0, tspan,p)
        println("Solving varying ODE")
        @time sol = solve(prob, QNDF(),abstol=1e-2,reltol=1e-2,saveat=tsfc)
        println("Varying ODE solved")
        
        #put sol into time x lon x lat x depth 
        sol_array = zeros((length(sol.t), 90,45,33))
        [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]
        
        #stability check
        stable = true ? NaNMath.maximum(sol_array) < 1.000001 && NaNMath.minimum(sol_array) > -0.000001 : false
        @test stable
        println("Varying case stable: ", stable)
       
    end

end
