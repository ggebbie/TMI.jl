"""
Surface oxygen saturation value and fraction of saturation value in field 
"""
function surface_oxygensaturation(file)
    # read temperature and o2.
    θ = readtracer(file,"θ")
    θsurface = view(θ,:,:,1)

    S = readtracer(file,"Sp")
    Ssurface = view(S,:,:,1)

    # GibbsSeaWater.jl for saturation value
    O₂sol = gsw_o2sol_sp_pt.(Ssurface, θsurface)

    O₂ = readtracer(file,"O₂")
    O₂surface = view(O₂,:,:,1)
    O₂fraction = O₂surface./O₂sol

    return O₂sol, O₂fraction 
end

"""
Reconstruct dissolved oxygen (that doesn't exist in TMI product)
by assuming same oxygen saturation fraction as modern
"""
function oxygen(version,O₂fraction)

    A, Alu, γ, file = config_from_nc(version)

    o2po4ratio = 170
    
    # read temperature and o2.
    θ = readtracer(file,"θ")
    θsurface = view(θ,:,:,1)

    S = readtracer(file,"Sp")
    Ssurface = view(S,:,:,1)

    # GibbsSeaWater.jl for saturation value
    O₂sol = gsw_o2sol_sp_pt.(Ssurface, θsurface)

    O₂surface = O₂sol.*O₂fraction

    # invert with stoichiometric ratio
    O₂ = tracerinit(γ.wet)
    qPO₄ = readtracer(file,"qPO₄")
    d = o2po4ratio*qPO₄

    #d = qO₂lgm
    d[:,:,1] = O₂surface;

    O₂[γ.wet] =  Alu\d[γ.wet]

    return O₂
end

