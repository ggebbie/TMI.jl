#=
 Example : reconstruct glacial oxygen
1. diagnose sat o2 percentage modern ocean
2. using sat o2 percentage and glacial T for glacial o2 boundary condition
3. propagate into interior with stoichiometric ratio
4. difference with modern
5. plot it
=#
using Revise
using TMI, BenchmarkTools, PyPlot, PyCall, GibbsSeaWater
pygui(true)
modernversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, modernfile, L, B = config_from_nc(modernversion)

O₂solmodern,O₂fractionmodern = surface_oxygensaturation(modernfile)

# modern section
O₂modern = readtracer(modernfile,"O₂")
Osection = section(O₂modern,lon_section,γ)
titlelabel = "O₂ Modern "*string(lon_section)*"°E"
# make a plot of dye in the ocean
dyeplot(γ.lat,γ.depth,Osection,lims,titlelabel)
figname = plotsdir("oxygen_30W_modern.eps")
savefig(figname)

# next get LGM temperature
lgmnames = ("G14","G14A","GPLS1","GPLS2","OG18")
for lgmname in lgmnames
    
    LGMversion = "LGM_90x45x33_"*lgmname
    O₂lgm = oxygen(LGMversion,O₂fractionmodern)

    # plot a section at 330 east longitude (i.e., 30 west)
    lon_section = 330;
    lims = 50:10:350

    # LGM section
    Osection = section(O₂lgm,lon_section,γ)
    titlelabel = "O₂ "*lgmname*string(lon_section)*"°E"
    # make a plot of dye in the ocean
    dyeplot(γ.lat,γ.depth,Osection,lims,title)
    figname = plotsdir("oxygen_30W_"*lgmname*".eps")
    savefig(figname)

    # LGM - modern
    lims = -100:10:100
    Osection = section(O₂lgm - O₂modern,lon_section,γ)
    titlelabel = "O₂ "*lgmname*"-modern "*string(lon_section)*"°E"
    # make a plot of dye in the ocean
    dyeplot(γ.lat,γ.depth,Osection,lims,titlelabel)
    figname = plotsdir("oxygen_30W_"*lgmname*"-Modern.png")
    savefig(figname)

end
