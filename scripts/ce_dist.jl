#= 
Following ex2,
 (1) make a watermassdistribution plot for specified longitude
 (2) at specified lat/lon values find % of GIN and LAB waters
 (3) plot %GIN, %LAB w.r.t. depth for each lat/lon position
=#
using Revise
using TMI, BenchmarkTools, PyPlot, PyCall

pygui(true) #needed for Atom, not sure what it will do in other places

TMIversion = "modern_180x90x33_GH11_GH12"
#TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# choose water mass (i.e., surface patch) of interest
# Enhancement: provide list of choices
list = ("GLOBAL","ANT","SUBANT",
            "NATL","NPAC","TROP","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

region = "GIN" 

# do numerical analysis
g = watermassdistribution(TMIversion,Alu,region,γ);

# plot a section at 330 east longitude (i.e., 30 west)
#lon_section = 329;
lon_section = 337;
gsection = section(g, lon_section, γ)
lims = 0:5:100

# make a plot of dye in the ocean @ lon_section
dyeplot(γ.lat,γ.depth[33:-1:1],100 * gsection[:,33:-1:1], lims)

#comparison w/ depth - janky! 
aind = findall(x->x ∈ (61,63,59), γ.lat)
oind = findall(x->x ∈ (339,337,335), γ.lon)   
gin = 100*watermassdistribution(TMIversion, Alu, "GIN", γ)
lab = 100*watermassdistribution(TMIversion, Alu, "LAB", γ)
pairs = hcat(aind, oind)

#make plot comparing
figure(figsize = (2, 10))
for i in 1:size(pairs)[1]
    subplot(size(pairs)[1],1,i)
    plot(γ.depth, gin[pairs[i,2], pairs[i,1], :])
    plot(γ.depth, lab[pairs[i,2], pairs[i,1], :])
    legend(["GIN", "LAB"])
    if i == 3
        xlabel("Depth [m]")
    end
    
    ylabel("%")
    title(string((γ.lat[pairs[i,1]], γ.lon[pairs[i,2]])))
end
