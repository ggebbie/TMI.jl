#=
% Example: Compare average δ¹³C and PO₄ values for different ocean basins %
=#

using Revise
using TMI, PyPlot, PyCall, Colors

# pygui(true) #needed for Atom, not sure what it will do in other places

# A, Alu, γ, inputfile = config(url,inputdir)
# ΔPO₄ = readtracer(inputfile,"qpo4")
TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
list = ("NATL","NPAC","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

# contains the percentage of volume sourced from each 
# surface ocean 
depth = getVar(TMIversion, "depth")

#initialize a grid containing the water mass distriubtion
#for each ocean in "list" 

which_water_mass = identifywatermass(TMIversion,Alu,γ, list)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330; # only works if exact
# wm_section = section(which_water_mass,lon_section,γ)

lims = 1:1.0:length(list)

cm = get_cmap(:tab20)
colorrange = (0:(length(list)-1)) ./ (length(list))

# oceans_mask = similar(wm_section[:,33:-1:1])

XC = γ.lon' .* ones(length(γ.lat))
YC = ones(length(γ.lon))' .*  γ.lat

fig, ax = subplots(1,1)
counter = []
for ff in eachindex(XC)
    val = which_water_mass[:, :, 20]'[ff]
    if isnan(val)
        # pass, dont plot anything
    else
        ax.scatter(XC[ff],YC[ff], color=cm(colorrange[Int(val)]))
        if Int(val) ∉ counter
            ax.scatter(XC[ff],YC[ff], color=cm(colorrange[Int(val)]), label = list[Int(val)])
            counter = push!(counter, Int(val))
        end
    end
#oceans_mask[:, :, 1]', lims, cmap=cmap_seismic
end
ax.legend()

#iterate through maskes for Ind, NPac, SPac, Na and Sa 
# possibley include equatorial waters 
# O13_values .* mask 
# PO_values .* mask 
# scatter(O13_values[:], PO_values[:])

basin_list = ["","SUBANTIND", "NPAC",
            "SUBANTPAC","NATL","SUBANTATL"]
            # "TROPATL","TROPPAC","TROPIND")

δ¹³C = getVar(TMIversion,"δ¹³C")
PO₄ = getVar(TMIversion,"PO₄")

fig, axs = plt.subplots(2, 3, figsize = (13, 8))
# lvl = 2 # right below surface
# lvl = 14
lvl = 19 
for i in 2:length(basin_list)
    ax = axs[i] 
    basin_mask = getRegion(TMIversion,"d_" * basin_list[i])
    og_pts = sum(filter(!isnan, basin_mask))
    println("Original mask has ",og_pts , " points")
    trop_basin_mask = getRegion(TMIversion,
                        "d_" * "TROP" * basin_list[i][end-2:end])
    
    if (basin_list[i][end-2:end] == "IND")
        basin_list[i] = "IND"
        hemisp = (γ.lat .> -1000)
    elseif basin_list[i][1] == 'N'
        hemisp = (γ.lat .> 0)
        basin_list[i] = "N" * basin_list[i][end-2:end]
    elseif basin_list[i][1] == 'S'
        hemisp = (γ.lat .< 0)
        basin_list[i] = "S" * basin_list[i][end-2:end] 
    end

    ax.set_title(basin_list[i])
    
    ax.set_xlim(0, 3.5)
    ax.set_xticks(collect(0:0.5:3.5))
    ax.set_ylim(-1, 3)
    ax.set_yticks(collect(-1:0.5:3))

    ax.set_xlabel("PO₄")
    ax.set_ylabel("δ¹³C")

    equatorpoints = findall(trop_basin_mask[:, hemisp] .== 1.0)
    hem_basin_mask = @view basin_mask[:, hemisp]
    hem_basin_mask[equatorpoints] .= 1.0
    println("Changed mask has ",sum(filter(!isnan, basin_mask)) - og_pts, " addtl points")

    mask_water_mass = which_water_mass[:, :, lvl:end] .* basin_mask
    good_pts = findall(x->(x!= NaN) && ( x >0.0),mask_water_mass)
    println("good points are: ", size(good_pts), "(should be geq to og + addtl pts)")
    mask_water_mass = mask_water_mass[good_pts]
    mask_δ¹³C = (δ¹³C[:, :, lvl:end] .* basin_mask)[good_pts]
    mask_PO₄ = (PO₄[:, :, lvl:end] .* basin_mask)[good_pts]
    counter = []

    for ff in eachindex(mask_water_mass)
        val = mask_water_mass[ff]
        ax.scatter(mask_PO₄[ff],mask_δ¹³C[ff], color=cm(colorrange[Int(val)]), s = 1)
        if Int(val) ∉ counter
            ax.scatter(mask_δ¹³C[ff], mask_δ¹³C[ff], color=cm(colorrange[Int(val)]), 
            label = list[Int(val)],s = 5)
            counter = push!(counter, Int(val))
        end
        
    end
    ax.legend()

    # good_idxs = findall( NaN .!= basin_mask .== 1.0)
    # findall(x->(x!= NaN)&&(x>1),x)
    # X = ones(length(δ¹³C[good_idxs]), 2)
    # X[:, 2] .= δ¹³C[good_idxs][:]
    # y = PO₄[good_idxs]
    # β̂ = (X'*X)\X'*y
    # ŷ=X*β̂
    # ax.plot(δ¹³C[good_idxs][:], ŷ, color = "black")
end
depth_range = string(depth[lvl]) * "m" * " to " * string(depth[end])
fig.suptitle(depth_range)
fig.tight_layout()
fig.savefig("PO₄vsδ¹³C_" * depth_range * ".png")
close("all")
