
import Pkg; Pkg.activate("./scripts")
Pkg.instantiate()
using Revise
using PlutoUI
using TMI
using Test
using PythonCall
using CondaPkg
using PythonPlot
const cartopy = pyimport("cartopy")
const matplotlib = pyimport("matplotlib")

ccrs = cartopy.crs

import Pkg;
using Interpolations
using Statistics
using LinearAlgebra
plotsdir(x) = "/Users/anthonymeza/Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology/Documents/GitHub/TMI.jl/plots/" * x
TMIversion = "modern_90x45x33_GH10_GH12";
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

# first guess of change to surface boundary conditions
# how many randomly sampled observations?
N = 50;
iters=25; 

# take synthetic, noisy observations
y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N);
y0 = zero(y) .+ mean(y)

function estimate_global_mean(yp)

	u = (;surface = zerosurfaceboundary(γ))
	uvec = vec(u)
	# make a silly first guess for surface
	b = (;surface = mean(yp) * onesurfaceboundary(γ))

	# assume temperature known ± 5°C
	σb = 5.0
	Dg = gaussiandistancematrix(γ,σb,1000.0);
	Q⁻ = inv(cholesky(Dg));
					
	out, f, fg, fg! = TMI.sparsedatamap(Alu,b,u,yp,W⁻,wis,locs,Q⁻,γ; iterations=iters);

	# reconstruct by hand to double-check.
	ũ = unvec(u,out.minimizer);
	b̃ = adjustboundarycondition(b,ũ);

	c̃  = steadyinversion(Alu,b̃,γ);

	return mean(c̃)
end


Δθ = 0.05
p_list = [-Δθ, 0.00, Δθ]
string(p_list[1])
c̃̄_dict = Dict()
for (i, p) in enumerate(p_list)
	println(p)
	nelm = (p == 0) ? 1 : N
	c̃̄_dict[p] = zeros(nelm)
	
	for j in 1:nelm
		yp = 1 .* y0 
		yp[j] += p

		c̃̄_dict[p][j]  = estimate_global_mean(yp);
	end
end

W = (c̃̄_dict[p_list[end]] .- c̃̄_dict[p_list[1]]) ./ (p_list[end] - p_list[1]);

println("mean of observations: ", mean(y))
println("observations: ", y)

println("estimated weights: ", W)

θ̄ = mean(readfield(TMIfile, "θ", γ));
θ̄̃0 =  c̃̄_dict[0.0][1];
θ̄̃ = (W' * (y - y0)) + θ̄̃0;
θ̄̃f = estimate_global_mean(y)

fig, ax = subplots(figsize = (5, 5))
labels = ["f(ŷ)", "f(ŷ₀)", "W(ŷ - ŷ₀) + f(ŷ₀)", "True θ̄"]
values = [θ̄̃f, θ̄̃0, θ̄̃, θ̄]
bar_container = ax.bar(labels, values, color = ["r", "b", "purple", "black"])
ax.bar_label(bar_container, fmt="{:,.2f}", padding = 0.4)
ax.set_ylabel("deg C")
ax.set_title("Estimates of Global Mean Temperature (θ̄)\n using " * string(N) * " observations, h = " * string(Δθ) * " and\n" * string(iters) * " iterations " * "to estimate W")
fig
fig.savefig(plotsdir("estimates.png"), dpi = 250, bbox_inches = "tight")

figl, axl = subplots(figsize = (10, 10), 
						subplot_kw=Dict("projection"=> ccrs.PlateCarree()),)
# Convert data using log base 2
base = 1/N
log_data = log.(base, abs.(W))
# Normalize this manually (no need for LogNorm now)
norm = matplotlib.colors.LogNorm(vmin = 1/(10 * N), vmax=10 / N)
cmap = matplotlib.cm.get_cmap("bwr")

# norm = matplotlib.colors.Normalize(vmin = log(1/(100 * N)), vmax = log(100 * N))
ticks = [1/(10 * N), 1/(3 * N), 1/N, 3/N, 10/N]
# norm = matplotlib.colors.BoundaryNorm(ticks, cmap.N)

cmll = []
for i in 1:N
	lon, lat = locs[i][1:2]
	println(lon, " ", lat)
	cml = axl.scatter(lon, lat, c = abs.(W)[i], cmap = "bwr", 
					  edgecolors="black", transform = ccrs.PlateCarree(),norm=norm)
	push!(cmll, cml)
	axl.coastlines(facecolor="grey")
end
cb = figl.colorbar(cmll[1], orientation = "horizontal", 
ticks=ticks, extend = "both")
cb.set_ticks(ticks)  # Force the ticks to appear
cb.set_ticklabels(["1/10N", "1/3N", "1/N", "3/N", "10/N"])  # Set the custom tick labels
figl
figl.savefig(plotsdir("planview.png"), dpi = 250, bbox_inches = "tight")

TMI.regionlist()
regions1 = ["TROP", "N", "SUBANT"] .* "PAC"
regions2 = ["TROP", "SUBANT"] .* "IND"
regions = vcat(regions1, regions2)
mask = sum(TMI.surfaceregion(TMIversion,region).tracer for region in regions)
mask = mask .* γ.wet
wbasin = W[iswetmask(locs,mask)]
locsbasin = locs[iswetmask(locs,mask)]

figl, axl = subplots(figsize = (7.5, 4))
axl.set_title("Indo-Pacific")
Nbasin = length(wbasin)
cmll = []
for i in 1:Nbasin
	lon, lat, depth = locsbasin[i]
	println(lon, " ", lat)
	println( wbasin[i])

	cml = axl.scatter(lat, -depth, c = wbasin[i], cmap = "bwr", 
	vmin = minimum(W), vmax = maximum(W), edgecolor = "k")
	push!(cmll, cml)
end
cb = figl.colorbar(cmll[1], orientation = "horizontal")
# cb.set_ticklabels(["1/10N", "1/3N", "1/N", "3N", "10N"])
figl.savefig(plotsdir("IndoPacific.png"), dpi = 250, bbox_inches = "tight")
figl

regions = ["TROP", "N", "SUBANT"] .* "ATL"
mask = sum(TMI.surfaceregion(TMIversion,region).tracer for region in regions)
mask = mask .* γ.wet
wbasin = W[iswetmask(locs,mask)]
locsbasin = locs[iswetmask(locs,mask)]

figl, axl = subplots(figsize = (7.5, 4))
axl.set_title("Atlantic")
Nbasin = length(wbasin)
cmll = []
for i in 1:Nbasin
	lon, lat, depth = locsbasin[i]
	println(lon, " ", lat)
	println( wbasin[i])
	cml = axl.scatter(lat, -depth, c = wbasin[i], cmap = "bwr", 
	vmin = minimum(W), vmax = maximum(W), edgecolor = "k")
	push!(cmll, cml)
end
cb = figl.colorbar(cmll[1], orientation = "horizontal")
# cb.set_ticklabels(["1/10N", "1/3N", "1/N", "3N", "10N"])
figl.savefig(plotsdir("Atlantic.png"), dpi = 250, bbox_inches = "tight")
figl

mask = 1 .* γ.wet
wbasin = W[iswetmask(locs,mask)]
locsbasin = locs[iswetmask(locs,mask)]
figl, axl = subplots(figsize = (7.5, 4))
axl.set_title("Global")
Nbasin = length(wbasin)
cmll = []
for i in 1:Nbasin
	lon, lat, depth = locsbasin[i]
	println(lon, " ", lat)
	println( wbasin[i])
	cml = axl.scatter(lat, -depth, c = wbasin[i], cmap = "bwr", 
	vmin = minimum(W), vmax = maximum(W), edgecolor = "k")
	push!(cmll, cml)
end
cb = figl.colorbar(cmll[1], orientation = "horizontal")
figl.savefig(plotsdir("Global.png"), dpi = 250, bbox_inches = "tight")
figl
