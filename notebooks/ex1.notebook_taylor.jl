### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 998a1ac6-6af9-11ef-2107-77e2241a6863
begin
	import Pkg; Pkg.activate("../scripts")
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

	TMIversion = "modern_90x45x33_GH10_GH12";
	A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

	# first guess of change to surface boundary conditions
	# how many randomly sampled observations?
	N = 100;
	iters=25; 
	
	# take synthetic, noisy observations
	y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N);
	y0 = zero(y) .+ mean(y)
end

# ╔═╡ 8e0ac86e-ee05-476f-a848-be1fbeef8356
begin
	using Dates 
	timenow =  Dates.format(now(), "HH:MM")
	md"Most recent plot updated on: $(timenow)"
end

# ╔═╡ 0b83d710-3cb7-47ef-ad2e-0d76391e8404
md"Given a set of surface observations, a first-guess of their points, their spatial covariance and a known circulation we can reconstruct the interior values using a composition of `sparsedatamap` and `steadyinversion`. Using the predicted interior points, we can then construct an estimate of mean global temperatures. Mathematically, we can write all of these steps as

$$\tilde{\bar{\theta} } = f(\hat{y})$$ 

Here, $\hat y$ represents our surface observations. $f$ is the function that takes observations to an estimate of mean global temperature, $\tilde{\bar{\theta} }$.

We will make the assumption that $f$ is linear. So that the previous equation can be equivalantly expressed as

$$\tilde{\bar{\theta} } = W \hat y$$ 

Here, $W$ is unknown. However, we can recover $W$ by taking the derivative of $\tilde{\bar{\theta} }$ with respect to $\hat y$. I.e. 

$$\frac{ \partial \tilde{\bar{\theta} }}{\partial \hat y} = W^T$$

Each component of $\frac{ \partial \tilde{\bar{\theta} }}{\partial \hat y}$ can be approximated as: 

$$\left (\frac{ \partial \tilde{\bar{\theta} }}{\partial \hat y}\right )_i = W_i \approx \frac{f(\hat{y} + h) - f(\hat{y} - h)}{2h}$$
"

# ╔═╡ 279164cd-9a72-4b34-b7fa-c537e2829293
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

# ╔═╡ a58b9af6-a6f3-46ab-83c6-023b48c2c35d
begin
	Cmin=0.0; Cmax=0.1;
	Cslider = @bind Δθ Slider(Cmin:0.01:Cmax, default=0.05);
	md"""
	h-value: $(Δθ) deg C $br $(Cmin) $(Cslider) $(Cmax)
	""";
	
end

# ╔═╡ 7040376c-5aeb-4392-87de-90943b60e8f5
begin
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
	
end

# ╔═╡ b8abae4f-3de7-40d0-b54a-69b4ec088d0a
begin 
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
end

# ╔═╡ 02930649-90e8-4419-8cde-78f855c8b96a
locs

# ╔═╡ b17873b3-c160-4f60-be25-71b372a4b3e3
begin 
	figl, axl = subplots(figsize = (10, 10), 
						 subplot_kw=Dict("projection"=> ccrs.PlateCarree()))
	cmll = []
	for i in 1:length(locs)
		lon, lat = locs[i][1:2]
		println(lon, " ", lat)
		cml = axl.scatter(lon, lat, c = abs(W[i]), cmap = "bwr", norm = matplotlib.colors.LogNorm(vmin = 1/(10*N), 
			vmax = 10*N), edgecolors="black", transform = ccrs.PlateCarree())
		push!(cmll, cml)
		axl.coastlines(facecolor="grey")
	end
	figl.colorbar(cmll[1])
	figl
end

# ╔═╡ 857cdc83-4519-4fee-a7c8-ddf8eb8c9da8


# ╔═╡ 623175a9-9ccf-410c-95da-753cce721e87


# ╔═╡ Cell order:
# ╠═998a1ac6-6af9-11ef-2107-77e2241a6863
# ╟─0b83d710-3cb7-47ef-ad2e-0d76391e8404
# ╠═279164cd-9a72-4b34-b7fa-c537e2829293
# ╠═7040376c-5aeb-4392-87de-90943b60e8f5
# ╠═a58b9af6-a6f3-46ab-83c6-023b48c2c35d
# ╠═8e0ac86e-ee05-476f-a848-be1fbeef8356
# ╠═b8abae4f-3de7-40d0-b54a-69b4ec088d0a
# ╠═02930649-90e8-4419-8cde-78f855c8b96a
# ╠═b17873b3-c160-4f60-be25-71b372a4b3e3
# ╠═857cdc83-4519-4fee-a7c8-ddf8eb8c9da8
# ╠═623175a9-9ccf-410c-95da-753cce721e87
