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
	using GeoPythonPlot
	import Pkg;
	using Interpolations
	using Statistics
	using LinearAlgebra

	TMIversion = "modern_90x45x33_GH10_GH12";
	A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

	# first guess of change to surface boundary conditions
	# how many randomly sampled observations?
	N = 5;
	
	# take synthetic, noisy observations
	y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N);
	
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

# ╔═╡ a58b9af6-a6f3-46ab-83c6-023b48c2c35d
begin
	Cmin=0.0; Cmax=0.1;
	Cslider = @bind Δθ Slider(Cmin:0.01:Cmax, default=0.05);
	md"$(Cmin) h-value (deg C) $(Cslider) $(Cmax)";
	
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
	
	    u = (;surface = zerosurfaceboundary(γ))
	    uvec = vec(u)
	    for j in 1:nelm
	        yp = 1 .* y 
	        yp[j] += p
	        # make a silly first guess for surface
	        b = (;surface = mean(yp) * onesurfaceboundary(γ))
	
	        # assume temperature known ± 5°C
	        σb = 5.0
	        Dg = gaussiandistancematrix(γ,σb,1000.0);
	        Q⁻ = inv(cholesky(Dg));
	                        
	        out, f, fg, fg! = TMI.sparsedatamap(Alu,b,u,yp,W⁻,wis,locs,Q⁻,γ);
	
	        # reconstruct by hand to double-check.
	        ũ = unvec(u,out.minimizer);
	        b̃ = adjustboundarycondition(b,ũ);
	
	        c̃  = steadyinversion(Alu,b̃,γ);
	        # reconstruct tracer map
	        c̃̄_dict[p][j]  = mean(c̃);
	    end
	end
	W = (c̃̄_dict[p_list[end]] .- c̃̄_dict[p_list[1]]) ./ (p_list[end] - p_list[1]);
end

# ╔═╡ f008bf51-ab2f-4d39-beb1-9600dd3f2e00
begin
	md"selected h-value: $(Δθ)";
end

# ╔═╡ 8e0ac86e-ee05-476f-a848-be1fbeef8356
W

# ╔═╡ aee9f81c-bf85-4143-9a27-cdd22e0fbf8e
y

# ╔═╡ b8abae4f-3de7-40d0-b54a-69b4ec088d0a
begin 
	θ̄ = mean(readfield(TMIfile, "θ", γ));
	θ̄̃ = W' * y;
	θ̄̃f =  c̃̄_dict[0.0][1];

	fig, ax = GeoPythonPlot.subplots(figsize = (5, 5))

	labels = ["f(ŷ)", "Wŷ", "True θ̄"]
	values = [θ̄̃f, θ̄̃, θ̄]
	bar_container = ax.bar(labels, values)
	ax.bar_label(bar_container, fmt="{:,.2f}", padding = 0.4)
	ax.set_ylabel("deg C")
	ax.set_title("Estimates of Global Mean Temperature (θ̄)\n using " * string(N) * "observations and h = " * string(Δθ) * " to estimate W" )
	fig
end

# ╔═╡ b17873b3-c160-4f60-be25-71b372a4b3e3
ax.bar_label

# ╔═╡ 857cdc83-4519-4fee-a7c8-ddf8eb8c9da8
θ̄̃

# ╔═╡ 623175a9-9ccf-410c-95da-753cce721e87
θ̄̃f + θ̄̃

# ╔═╡ Cell order:
# ╠═998a1ac6-6af9-11ef-2107-77e2241a6863
# ╟─0b83d710-3cb7-47ef-ad2e-0d76391e8404
# ╠═7040376c-5aeb-4392-87de-90943b60e8f5
# ╠═a58b9af6-a6f3-46ab-83c6-023b48c2c35d
# ╠═f008bf51-ab2f-4d39-beb1-9600dd3f2e00
# ╠═8e0ac86e-ee05-476f-a848-be1fbeef8356
# ╠═aee9f81c-bf85-4143-9a27-cdd22e0fbf8e
# ╠═b8abae4f-3de7-40d0-b54a-69b4ec088d0a
# ╠═b17873b3-c160-4f60-be25-71b372a4b3e3
# ╠═857cdc83-4519-4fee-a7c8-ddf8eb8c9da8
# ╠═623175a9-9ccf-410c-95da-753cce721e87
