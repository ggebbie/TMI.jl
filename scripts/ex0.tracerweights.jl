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

# ╔═╡ 7077e9c0-6a5d-11ef-0afa-65c9d0ddf101
begin
	#=
	     Find the distribution of a tracer given:
	     (a) the pathways described by A or its LU decomposition Alu,
	     (b) first-guess boundary conditions and interior sources given by d₀,
	     (c) perturbations to the surface boundary condition u
	    that best fits observations, y,
	    according to the cost function,
	    J = (ỹ - Ec)ᵀ W⁻¹ (ỹ - Ec)
	    subject to Ac = d₀ + Γ u₀.                 
	    W⁻¹ is a (sparse) weighting matrix.
	    See Supplementary Section 2, Gebbie & Huybers 2011.
	# Arguments
	- `u`: control vector of surface tracer perturbations
	- `Alu`: LU decomposition of water-mass matrix A
	- `y`: observations on 3D grid
	- `d₀`: first guess of boundary conditions and interior sources
	- `W⁻`: weighting matrix best chosen as inverse error covariance matrix
	- `wet`: BitArray mask of ocean points
	=#
	#using Revise, TMI, Interpolations, Statistics
	
	import Pkg; Pkg.activate(".")
	using PlutoUI
	using Revise
	using TMI
	using Test
	using GeoPythonPlot
	using Interpolations
	using Statistics
	using LinearAlgebra
	
	TMIversion = "modern_90x45x33_GH10_GH12"
	A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);
end

# ╔═╡ 6464fa7c-eac8-4673-9099-be8868d3046b
begin 
	@bind σb Slider(0:0.1:10, default=5)
end

# ╔═╡ 71f42840-06cd-4906-9f35-8da36fed1f77
println("Uncertainty of spatial covariance: ", σb)

# ╔═╡ 5cef02fc-7bfe-47ac-9a22-21a8a9db6c1d
@bind l Slider(0:10:5000, default=1000)

# ╔═╡ 73d52282-15ff-4fd9-8d9a-156fd1ed7092
begin
	# assume temperature known ± 5°C
	Dg = gaussiandistancematrix(γ,σb,l)
	# Dg[Dg .≈ 0] .= 0
	Q⁻ = inv(cholesky(Dg));

	fig, ax = GeoPythonPlot.subplots(1, 2)
	cm = ax[0].contourf(Q⁻[1:500, 1:500])
	fig.colorbar(cm, ax = ax[0])
	
	cm = ax[1].contourf(Dg[1:500, 1:500])
	fig.colorbar(cm, ax = ax[1])
	fig
end

# ╔═╡ 816180f8-66df-4bd2-ab76-785760199904
Dg

# ╔═╡ b8984986-7d4e-4e3d-8032-b398327ffc6f
println("Decorrelation length scale: ", l, " meters")

# ╔═╡ 26769ada-75c9-419c-998f-a02afd097c0f
@bind N Slider(0:1:50, default=5)

# ╔═╡ bc3ebfc8-b442-457f-a32e-917789aaf3fa
begin
	y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N);
	figy, axy = GeoPythonPlot.subplots()
	cmy = axy.pcolormesh(W⁻)
	figy.colorbar(cmy, ax = axy)
	figy
end

# ╔═╡ 0910605d-cdfa-43be-9dae-cd4d49c5ab91
inv(W⁻)

# ╔═╡ 4756be05-042e-4d76-8c07-05217c2e4319
σθ = readfield(TMIfile,"σ"*"θ",γ)

# ╔═╡ a0b1462a-bb62-4b06-b62d-20c1444cb4b8
extrema(σθ.tracer)

# ╔═╡ Cell order:
# ╠═7077e9c0-6a5d-11ef-0afa-65c9d0ddf101
# ╠═73d52282-15ff-4fd9-8d9a-156fd1ed7092
# ╠═6464fa7c-eac8-4673-9099-be8868d3046b
# ╠═816180f8-66df-4bd2-ab76-785760199904
# ╠═71f42840-06cd-4906-9f35-8da36fed1f77
# ╠═5cef02fc-7bfe-47ac-9a22-21a8a9db6c1d
# ╠═b8984986-7d4e-4e3d-8032-b398327ffc6f
# ╠═bc3ebfc8-b442-457f-a32e-917789aaf3fa
# ╠═26769ada-75c9-419c-998f-a02afd097c0f
# ╠═0910605d-cdfa-43be-9dae-cd4d49c5ab91
# ╠═4756be05-042e-4d76-8c07-05217c2e4319
# ╠═a0b1462a-bb62-4b06-b62d-20c1444cb4b8
