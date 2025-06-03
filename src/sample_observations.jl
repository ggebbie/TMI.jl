export estimate_global_mean
export estimate_global_mean_from_surface

function estimate_global_mean(Alu, W⁻, γ, locs, wis, yp; iters = 25)

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

function estimate_global_mean_from_surface(Alu, TMIversion; 
										   N = 20, iters = 50, Δθ, maxdepth = 50, σ = nothing)

	y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N, maxdepth, σ);
	y0 = zero(y) .+ mean(y)

    p_list = [-Δθ, 0.00, Δθ]
	string(p_list[1])
	c̃̄_dict = Dict()
	for (i, p) in enumerate(p_list)
		nelm = (p == 0) ? 1 : N
		c̃̄_dict[p] = zeros(nelm)		
		for j in 1:nelm
			yp = 1 .* y0 
			yp[j] += p
			c̃̄_dict[p][j]  = estimate_global_mean(Alu, W⁻, γ, locs, wis, yp; iters = iters);
		end
	end

	W = (c̃̄_dict[p_list[end]] .- c̃̄_dict[p_list[1]]) ./ (p_list[end] - p_list[1]);

	θ̄ = mean(readfield(TMIfile, "θ", γ)); #true mean

	θ̄̃0 =  c̃̄_dict[0.0][1]; #initial guess (just mean of observations)
	θ̄̃ = (W' * (y - y0)) + θ̄̃0; #posterior estimate using finite purturbations
	θ̄̃f = estimate_global_mean(Alu, W⁻, γ, locs, wis, y; iters = iters) #posterior estimate from actual function 

	labels = ["f(ŷ)", "f(ŷ₀)", "W(ŷ - ŷ₀) + f(ŷ₀)", "True θ̄"]
	values = [θ̄̃f, θ̄̃0, θ̄̃, θ̄]

	est_dict = Dict( labels .=> values)
	sample_dict = Dict( "y" => y, "locs" =>locs, "W" => W)
	println(sample_dict)
	return (est_dict, sample_dict)
end
