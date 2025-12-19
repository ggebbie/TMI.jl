"""
    cache_f_and_g!(f_and_g!, dimension) -> f, grad

Wrap a combined objective/gradient callback into separate objective and gradient
functions with caching. Many solvers (e.g. `Optim` or `Ipopt`) ask for the
objective and gradient in separate calls; this wrapper avoids recomputing when
those calls use the same decision vector `x`.

# Arguments
- `f_and_g!`: Function with signature `f_and_g!(F_val, G_buf, x)` that returns
  the objective value and, when `G_buf` is not `nothing`, writes the gradient
  in-place to `G_buf`.
- `dimension::Integer`: Length of `x` (used to preallocate the cache).

# Returns
- `f`: Objective function `f(x)::Float64`.
- `grad`: Gradient function `grad(x, G)` that fills `G` in-place.
"""
function cache_f_and_g!(f_and_g!, dimension)
    last_x = nothing
    last_f = 0.0 # Will store the objective value
    last_grad = zeros(Float64, dimension) # Will store the gradient

    function update!(x_current)
        # Check if x has changed from the last time.
        # Using `!all(x .== last_x)` for array comparison.
        if last_x === nothing || !all(x_current .== last_x)
            # println("Updating cache") # Debugging line, uncomment if needed

            # Call the `f_and_g!` function. It updates `last_grad` in-place
            # and returns the objective value, which we store in `last_f`.
            computed_f_val = f_and_g!(last_f, last_grad, x_current) 
            
            # If f_and_g! returns the objective value (as expected for caching), use it.
            # Otherwise, last_f will retain its value from the prior call to f_and_g!
            # where F was a valid float.
            if !(computed_f_val === nothing || isnan(computed_f_val))
                last_f = computed_f_val
            end
            
            last_x = copy(x_current) # Store a copy of x_current
        end
    end

    function f(x)
        update!(x)
        return last_f
    end

    function grad(x, G)
        update!(x)
        G .= last_grad # Populate G with the cached gradient
    end

    return f, grad
end


function check_gradient(x, idx, f::Function, g::Function)
    ee = 1e-3
    x1 = deepcopy(x); x1[idx] += ee
    x2 = deepcopy(x); x2[idx] -= ee

    fp1 = f(x1)
    fp2 = f(x2)
    g_finite = (fp1 - fp2) / (2ee)

    g_analytic = g(x)
    return 100 * abs.((g_finite - g_analytic[idx]) ./ g_analytic[idx])
end
