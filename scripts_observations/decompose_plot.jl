# plot_decomposition.jl

using PythonPlot
include("decompose_discrete_full.jl")  # make sure this file is on your LOAD_PATH

"""
    plot_decomposition(f, g, D_old, D_new, w_old, w_new)

Compute the six‐term discrete decomposition and draw two figures:
1. A bar chart of the five individual terms.
2. A line plot comparing total_direct vs total_decomposed.
"""
function plot_decomposition(
    f::AbstractVector, g::AbstractVector,
    D_old::AbstractVector{Bool}, D_new::AbstractVector{Bool},
    w_old::AbstractVector, w_new::AbstractVector
)
    # call the decomposition
    res = decompose_discrete_full(f, g, D_old, D_new, w_old, w_new)

    show()   # display both figures
end

# Example usage:
# ---------------------------------------------------
# using DelimitedFiles
# f = rand(100); g = rand(100)
# w_old = rand(100); w_new = rand(100)
# D_old = rand(Bool, 100); D_new = rand(Bool, 100)
# plot_decomposition(f, g, D_old, D_new, w_old, w_new)
# ---------------------------------------------------
