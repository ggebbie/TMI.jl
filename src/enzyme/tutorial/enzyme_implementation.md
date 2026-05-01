# AD Implementation with Enzyme.jl in TMI

This document covers the practical application of Automatic Differentiation (AD) in TMI.jl using the Enzyme.jl library. It assumes you are familiar with the concepts of reverse-mode AD as explained in `analytical_ad.md`.

## Part 1: Practical AD with Enzyme.jl

Enzyme is the tool that performs reverse-mode AD on Julia code. To make the connection between the theory and the code explicit, let's look at the terminology.

The **adjoint** (or sensitivity) from the theory, denoted with a bar, represents the partial derivative of the final output with respect to an intermediate variable. In Enzyme, this concept is realized as a mutable storage location (a variable) that accumulates the gradient value.

For our example function $L(x, y) = \sin(x) + x \cdot y$:
- The theoretical adjoint $\bar{x} = \frac{\partial L}{\partial x}$ corresponds to the Julia variable `dx`.
- The theoretical adjoint $\bar{y} = \frac{\partial L}{\partial y}$ corresponds to the Julia variable `dy`.

Enzyme's `autodiff` function calculates the values of these adjoints and places them in the corresponding storage variables.

Here is the code to differentiate the function from the theory document:

```julia
using Enzyme
import .sin, .cos

# The function to differentiate
h(x, y) = sin(x) + x*y

# Input values
x = π/2
y = 2.0

# Gradient storage variables, corresponding to x̄ and ȳ
dx = 0.0
dy = 0.0

# Perform reverse-mode AD.
# `Duplicated(x, dx)` marks x as an active variable whose gradient
# (adjoint) will be accumulated into `dx`.
autodiff(Reverse, h, Duplicated(x, dx), Duplicated(y, dy))

println("∂L/∂x = ", dx) # Expected: 2.0
println("∂L/∂y = ", dy) # Expected: 1.5707... (π/2)
```

## Part 2: Custom Derivatives for TMI's Key Operations

Sometimes, a function's implementation is too complex for Enzyme to differentiate automatically (e.g., due to data copies, complex memory layouts, or non-trivial linear algebra). In these cases, we give Enzyme a **custom rule** that explicitly provides the forward and reverse passes, performing the same logic derived from the chain rule in the theory document.

TMI uses custom rules for a few key mathematical patterns.

### Pattern 1: Scatter/Gather (Vector-to-Structure)

A common operation in TMI is to take a flat vector of values and **scatter** it into a structured object, like a gridded field.

- **Forward Operation:** $y = 	ext{unvec}(	ext{template}, 	ext{uvec})$. A vector `uvec` containing values for "wet" points is scattered into a full structured field `y`. Non-wet points in `y` are filled with default values from the `template`.
- **Reverse Operation (Adjoint):** The derivative must do the opposite. It **gathers** sensitivities from the corresponding wet points of the structure's adjoint, $\bar{y}$ (also denoted $g_y$), and accumulates them into the vector's adjoint, $\bar{uvec}$ (or $g_{uvec}$).

### Pattern 2: The Sparse Linear Solve

The most important custom rule is for the sparse linear solve, $A \cdot c = d$. The reverse pass for this operation is a perfect example of the **adjoint method**.

- **Forward Operation:** Solve $A \cdot c = d$ for $c$ and save it for the reverse pass.
- **Reverse Operation (Adjoint):** The reverse pass receives $g_c$ (the adjoint of $c$). To compute the adjoints of $A$ and $d$, it:
    1.  Solves the **adjoint system**: $A^T \cdot \lambda = g_c$ to find the adjoint variable $\lambda$. This is the key vector-Jacobian product for this operation.
    2.  The resulting gradients are then $g_d = \lambda$ and the contribution to the gradient of $A$ is $g_A = -\lambda \cdot c^T$.

## Part 3: How the Rules Come Together

Each of these patterns is implemented as a custom rule in a dedicated `.jl` file (e.g., `ldiv_field.jl`, `watermassmatrix.jl`).

The file `TMIEnzymeRules.jl` is a module that simply collects and registers all of these custom rules, making them available to Enzyme when it differentiates the top-level TMI functions.
