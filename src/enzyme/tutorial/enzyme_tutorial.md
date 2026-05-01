# A Guide to Automatic Differentiation in TMI with Enzyme.jl

Welcome! This document is your complete guide to using Automatic Differentiation (AD) with Enzyme.jl within the TMI (Tracer Model I) codebase. It's written like a lecture, starting with the mathematical first principles and building up to the specific, real-world code used in TMI.

## Part 1: The Theory of Automatic Differentiation

### Motivation: Why We Need Gradients

Many scientific problems are optimization problems. In TMI, for instance, we want to find the ocean circulation state that best fits a set of observations. To do this efficiently, we use gradient-based optimization algorithms, which require us to compute the **gradient** of our model's cost function with respect to its many input parameters.

Manually deriving these gradients for a complex model is tedious and error-prone. **Automatic Differentiation (AD)** is a technique for computing gradients of functions written as computer programs, accurately and automatically. It is the engine that makes large-scale optimization of complex models possible.

### The Mathematics of AD: The Jacobian Chain Rule

At its heart, AD is a mechanical application of the chain rule. Let's consider a function `L(x, y) = sin(x) + x*y` evaluated at `(x, y) = (π/2, 2)`. We can express this as a composition of vector functions:

1. Let the inputs be a vector `v = [x, y]ᵀ`.
2. The first stage computes intermediate variables: `w = g(v) = [sin(x), x*y]ᵀ`.
3. The final stage computes the result: `L = h(w) = w₁ + w₂`.

The chain rule in matrix form states that the overall gradient `∂L/∂v` is the product of the Jacobians of the individual stages: `∂L/∂v = (∂h/∂w) * (∂g/∂v)`.

#### Forward-Mode AD: Jacobian-Vector Products

Forward-mode AD computes the full gradient by evaluating the chain rule one input at a time. This is equivalent to performing a **Jacobian-vector product** `J*v` at each step.

Let's find `∂L/∂x`. We define a "seed" vector `v̇ = [∂x/∂x, ∂y/∂x]ᵀ = [1, 0]ᵀ`.

1. **Stage 1:** `ẇ = (∂g/∂v) * v̇`
   The Jacobian `∂g/∂v` is:
   `[[∂g₁/∂x, ∂g₁/∂y], [∂g₂/∂x, ∂g₂/∂y]] = [[cos(x), 0], [y, x]]`
   At our evaluation point, this is `[[0, 0], [2, π/2]]`.
   So, `ẇ = [[0, 0], [2, π/2]] * [1, 0]ᵀ = [0, 2]ᵀ`.

2. **Stage 2:** `L̇ = (∂h/∂w) * ẇ`
   The Jacobian `∂h/∂w` is `[∂h/∂w₁, ∂h/∂w₂] = [1, 1]`.
   So, `L̇ = [1, 1] * [0, 2]ᵀ = 2`.

The result is `∂L/∂x = 2`. To find `∂L/∂y`, we would need a separate pass with a different seed vector `v̇ = [0, 1]ᵀ`. For a model with thousands of inputs like TMI, this is prohibitively expensive.

#### Reverse-Mode AD: Vector-Jacobian Products

Reverse-mode AD computes the full gradient in a single pass by propagating sensitivities backward. This is equivalent to performing a **vector-Jacobian product** `vᵀ*J` (or `Jᵀ*v` for column vectors) at each step.

We start with the sensitivity of the output with respect to itself, `L̄ = 1`. The bar (⁻) denotes the sensitivity, or **adjoint**. The rule is `input̄ = (Jacobian)ᵀ * output̄`.

1. **Stage 2:** `w̄ = (∂h/∂w)ᵀ * L̄`
   The Jacobian `∂h/∂w` was `[1, 1]`. Its transpose is `[1, 1]ᵀ`.
   `w̄ = [1, 1]ᵀ * 1 = [1, 1]ᵀ`. This gives us the sensitivities of the intermediate variables.

2. **Stage 1:** `v̄ = (∂g/∂v)ᵀ * w̄`
   The Jacobian `∂g/∂v` was `[[0, 0], [2, π/2]]`. Its transpose is `[[0, 2], [0, π/2]]`.
   `v̄ = [[0, 2], [0, π/2]] * [1, 1]ᵀ = [0*1 + 2*1, 0*1 + (π/2)*1]ᵀ = [2, π/2]ᵀ`.

The final vector of adjoints is `v̄ = [x̄, ȳ]ᵀ = [∂L/∂x, ∂L/∂y]ᵀ = [2, π/2]ᵀ`.

With a **single backward pass**, we computed **all** partial derivatives. This is why for a model with a vast number of inputs and a single scalar output, only **reverse-mode AD** is computationally feasible. The concept of multiplying by the **transposed Jacobian** is the mathematical foundation for the "adjoint method" used throughout TMI.

## Part 2: Practical AD with Enzyme.jl

Enzyme is the tool that performs this reverse-mode AD on Julia code. The `v̄ = [x̄, ȳ]ᵀ` adjoints from our analytical example map directly to the gradient storage variables `dx` and `dy`.

Here is the code to differentiate `h(x, y) = sin(x) + x*y`, which we analyzed in Part 1.

```julia
using Enzyme
import .sin, .cos

# The function to differentiate
h(x, y) = sin(x) + x*y

# Input values
x = π/2
y = 2.0

# Gradient storage (maps to x̄ and ȳ from the theory)
dx = 0.0
dy = 0.0

# Perform reverse-mode AD.
# Duplicated(x, dx) marks x as an active variable
# and tells Enzyme to store its gradient in dx.
autodiff(Reverse, h, Duplicated(x, dx), Duplicated(y, dy))

println("∂L/∂x = ", dx) # Expected: 2.0
println("∂L/∂y = ", dy) # Expected: 1.5707... (π/2)
```

## Part 3: Custom Derivatives for TMI's Key Operations

Sometimes, a function's implementation is too complex for Enzyme to differentiate automatically (e.g., due to data copies or non-trivial linear algebra). In these cases, we give Enzyme a **custom rule** that explicitly provides the forward and reverse passes, which perform the Jacobian-vector products described in Part 1.

TMI uses custom rules for a few key mathematical patterns.

### Pattern 1: Scatter/Gather (Vector-to-Structure)
- **Forward:** `y = unvec(template, uvec)` maps a vector `uvec` to the wet points of a larger structure `y`.
- **Reverse:** The derivative must do the opposite. It **gathers** the sensitivities from the structure `g_y` back into the vector `g_uvec`. This is the **adjoint** of the scatter, equivalent to a transposed-Jacobian-vector product: `g_uvec = Jᵀ * g_y`.

### Pattern 2: The Sparse Linear Solve

The most important rule is for the sparse linear solve `A * c = d`. Here, the transposed-Jacobian-vector product is embodied by the **adjoint method**.
- **Forward:** Solve `A * c = d` and store the result `c`.
- **Reverse:** The reverse pass gets `g_c` (the sensitivity of the objective to `c`). It then:
    1.  Solves the **adjoint system**: `Aᵀ * λ = g_c` to find the **adjoint variable** `λ`. This is the key transposed-Jacobian-vector product.
    2.  Calculates the gradients: `g_d = λ` and `g_A = -λ * cᵀ`.

#### In-Depth: The Derivation
1. Primal: `A * c = d`.
2. Differential: `dA * c + A * dc = dd`.
3. Rearrange: `dc = A⁻¹(dd - dA * c)`.
4. Objective change `δJ = g_cᵀ * dc = g_cᵀ * A⁻¹ * dd - g_cᵀ * A⁻¹ * dA * c`.
5. Define adjoint `λ = (A⁻¹)ᵀ * g_c` (or `Aᵀ * λ = g_c`). This is the vector-Jacobian product.
6. Substitute `λ`: `δJ = λᵀ * dd - λᵀ * dA * c`.
7. By inspection, `g_d = λ` and the update to `g_A` is `-λ * cᵀ`.

## Part 4: How the Rules Come Together

Each of these patterns is implemented as a custom rule file (`ldiv_field.jl`, `watermassmatrix.jl`, etc.).

The file `TMIEnzymeRules.jl` is the module that collects and registers all these custom rules, making them available to Enzyme when it differentiates the top-level TMI functions. This keeps the implementation clean, organized, and performant.
