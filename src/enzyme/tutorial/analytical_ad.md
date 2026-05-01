# The Theory of Automatic Differentiation

This document provides a textbook-style introduction to the theory of Automatic Differentiation (AD), starting from the first principles of calculus.

## Motivation: Why We Need Gradients

Many scientific problems are optimization problems. In TMI, for instance, we want to find the ocean circulation state that best fits a set of observations. To do this efficiently, we use gradient-based optimization algorithms, which require us to compute the **gradient**: the vector of all partial derivatives of a cost function with respect to its many input parameters.

For a complex computer model, this cost function is an incredibly long and nested series of simple operations. Manually deriving the symbolic gradient is practically impossible. **Automatic Differentiation (AD)** is a technique for computing the exact numerical value of a function's gradient without ever deriving the full symbolic expression.

## The Core Idea: The Chain Rule on a Computer Program

Let's consider a simple function that mimics a computer program: $L(x, y) = \sin(x) + x \cdot y$. This can be seen as a sequence of nested operations:
- Let $a(x) = \sin(x)$
- Let $b(x, y) = x \cdot y$
- Then $L(x,y) = a(x) + b(x,y)$

From multivariable calculus, the chain rule tells us how to find the partial derivative of $L$ with respect to $x$:
$$ \frac{\partial L}{\partial x} = \frac{\partial L}{\partial a}\frac{\partial a}{\partial x} + \frac{\partial L}{\partial b}\frac{\partial b}{\partial x} $$
This expression forms the basis for both modes of Automatic Differentiation. The difference between them is the order in which we evaluate and group the terms.

Let's evaluate the gradient at $(x, y) = (\frac{\pi}{2}, 2)$.

### Strategy 1: Forward Accumulation (Forward-Mode AD)

The first strategy is to evaluate the chain rule from left to right, in the same order as the program's execution. This approach calculates how a change in a single **input** propagates forward to affect all intermediate and final values.

Let's compute $\frac{\partial L}{\partial x}$. To do this, we compute the derivatives of all intermediate variables ($a, b$) with respect to $x$ as we go.
1.  **Input Derivatives**: We seed our calculation with the derivative of the inputs with respect to $x$:
    - $\frac{\partial x}{\partial x} = 1$
    - $\frac{\partial y}{\partial x} = 0$
2.  **Step-by-step Evaluation**: We compute the values and derivatives for each step of the program:
    - $a = \sin(x) = \sin(\frac{\pi}{2}) = 1$
    - $\frac{\partial a}{\partial x} = \cos(x) \cdot \frac{\partial x}{\partial x} = \cos(\frac{\pi}{2}) \cdot 1 = 0$
    - $b = x \cdot y = \frac{\pi}{2} \cdot 2 = \pi$
    - $\frac{\partial b}{\partial x} = y \cdot \frac{\partial x}{\partial x} + x \cdot \frac{\partial y}{\partial x} = 2 \cdot 1 + \frac{\pi}{2} \cdot 0 = 2$
3.  **Final Combination**: Now we can evaluate the final derivative:
    - $L = a + b = 1 + \pi$
    - $\frac{\partial L}{\partial x} = \frac{\partial a}{\partial x} + \frac{\partial b}{\partial x} = 0 + 2 = 2$

This process of evaluating the primal values and their derivatives together in one pass from inputs to outputs is **Forward-Mode AD**. It efficiently computes the derivative of the output with respect to *one* input. To find $\frac{\partial L}{\partial y}$, we would need a completely separate pass. For a function with thousands of inputs, this would require thousands of passes, making it prohibitively expensive.

### Strategy 2: Reverse Accumulation (Reverse-Mode AD)

The second strategy is to evaluate the chain rule from right to left, starting from the final output. This approach calculates how a change in the final **output** is affected by every intermediate and input variable.

This strategy consists of two stages:
1.  **Forward Pass**: First, run the original program from start to finish to compute the value of every intermediate variable. This is the "primal" pass.
    - $x = \pi/2$, $y = 2$
    - $a = \sin(x) = 1$
    - $b = x \cdot y = \pi$
    - $L = a + b = 1 + \pi$
2.  **Backward Pass**: Now, we propagate sensitivities backward. This sensitivity is called the **adjoint**, denoted with a bar ($\bar{u}$), and is defined as the partial derivative of the final output with respect to that variable, $\bar{u} = \frac{\partial L}{\partial u}$.
    - **Start at the end**: The sensitivity of the output with respect to itself is $\bar{L} = \frac{\partial L}{\partial L} = 1$.
    - **Propagate backward to $a$ and $b$**:
      $\bar{a} = \frac{\partial L}{\partial a} = \bar{L} \cdot \frac{\partial(a+b)}{\partial a} = 1 \cdot 1 = 1$
      $\bar{b} = \frac{\partial L}{\partial b} = \bar{L} \cdot \frac{\partial(a+b)}{\partial b} = 1 \cdot 1 = 1$
    - **Propagate backward to $x$ and $y$**: We use the chain rule again. The total sensitivity of an input is the sum of sensitivities from all paths it affects.
      $\bar{y} = \frac{\partial L}{\partial y} = \bar{b} \cdot \frac{\partial b}{\partial y} = 1 \cdot x = \frac{\pi}{2}$
      $\bar{x} = \frac{\partial L}{\partial x} = (\bar{a} \cdot \frac{\partial a}{\partial x}) + (\bar{b} \cdot \frac{\partial b}{\partial x}) = (1 \cdot \cos(x)) + (1 \cdot y) = 0 + 2 = 2$

This process, consisting of a forward pass for values followed by a backward pass for sensitivities, is **Reverse-Mode AD**. With just one forward and one backward pass, we computed the entire gradient vector $[\bar{x}, \bar{y}] = [2, \frac{\pi}{2}]$. This remarkable efficiency is why reverse-mode is the only feasible choice for large-scale optimization.

### The Matrix Formulation (for completeness)

The reverse-mode process can be written more formally using **Jacobian matrices**. The rule $\bar{	ext{input}} = \bar{	ext{output}} \cdot \frac{\partial 	ext{output}}{\partial 	ext{input}}$ is equivalent to a **transposed-Jacobian-vector product**. This is the mathematical foundation of the "adjoint method" used in TMI's custom rules. For example, solving $A^T\lambda=g_c$ is the high-performance way to compute the vector-Jacobian product for the linear solve operation.
