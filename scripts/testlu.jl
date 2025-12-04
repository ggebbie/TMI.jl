using LinearSolve, IncompleteLU, SparseArrays, LinearAlgebra
using BenchmarkTools

# Function to create a random sparse matrix
function create_random_sparse_matrix(N, density=0.01)
    # Create random sparse matrix with given density
    # Make it diagonally dominant to ensure it's invertible
    A = sprand(N, N, density)
    
    # Make it diagonally dominant (ensures invertibility)
    for i in 1:N
        A[i, i] = sum(abs, A[i, :]) + 1.0
    end
    
    return A
end

println("="^70)
println("Sparse LU Factorization Performance Comparison")
println("="^70)

# Test with 10,000 x 10,000 matrix
N = 5000
density = 0.01  # 1% density (adjust if needed)

println("\n" * "="^70)
println("Matrix size: $N × $N")
println("Density: $(density * 100)%")
println("="^70)

# Create test matrix and RHS
println("\nGenerating random sparse matrix...")
A = create_random_sparse_matrix(N, density)
b = rand(N)
    
# Display sparsity info
nnz_count = nnz(A)
sparsity = 100 * (1 - nnz_count / N^2)
println("Nonzeros: $nnz_count")
println("Sparsity: $(round(sparsity, digits=2))%")
println("-"^70)

# ==================================================================
# Method 1: Full LU factorization (UMFPACK)
# ==================================================================
println("\n1. Full LU Factorization (UMFPACK):")
prob = LinearProblem(A, b)

# Benchmark factorization + solve
result_lu = @benchmark solve($prob, UMFPACKFactorization()) samples=5
t_lu = median(result_lu.times) / 1e6  # Convert to ms
println("   Time: $(round(t_lu, digits=2)) ms")

# Get solution for accuracy comparison
sol_lu = solve(prob, UMFPACKFactorization())
x_lu = sol_lu.u

# ==================================================================
# Method 2: ILU + GMRES
# ==================================================================
println("\n2. ILU + GMRES:")

# Benchmark ILU preconditioner construction
result_ilu = @benchmark ilu($A, τ=0.01) samples=5
t_ilu = median(result_ilu.times) / 1e6
println("   ILU setup: $(round(t_ilu, digits=2)) ms")

P = ilu(A, τ=0.01)

# Benchmark GMRES solve
result_gmres = @benchmark solve($prob, KrylovJL_GMRES(); 
                                 Pl=$P, reltol=1e-6, maxiter=1000) samples=5
t_gmres = median(result_gmres.times) / 1e6
println("   GMRES solve: $(round(t_gmres, digits=2)) ms")
println("   Total: $(round(t_ilu + t_gmres, digits=2)) ms")

# Get solution
sol_gmres = solve(prob, KrylovJL_GMRES(); Pl=P, reltol=1e-6, maxiter=1000)
x_gmres = sol_gmres.u

# ==================================================================
# Method 3: ILU + BiCGStab (often faster than GMRES)
# ==================================================================
println("\n3. ILU + BiCGStab:")

result_bicg = @benchmark solve($prob, KrylovJL_BICGSTAB(); 
                                Pl=$P, reltol=1e-6, maxiter=1000) samples=5
t_bicg = median(result_bicg.times) / 1e6
println("   BiCGStab solve: $(round(t_bicg, digits=2)) ms")
println("   Total: $(round(t_ilu + t_bicg, digits=2)) ms")

# Get solution
sol_bicg = solve(prob, KrylovJL_BICGSTAB(); Pl=P, reltol=1e-6, maxiter=1000)
x_bicg = sol_bicg.u

# ==================================================================
# Method 4: ILU with higher drop tolerance (faster setup)
# ==================================================================
println("\n4. ILU(τ=0.1) + BiCGStab (faster setup, weaker preconditioner):")

# Benchmark ILU preconditioner with higher drop tolerance
result_ilu_fast = @benchmark ilu($A, τ=0.1) samples=5
t_ilu_fast = median(result_ilu_fast.times) / 1e6
println("   ILU(τ=0.1) setup: $(round(t_ilu_fast, digits=2)) ms")

P_fast = ilu(A, τ=0.1)

# Benchmark BiCGStab solve with faster ILU
result_bicg_fast = @benchmark solve($prob, KrylovJL_BICGSTAB(); 
                                     Pl=$P_fast, reltol=1e-6, maxiter=1000) samples=5
t_bicg_fast = median(result_bicg_fast.times) / 1e6
println("   BiCGStab solve: $(round(t_bicg_fast, digits=2)) ms")
println("   Total: $(round(t_ilu_fast + t_bicg_fast, digits=2)) ms")

sol_bicg_fast = solve(prob, KrylovJL_BICGSTAB(); Pl=P_fast, reltol=1e-6, maxiter=1000)
x_bicg_fast = sol_bicg_fast.u

# ==================================================================
# Accuracy comparison
# ==================================================================
println("\n" * "-"^70)
println("Accuracy (relative error vs LU):")
println("   GMRES:         $(norm(x_gmres - x_lu) / norm(x_lu))")
println("   BiCGStab:      $(norm(x_bicg - x_lu) / norm(x_lu))")
println("   ILU(0.1)+BiCG: $(norm(x_bicg_fast - x_lu) / norm(x_lu))")

# ==================================================================
# Speedup summary
# ==================================================================
println("\n" * "-"^70)
println("Speedup vs Full LU:")
println("   ILU(0.01)+GMRES:   $(round(t_lu / (t_ilu + t_gmres), digits=2))x")
println("   ILU(0.01)+BiCGStab: $(round(t_lu / (t_ilu + t_bicg), digits=2))x")
println("   ILU(0.1)+BiCGStab:  $(round(t_lu / (t_ilu_fast + t_bicg_fast), digits=2))x")

# ==================================================================
# Multiple solves test (reusing factorization/preconditioner)
# ==================================================================
println("\n" * "-"^70)
println("Multiple Solves Test (100 different RHS, reusing factorization):")

# Full LU - factorize once, solve many
F_lu = lu(A)
result_multi_lu = @benchmark begin
    for i in 1:100
        x = $F_lu \ rand($N)
    end
end samples=3
t_multi_lu = median(result_multi_lu.times) / 1e6
println("   Full LU: $(round(t_multi_lu, digits=2)) ms")

# ILU+BiCGStab - compute ILU once, solve many
prob_temp = LinearProblem(A, rand(N))
cache = init(prob_temp, KrylovJL_BICGSTAB(); Pl=P, reltol=1e-6)
result_multi_ilu = @benchmark begin
    for i in 1:100
        $cache.b = rand($N)
        solve!($cache)
    end
end samples=3
t_multi_ilu = median(result_multi_ilu.times) / 1e6
println("   ILU+BiCGStab: $(round(t_multi_ilu, digits=2)) ms")
println("   Speedup: $(round(t_multi_lu / t_multi_ilu, digits=2))x")

println("\n" * "="^70)
println("Summary:")
println("="^70)
println("• Full LU (UMFPACK): Most accurate, good baseline")
println("• ILU(τ=0.01)+BiCGStab: Best balance of speed and accuracy")
println("• ILU(τ=0.1): Faster setup but weaker preconditioning")
println("• For multiple solves: ILU+iterative wins by even larger margins")
println("• Recommendation: Use ILU(τ=0.01)+BiCGStab for best balance")
println("="^70)