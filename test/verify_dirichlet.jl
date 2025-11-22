using MultiGridBarrier3d
using LinearAlgebra
using Test

function verify_dirichlet()
    println("Verifying Dirichlet Projection...")
    
    # Create geometry
    k = 3
    L = 2
    geo = fem3d(Float64; L=L, k=k)
    
    # Get Dirichlet subspace matrix D
    # D maps interior DOFs to full vector
    D = geo.subspaces[:dirichlet][L] # Fine grid
    
    # Check orthogonality: D'D should be Identity
    DtD = D' * D
    diff_I = norm(DtD - I, Inf)
    println("||D'D - I||_inf = $diff_I")
    
    # Define u that is zero on boundary of [-1, 1]^3
    # User suggested cos(pi*x)*... but that is non-zero at +/- 1.
    # We use cos(pi*x/2)*cos(pi*y/2)*cos(pi*z/2) which is zero at +/- 1.
    
    u = zeros(size(geo.x, 1))
    for i in 1:size(geo.x, 1)
        x, y, z = geo.x[i, :]
        u[i] = cos(π*x/2) * cos(π*y/2) * cos(π*z/2)
    end
    
    # Check boundary values are indeed small
    boundary_indices, _, node_map = MultiGridBarrier3d.get_boundary_nodes(geo.x, k)
    boundary_set = Set(boundary_indices)
    boundary_element_indices = [i for i in 1:length(u) if node_map[i] in boundary_set]
    
    max_boundary = norm(u[boundary_element_indices], Inf)
    println("Max u on boundary: $max_boundary")
    
    # Project: P = D * (D'D)^-1 * D'
    # D'D is diagonal matrix of multiplicities.
    
    # Compute inv(D'D) efficiently
    # Since D'D is diagonal, we can just invert the diagonal elements.
    # But for this test, explicit inverse is fine as size is small.
    
    invDtD = inv(Matrix(DtD)) # Convert sparse to dense for inv, or use sparse solve
    
    P = D * invDtD * D'
    
    u_proj = P * u
    
    # Check error
    err = norm(u - u_proj, Inf)
    println("||u - P*u||_inf = $err")
    
    if err < 1e-14
        println("VERIFICATION SUCCESSFUL")
    else
        println("VERIFICATION FAILED")
    end
end

verify_dirichlet()
