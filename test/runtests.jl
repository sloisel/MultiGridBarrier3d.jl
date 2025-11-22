using MultiGridBarrier3d
using LinearAlgebra
using SparseArrays
using Test

function test_polynomial_reproduction()
    k = 3
    L = 2
    g = fem3d(L=L, k=k)
    
    println("Testing Polynomial Reproduction (k=$k, L=$L)")
    
    # Random tricubic polynomial
    # P(x,y,z) = sum_{i,j,l=0..3} c_{ijl} x^i y^j z^l
    coeffs = rand(4, 4, 4)
    
    function poly(x, y, z)
        val = 0.0
        for l in 1:4, j in 1:4, i in 1:4
            val += coeffs[i,j,l] * x^(i-1) * y^(j-1) * z^(l-1)
        end
        return val
    end
    
    # Evaluate on coarse grid (L=1)
    # We need to access the coarse mesh. It's not stored in g directly, 
    # but we can recreate it or access via g.refine?
    # Actually, g only stores fine mesh.
    # But create_geometry builds the hierarchy.
    # Let's modify create_geometry to store all meshes? 
    # Or just test on fine grid for now.
    
    # Evaluate on fine grid
    x_fine = g.x
    u_fine = zeros(size(x_fine, 1))
    for i in 1:size(x_fine, 1)
        u_fine[i] = poly(x_fine[i, 1], x_fine[i, 2], x_fine[i, 3])
    end
    
    # Check interpolation at random points
    n_points = 10
    max_err = 0.0
    for _ in 1:n_points
        pt = 2.0 .* rand(3) .- 1.0 # Random point in [-1, 1]^3
        val_exact = poly(pt[1], pt[2], pt[3])
        val_interp, found = MultiGridBarrier3d.evaluate_field(g, u_fine, pt)
        
        if found
            err = abs(val_interp - val_exact)
            max_err = max(max_err, err)
        else
            println("Point $pt not found in mesh!")
        end
    end
    
    println("Max interpolation error: $max_err")
    @test max_err < 1e-10
end

function test_multigrid_ops()
    k = 3
    L = 2
    g = fem3d(L=L, k=k)
    
    println("Testing Multigrid Operators")
    
    # Check refine * coarsen approx Identity? No.
    # Check coarsen * refine = Identity (on coarse grid space)
    # P: coarse -> fine
    # R: fine -> coarse
    # R * P should be Identity on coarse grid
    
    P = g.refine[1] # Level 1 to 2
    R = g.coarsen[1] # Level 2 to 1
    
    I_coarse = R * P
    
    # Size of coarse grid
    n_coarse = size(I_coarse, 1)
    
    # Should be identity
    err = norm(I_coarse - I, Inf)
    println("||R*P - I||_inf = $err")
    @test err < 1e-10
end

function test_dirichlet()
    k = 3
    L = 1
    g = fem3d(L=L, k=k)
    
    println("Testing Dirichlet Boundary")
    
    D = g.subspaces[:dirichlet][1]
    
    # D projects from interior parameters to full vector
    # u = D * v should be zero on boundary
    
    n_interior = size(D, 2)
    v = rand(n_interior)
    u = D * v
    
    # Identify boundary nodes
    boundary_unique_indices, _, node_map = MultiGridBarrier3d.get_boundary_nodes(g.x, k)
    boundary_set = Set(boundary_unique_indices)
    
    # Find indices in u (element-wise) that correspond to boundary nodes
    boundary_element_indices = [i for i in 1:length(u) if node_map[i] in boundary_set]
    
    boundary_vals = u[boundary_element_indices]
    err = norm(boundary_vals, Inf)
    println("Max boundary value: $err")
    @test err < 1e-14
end

function test_projection()
    k = 3
    L = 2
    g = fem3d(Float64; L=L, k=k)
    
    println("Testing Projection (k=$k, L=$L)")
    
    # Function u with zero boundary data
    # u = cos(pi*x/2)*cos(pi*y/2)*cos(pi*z/2) (zero at +/- 1)
    
    u = zeros(size(g.x, 1))
    for i in 1:size(g.x, 1)
        x, y, z = g.x[i, :]
        u[i] = cos(pi*x/2) * cos(pi*y/2) * cos(pi*z/2)
    end
    
    D = g.subspaces[:dirichlet][L]
    
    # Project u onto range(D)
    # P = D * (D'D)^-1 * D'
    # D'D is diagonal matrix of multiplicities.
    
    DtD = D' * D
    invDtD = inv(Matrix(DtD))
    
    P = D * invDtD * D'
    u_proj = P * u
    
    # u should be in range(D) because it is zero on boundary
    # But wait, u is zero on boundary analytically.
    # Numerically, nodes on boundary should have u=0 (approx).
    
    # Check boundary values of u
    boundary_unique_indices, _, node_map = MultiGridBarrier3d.get_boundary_nodes(g.x, k)
    boundary_set = Set(boundary_unique_indices)
    boundary_element_indices = [i for i in 1:length(u) if node_map[i] in boundary_set]
    
    println("Max u on boundary: $(norm(u[boundary_element_indices], Inf))")
    
    # u_proj should match u
    err = norm(u - u_proj, Inf)
    println("Projection error: $err")
    @test err < 1e-14
end

test_polynomial_reproduction()
test_multigrid_ops()
test_dirichlet()
test_projection()

include("test_plotting_api.jl")
include("test_slicing.jl")
include("test_smart_volume.jl")
include("test_fem3d.jl")
include("test_derivatives.jl")
include("test_quadrature.jl")
include("test_point_location.jl")
include("test_solver.jl")
