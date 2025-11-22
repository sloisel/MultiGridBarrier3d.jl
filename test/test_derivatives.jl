using MultiGridBarrier3d
using Test
using LinearAlgebra

@testset "Derivative Tests" begin
    # Use k=3 (cubic elements) so that tricubic polynomials are exactly represented
    k = 3
    L = 2
    geo = fem3d(Float64; L=L, k=k)
    
    println("Testing Derivatives (k=$k, L=$L)")
    
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
    
    function poly_dx(x, y, z)
        val = 0.0
        for l in 1:4, j in 1:4, i in 2:4
            val += coeffs[i,j,l] * (i-1) * x^(i-2) * y^(j-1) * z^(l-1)
        end
        return val
    end
    
    function poly_dy(x, y, z)
        val = 0.0
        for l in 1:4, j in 2:4, i in 1:4
            val += coeffs[i,j,l] * (j-1) * x^(i-1) * y^(j-2) * z^(l-1)
        end
        return val
    end
    
    function poly_dz(x, y, z)
        val = 0.0
        for l in 2:4, j in 1:4, i in 1:4
            val += coeffs[i,j,l] * (l-1) * x^(i-1) * y^(j-1) * z^(l-2)
        end
        return val
    end
    
    # Evaluate on fine grid
    x_fine = geo.x
    n_nodes = size(x_fine, 1)
    
    u = zeros(n_nodes)
    u_exact_dx = zeros(n_nodes)
    u_exact_dy = zeros(n_nodes)
    u_exact_dz = zeros(n_nodes)
    
    for i in 1:n_nodes
        x, y, z = x_fine[i, :]
        u[i] = poly(x, y, z)
        u_exact_dx[i] = poly_dx(x, y, z)
        u_exact_dy[i] = poly_dy(x, y, z)
        u_exact_dz[i] = poly_dz(x, y, z)
    end
    
    # Apply operators
    Dx = geo.operators[:dx]
    Dy = geo.operators[:dy]
    Dz = geo.operators[:dz]
    
    u_dx = Dx * u
    u_dy = Dy * u
    u_dz = Dz * u
    
    # Compare
    err_dx = norm(u_dx - u_exact_dx, Inf)
    err_dy = norm(u_dy - u_exact_dy, Inf)
    err_dz = norm(u_dz - u_exact_dz, Inf)
    
    println("Error dx: $err_dx")
    println("Error dy: $err_dy")
    println("Error dz: $err_dz")
    
    # Tolerances
    # Since k=3 matches the polynomial degree, error should be close to machine epsilon
    # However, there might be some conditioning issues or floating point accumulation.
    # 1e-12 seems safe.
    
    @test err_dx < 1e-12
    @test err_dy < 1e-12
    @test err_dz < 1e-12
end

@testset "Distorted Element Tests" begin
    println("Testing Distorted Element (Non-constant Jacobian)")
    
    # Define a distorted Q1 element (8 vertices)
    # Cube is [-1, 1]^3. Let's perturb it.
    # Vertices ordered as tensor product:
    # (-1,-1,-1), (1,-1,-1), (-1,1,-1), (1,1,-1), ...
    
    K_q1 = Float64[
        -1.0 -1.0 -1.0; # 1
         1.2 -0.9 -1.1; # 2 (perturbed)
        -0.9  1.1 -0.9; # 3
         1.0  1.0 -1.0; # 4
        -1.1 -1.1  1.0; # 5
         1.0 -1.0  1.0; # 6
        -1.0  1.0  1.1; # 7
         0.9  0.9  0.9  # 8
    ]
    
    k = 3
    x_fine = MultiGridBarrier3d.promote_to_Qk(K_q1, k)
    
    # Build operators manually since we don't have a full Geometry object
    # But build_operators expects matrix x.
    ops = MultiGridBarrier3d.build_operators(x_fine, k)
    Dx = ops[:dx]
    Dy = ops[:dy]
    Dz = ops[:dz]
    
    # Test 1: Linear function in physical space
    # u = c1*x + c2*y + c3*z + c0
    # This is exactly representable in Qk (k>=1) even on distorted mesh
    # because x(xi), y(xi), z(xi) are trilinear, so u(xi) is trilinear.
    
    c1, c2, c3, c0 = 2.0, -3.0, 0.5, 1.0
    
    u = c1 .* x_fine[:, 1] .+ c2 .* x_fine[:, 2] .+ c3 .* x_fine[:, 3] .+ c0
    
    u_dx = Dx * u
    u_dy = Dy * u
    u_dz = Dz * u
    
    err_dx = norm(u_dx .- c1, Inf)
    err_dy = norm(u_dy .- c2, Inf)
    err_dz = norm(u_dz .- c3, Inf)
    
    println("Linear function errors: dx=$err_dx, dy=$err_dy, dz=$err_dz")
    
    @test err_dx < 1e-11
    @test err_dy < 1e-11
    @test err_dz < 1e-11
    
    # Test 2: Random function in reference space
    # u(xi) is a random polynomial.
    # We verify that Dx computes J^{-T} * grad_xi u
    
    # We need to access the internal Jacobian logic or replicate it.
    # Let's replicate the chain rule calculation for one point.
    
    # Random vector u (represents random polynomial in reference space)
    u_ref = rand(size(x_fine, 1))
    
    # Compute reference derivatives
    nodes_1d = MultiGridBarrier3d.chebyshev_nodes(k)
    D_1d = MultiGridBarrier3d.derivative_matrix_1d(nodes_1d)
    I_1d = I(k+1)
    
    D_xi = kron(I_1d, kron(I_1d, D_1d))
    D_eta = kron(I_1d, kron(D_1d, I_1d))
    D_zeta = kron(D_1d, kron(I_1d, I_1d))
    
    u_xi = D_xi * u_ref
    u_eta = D_eta * u_ref
    u_zeta = D_zeta * u_ref
    
    # Physical derivatives from operator
    u_x_op = Dx * u_ref
    u_y_op = Dy * u_ref
    u_z_op = Dz * u_ref
    
    # Check at each point
    max_err = 0.0
    
    # Re-compute metrics locally to verify
    # We need the mapping x(xi).
    # x_fine contains the mapped points.
    # We can compute x_xi, x_eta, x_zeta using D_xi, etc. on x_fine.
    
    x_coords = x_fine[:, 1]
    y_coords = x_fine[:, 2]
    z_coords = x_fine[:, 3]
    
    x_xi = D_xi * x_coords
    x_eta = D_eta * x_coords
    x_zeta = D_zeta * x_coords
    
    y_xi = D_xi * y_coords
    y_eta = D_eta * y_coords
    y_zeta = D_zeta * y_coords
    
    z_xi = D_xi * z_coords
    z_eta = D_eta * z_coords
    z_zeta = D_zeta * z_coords
    
    for i in 1:size(x_fine, 1)
        J = [x_xi[i] x_eta[i] x_zeta[i];
             y_xi[i] y_eta[i] y_zeta[i];
             z_xi[i] z_eta[i] z_zeta[i]]
             
        # grad_x = J^{-T} * grad_xi
        # We computed grad_x_op using the operator.
        # We want to check if grad_x_op matches J^{-T} * grad_xi.
        
        grad_xi = [u_xi[i], u_eta[i], u_zeta[i]]
        grad_x_expected = inv(J') * grad_xi
        
        grad_x_op = [u_x_op[i], u_y_op[i], u_z_op[i]]
        
        err = norm(grad_x_expected - grad_x_op)
        max_err = max(max_err, err)
    end
    
    println("Max chain rule consistency error: $max_err")
    @test max_err < 1e-11
end
