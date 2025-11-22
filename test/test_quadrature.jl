using MultiGridBarrier3d
using Test
using LinearAlgebra

@testset "Quadrature Tests" begin
    k = 3
    L = 2 # We need at least 2 levels to trigger the refinement loop where weights are computed
    geo = fem3d(Float64; L=L, k=k)
    
    println("Testing Quadrature (k=$k, L=$L)")
    
    # Fine grid weights and nodes
    w = geo.w
    x = geo.x
    
    # Test 1: Integrate 1 on the domain
    # Domain is [-1, 1]^3, volume = 8
    vol = sum(w)
    println("Volume: $vol (expected 8.0)")
    @test abs(vol - 8.0) < 1e-12
    
    # Test 2: Integrate x
    # Integral of x on [-1, 1] is 0.
    int_x = sum(w .* x[:, 1])
    println("Integral x: $int_x (expected 0.0)")
    @test abs(int_x) < 1e-12
    
    # Test 3: Integrate x^2
    # Integral x^2 on [-1, 1] is 2/3.
    # Integral x^2 on [-1, 1]^3 is (2/3) * 2 * 2 = 8/3 = 2.666...
    int_x2 = sum(w .* x[:, 1].^2)
    println("Integral x^2: $int_x2 (expected $(8/3))")
    @test abs(int_x2 - 8/3) < 1e-12
    
    # Test 4: Distorted Element
    println("Testing Distorted Element Quadrature")
    
    # Define a distorted Q1 element (8 vertices)
    # We'll use a simple shear/stretch to make it easy to compute exact volume.
    # Map: x = xi + 0.5*eta, y = eta, z = zeta
    # Jacobian = [1 0.5 0; 0 1 0; 0 0 1], det(J) = 1.
    # Volume should still be 8.
    
    K_q1 = Float64[
        -1.5 -1.0 -1.0; # (-1,-1,-1) -> x = -1 - 0.5 = -1.5
         0.5 -1.0 -1.0; # ( 1,-1,-1) -> x =  1 - 0.5 =  0.5
        -0.5  1.0 -1.0; # (-1, 1,-1) -> x = -1 + 0.5 = -0.5
         1.5  1.0 -1.0; # ( 1, 1,-1) -> x =  1 + 0.5 =  1.5
        -1.5 -1.0  1.0; # ... z=1
         0.5 -1.0  1.0;
        -0.5  1.0  1.0;
         1.5  1.0  1.0
    ]
    
    x_fine = MultiGridBarrier3d.promote_to_Qk(K_q1, k)
    
    # We need to compute weights manually for this single element mesh
    # Re-use logic from create_geometry
    
    n_nodes_per_elem = (k+1)^3
    _, w_ref_elem = MultiGridBarrier3d.cube_mesh(k)
    
    nodes_1d = MultiGridBarrier3d.chebyshev_nodes(k)
    D_1d = MultiGridBarrier3d.derivative_matrix_1d(nodes_1d)
    I_1d = I(k+1)
    
    D_xi = kron(I_1d, kron(I_1d, D_1d))
    D_eta = kron(I_1d, kron(D_1d, I_1d))
    D_zeta = kron(D_1d, kron(I_1d, I_1d))
    
    w_distorted = zeros(size(x_fine, 1))
    
    x_elem = x_fine[:, 1]
    y_elem = x_fine[:, 2]
    z_elem = x_fine[:, 3]
    
    x_xi = D_xi * x_elem
    x_eta = D_eta * x_elem
    x_zeta = D_zeta * x_elem
    
    y_xi = D_xi * y_elem
    y_eta = D_eta * y_elem
    y_zeta = D_zeta * y_elem
    
    z_xi = D_xi * z_elem
    z_eta = D_eta * z_elem
    z_zeta = D_zeta * z_elem
    
    for i in 1:n_nodes_per_elem
        J = [x_xi[i] x_eta[i] x_zeta[i];
             y_xi[i] y_eta[i] y_zeta[i];
             z_xi[i] z_eta[i] z_zeta[i]]
        
        detJ = abs(det(J))
        w_distorted[i] = w_ref_elem[i] * detJ
    end
    
    vol_distorted = sum(w_distorted)
    println("Distorted Volume: $vol_distorted (expected 8.0)")
    @test abs(vol_distorted - 8.0) < 1e-12
    
    # Test integration of constant function on distorted mesh
    # Should be volume
    int_1_dist = sum(w_distorted)
    @test abs(int_1_dist - 8.0) < 1e-12
end
