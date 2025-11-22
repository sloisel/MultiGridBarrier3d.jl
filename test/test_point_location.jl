using MultiGridBarrier3d
using Test
using LinearAlgebra

@testset "Point Location Tests" begin
    println("Testing Point Location")
    
    k = 3
    L = 1
    
    # Test 1: Standard Cube
    geo_cube = fem3d(Float64; L=L, k=k)
    
    # Define a function u(x,y,z) = x + 2y + 3z
    u_cube = zeros(size(geo_cube.x, 1))
    for i in 1:size(geo_cube.x, 1)
        x, y, z = geo_cube.x[i, :]
        u_cube[i] = x + 2*y + 3*z
    end
    
    # Test point inside
    pt_in = [0.5, -0.2, 0.1]
    val, found = MultiGridBarrier3d.evaluate_field(geo_cube, u_cube, pt_in)
    @test found
    @test isapprox(val, pt_in[1] + 2*pt_in[2] + 3*pt_in[3], atol=1e-10)
    
    # Test point outside
    pt_out = [2.0, 0.0, 0.0]
    val, found = MultiGridBarrier3d.evaluate_field(geo_cube, u_cube, pt_out)
    @test !found
    @test val == 0.0
    
    # Test 2: Distorted Element
    # Use the same distorted element as in test_derivatives.jl
    # Map: x = xi + 0.5*eta, y = eta, z = zeta
    # Jacobian = [1 0.5 0; 0 1 0; 0 0 1]
    
    K_distorted = Float64[
        -1.5 -1.0 -1.0; # (-1,-1,-1) -> x = -1 - 0.5 = -1.5
         0.5 -1.0 -1.0; # ( 1,-1,-1) -> x =  1 - 0.5 =  0.5
        -0.5  1.0 -1.0; # (-1, 1,-1) -> x = -1 + 0.5 = -0.5
         1.5  1.0 -1.0; # ( 1, 1,-1) -> x =  1 + 0.5 =  1.5
        -1.5 -1.0  1.0; # ... z=1
         0.5 -1.0  1.0;
        -0.5  1.0  1.0;
         1.5  1.0  1.0
    ]
    
    geo_dist = fem3d(Float64; L=L, K=K_distorted, k=k)
    
    # Function u(x,y,z) = x + y + z
    u_dist = zeros(size(geo_dist.x, 1))
    for i in 1:size(geo_dist.x, 1)
        x, y, z = geo_dist.x[i, :]
        u_dist[i] = x + y + z
    end
    
    # Test point inside distorted element
    # Reference point xi = [0.5, 0.5, 0.5]
    # Physical point x = 0.5 + 0.5*0.5 = 0.75, y = 0.5, z = 0.5
    pt_dist_in = [0.75, 0.5, 0.5]
    
    val, found = MultiGridBarrier3d.evaluate_field(geo_dist, u_dist, pt_dist_in)
    @test found
    @test isapprox(val, pt_dist_in[1] + pt_dist_in[2] + pt_dist_in[3], atol=1e-10)
    
    # Test point that would be inside bounding box but outside element?
    # Bounding box of distorted element:
    # x: [-1.5, 1.5], y: [-1, 1], z: [-1, 1]
    # Point x=-1.4, y=1.0, z=0.0
    # At y=1, x range is [-0.5, 1.5]. So -1.4 is outside element but inside global bbox.
    
    pt_dist_out = [-1.4, 1.0, 0.0]
    val, found = MultiGridBarrier3d.evaluate_field(geo_dist, u_dist, pt_dist_out)
    @test !found
    
end
