using MultiGridBarrier3d
using Test

@testset "FEM3D Struct Tests" begin
    k = 2
    L = 2
    # Test default creation
    geo = fem3d(Float64; L=L, k=k)
    
    @test geo.discretization isa FEM3D{k, Float64}
    @test geo.discretization.L == L
    
    # Check K is the standard cube
    K_expected, _ = MultiGridBarrier3d.cube_mesh(1)
    @test geo.discretization.K ≈ K_expected
    @test size(geo.discretization.K, 2) == 3
    
    # Check promote_to_Qk
    # promote_to_Qk(K, k) should match cube_mesh(k) for a single cube
    
    K_q1 = geo.discretization.K
    x_Qk = MultiGridBarrier3d.promote_to_Qk(K_q1, k)
    
    x_cube, _ = MultiGridBarrier3d.cube_mesh(k)
    
    @test size(x_Qk) == size(x_cube)
    @test x_Qk ≈ x_cube
end
