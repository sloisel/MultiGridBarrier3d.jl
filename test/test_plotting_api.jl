using MultiGridBarrier3d
using MultiGridBarrier
using PyPlot
using Test

@testset "Plotting API Tests" begin
    println("Testing Plotting API...")
    
    # Create dummy geometry and solution
    geo = fem3d(L=1, k=1)
    n_nodes = size(geo.x, 1)
    u = rand(n_nodes)
    
    # Test plot(geo, u)
    println("Testing plot(geo, u)...")
    fig = plot(geo, u; volume=true, show_grid=false)
    @test fig isa MultiGridBarrier3d.Plotting.MGB3DFigure
    
    # Test savefig
    println("Testing savefig...")
    filename = "test_plot.png"
    savefig(fig, filename)
    @test isfile(filename)
    rm(filename)
    
    # Test plot(sol) wrapper
    println("Testing plot(sol)...")
    # Mock AMGBSOL
    # struct AMGBSOL{T,X,W,Mat,Discretization}
    #    z::Matrix{T} ...
    # end
    # We can't easily mock it without importing internal types or running the solver.
    # Let's run the solver quickly.
    
    sol = fem3d_solve(L=1, k=1, maxiter=1)
    fig_sol = plot(sol)
    @test fig_sol isa MultiGridBarrier3d.Plotting.MGB3DFigure
    
    println("Plotting API tests passed!")
end
