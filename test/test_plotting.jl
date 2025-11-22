using MultiGridBarrier3d
using Test
using FileIO

@testset "Plotting Tests" begin
    # Create a simple geometry
    k = 2
    L = 2
    geo = fem3d(Float64; L=L, k=k)
    
    # Create a dummy solution vector
    # u(x,y,z) = x^2 + y^2 + z^2sin(pi*y) * sin(pi*z)
    u = zeros(size(geo.x, 1))
    for i in 1:size(geo.x, 1)
        x, y, z = geo.x[i, :]
        u[i] = sin(pi * x) * sin(pi * y) * sin(pi * z)
    end

    # Plot
    println("Generating plot...")
    fig = plot_solution(geo, u; 
                        volume=true, 
                        isosurfaces=5, 
                        slices=:x, 
                        resolution=(800, 600),
                        show_grid=true)

    @test fig isa MGB3DFigure
    println("Figure type: ", typeof(fig))

    # Save figure
    filename = "test_plot.png"
    if isfile(filename)
        rm(filename)
    end
    
    save_figure(fig, filename)
    
    @test isfile(filename)
    println("Saved test_plot.png")
    
    # Cleanup
    if isfile(filename)
        rm(filename)
    end
end
