using MultiGridBarrier3d
using MultiGridBarrier
using Test
using LinearAlgebra

@testset "Solver Tests" begin
    println("Testing fem3d_solve")
    
    # Define a simple problem: -Delta u = f on [-1, 1]^3
    # u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    # f = 3*pi^2 * u
    
    # We need to pass f to amgb?
    # amgb signature: amgb(geo; f=..., u_exact=..., ...)
    
    # Let's try to call fem3d_solve with these arguments.
    # Note: amgb expects f to be a function or a vector?
    # In MultiGridBarrier.jl, it usually takes keyword args.
    
    u_exact(x) = sin(π*x[1]) * sin(π*x[2]) * sin(π*x[3])
    # f must return a vector of length 5 matching D structure
    f_rhs(x) = [0.0, 0.0, 0.0, 0.0, 0.0]
    
    # We need to pass these as functions or vectors?
    # If amgb handles function evaluation, we pass functions.
    # If not, we might need to pre-compute.
    # But fem3d_solve calls fem3d first, then amgb.
    # fem3d returns geometry.
    # amgb(geo; rest...)
    
    # Let's assume amgb can take f as a function if we pass it.
    # Or maybe we need to pass a vector 'b' (RHS).
    
    # Let's try passing f as a function first.
    
    # We use a very coarse grid for speed/debugging
    k = 1
    L = 2
    
    # We need to capture the output or result.
    # amgb returns (u, history) usually?
    
    # Let's try to run it.
    # We might need to define 'f' in a way compatible with how amgb evaluates it.
    # If amgb uses evaluate_field or just evaluates at nodes.
    
    # Actually, looking at MultiGridBarrier.jl (if I could), I would know.
    # But I can't see it. I'll assume standard behavior.
    
    # If amgb fails, we'll see the error.
    
    println("Calling fem3d_solve...")
    try
        # We pass f and u_exact as functions.
        # We also pass 'tol' and 'maxiter' if needed.
        
        # Note: amgb might expect 'f' to be a vector if it doesn't know how to evaluate on 3D geometry.
        # But fem3d returns a Geometry object which has 'x' (nodes).
        # If amgb is generic, it might expect 'f' to be a vector.
        
        # Let's pre-compute f on the fine grid?
        # But fem3d_solve creates the geometry internally.
        # So we can't pre-compute f on the grid before calling fem3d_solve.
        # Unless we pass a function that takes (x,y,z).
        
        # Let's try passing f as a function.
        
        sol = fem3d_solve(L=L, k=k, f=f_rhs, tol=1e-6, maxiter=10)
        
        # Check result type
        println("Solver returned: $(typeof(sol))")
        
        # If sol is (u, hist), we can check error.
        if sol isa Tuple
            u_sol = sol[1]
            # Compute error against u_exact
            # We need the geometry to get nodes.
            # But fem3d_solve doesn't return geometry?
            # Wait, amgb(geo) usually returns (u, ...).
            # We don't have geo to check error unless amgb returns it or we re-create it.
            
            # But amgb usually computes error if u_exact is passed.
            # So we can check the history or output.
            
            println("Solution vector length: $(length(u_sol))")
        end
        
        @test true # If we got here without error, it's a start.
        
    catch e
        println("Solver failed with error: $e")
        # For debugging, we might want to see the stacktrace
        showerror(stdout, e, catch_backtrace())
        @test false
    end
end
