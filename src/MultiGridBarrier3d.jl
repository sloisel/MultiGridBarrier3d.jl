module MultiGridBarrier3d

using LinearAlgebra
using SparseArrays
using MultiGridBarrier

include("Geometry.jl")

include("MeshGen.jl")
include("Operators.jl")
include("Plotting.jl")
using .Plotting

export FEM3D, plot, savefig, fem3d, fem3d_solve

"""
    fem3d_solve(::Type{T}=Float64; 
                D = [:u :id; :u :dx; :u :dy; :u :dz; :s :id],
                f = (x) -> T[0.5, 0.0, 0.0, 0.0, 1.0],
                g = (x) -> T[x[1]^2 + x[2]^2 + x[3]^2, 100.0],
                rest...) where {T}

Solve a 3D PDE using the Spectral Barrier Method.

# Arguments
- `T`: Numeric type (default `Float64`).
- `D`: Operator structure definition.
- `f`: Source term function. Must return a Vector of length matching `D` (default 5).
- `g`: Boundary condition function.
- `rest...`: Additional arguments passed to `fem3d` and `amgb` (e.g., `L`, `k`, `maxiter`).

# Returns
- `AMGBSOL`: Solution object containing the result `z` and convergence history.
"""
function fem3d_solve(::Type{T}=Float64; 
    D = [:u :id; :u :dx; :u :dy; :u :dz; :s :id],
    f = (x) -> T[0.5, 0.0, 0.0, 0.0, 1.0],
    g = (x) -> T[x[1]^2 + x[2]^2 + x[3]^2, 100.0],
    rest...) where {T}
    
    # Create geometry
    geo = fem3d(T; rest...)
    
    # Call amgb
    return amgb(geo; D=D,f=f,g=g,rest...)
end

end
