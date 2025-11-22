# MultiGridBarrier3d.jl

A Julia package for solving 3D PDEs using the Spectral Barrier Method with hexahedral finite elements.

## Demo

Here is a simple example of solving a 3D problem and plotting the result.

```@example 1
using MultiGridBarrier3d

# Solve a simple problem on [-1,1]^3
# -Delta u = f, u = 0 on boundary
# Default f and g are used.
sol = fem3d_solve(L=2, k=1, verbose=false)

# Plot the solution
# This uses PyVista backend but extends PyPlot API
fig = plot(sol, volume=(;), show_grid=true)
savefig(fig, "fem3d_demo.png"); nothing # hide
```

![](fem3d_demo.png)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sloisel/MultiGridBarrier3d.jl")
```
