# MultiGridBarrier3d.jl

A Julia package for solving 3D PDEs using the Spectral Barrier Method with hexahedral finite elements.

## Demo

Here is a simple example of solving a 3D problem and plotting the result.

```@example 1
using MultiGridBarrier3d

# Solve a 1-Laplacian problem on [-1,1]^3
sol = fem3d_solve(L=2, verbose=false)

# Plot the solution using PyVista
fig = plot(sol, show_grid=true, isosurfaces=[1.1,1.2,1.5])
savefig(fig, "fem3d_demo.png"); nothing # hide
```

![](fem3d_demo.png)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sloisel/MultiGridBarrier3d.jl")
```
