using MultiGridBarrier3d
using PyCall
using Test

# Create dummy geometry
geo = fem3d(L=1, k=1)
u = rand(size(geo.x, 1))

println("Testing slice_orthogonal...")
# slice_orthogonal(x=0.2, y=0.7)
plot(geo, u; slice_orthogonal=Dict("x" => 0.2, "y" => 0.7))

println("Testing slice (normal='z')...")
# slice(normal='z')
plot(geo, u; slice=Dict("normal" => "z"))

println("Testing slice_along_axis...")
# slice_along_axis(n=3, axis='x')
plot(geo, u; slice_along_axis=Dict("n" => 3, "axis" => "x"))

println("Slicing tests passed (if no errors).")
