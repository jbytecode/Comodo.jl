using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

#=
This demo is for the `kabsch_rot` function. The Kabsch method allows one to 
determine the rotation that occurred between two meshes (or point sets) with 
point-to-point correspondence. In this demo a triangulated surface is loaded 
from an STL file. Next the coordinates are rotated, and the Kabsch method is 
used to retrieve the rotation performed, allowing one to "unrotate" the data.  
=#

# Loading a mesh

fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = tofaces(faces(M))
V1 = topoints(coordinates(M))
F,V1,_,_ = mergevertices(F,V1)

# Define a rotation tensor using Euler angles
Q = RotXYZ(0.0,0.25*π,0.25*π)
V2 = [GeometryBasics.Point{3, Float64}(Q*v) for v ∈ V1] 

R = kabsch_rot(V2,V1)
V3 = [GeometryBasics.Point{3, Float64}(R*v) for v ∈ V2] 

# Rotate the coordinates
fig = Figure(size = (800,800))
ax = Axis3(fig[1, 1], aspect = :data)

hp1 = poly!(ax, GeometryBasics.Mesh(V1,F), color=:green,transparency=false,shading = FastShading)
hp2 = poly!(ax, GeometryBasics.Mesh(V2,F), strokewidth=2,color=:red,shading = FastShading)
hp3 = wireframe!(ax, GeometryBasics.Mesh(V3,F), color=:red,linewidth=2)

Legend(fig[1, 2],[hp1,hp2,hp3],["Initial","Rotated","Back rotated"])

fig


