using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics
using FileIO

GLMakie.closeall()

GLMakie.closeall()

for testCase in 1:4
    if testCase == 1
        F, V = geosphere(2, 15.0)
        nSteps = 25
        w = 16
        xr,yr,zr = ntuple(_->range(-w, w,nSteps),3)
        pointSpacing = pointspacingmean(F,V)
        voxelSize = ((2.0*w)/(nSteps-1), (2.0*w)/(nSteps-1), (2.0*w)/(nSteps-1))
        origin = Point{3,Float64}(-w, -w, -w)
    elseif testCase == 2
        r = 1.0
        nc = 16
        t = range(2.0*π-(2.0*π/nc),0,nc)
        Vc = [Point3{Float64}(2.0+cos(tt),0.0,sin(tt)) for tt ∈ t]
        n = Vec{3, Float64}(0.0,0.0,1.0)
        num_steps = 24
        close_loop = true

        F,V = revolvecurve(Vc; extent = (2*pi - (2*pi/num_steps)), direction = :negative, 
        n = n, num_steps = num_steps, periodicity = (true,true), face_type = :tri)

        pointSpacing = pointspacingmean(F,V)
        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (r/4.0, r/4.0, r/4.0)
        xr = minV[1]:voxelSize[2]:maxV[1]
        yr = minV[2]:voxelSize[1]:maxV[2]
        zr = minV[3]:voxelSize[3]:maxV[3] 
        origin = minV
    elseif testCase == 3
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)

        pointSpacing = pointspacingmean(F,V)
        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (pointSpacing/3.0, pointSpacing/3.0, pointSpacing/3.0)
        xr = minV[1]:voxelSize[2]:maxV[1]
        yr = minV[2]:voxelSize[1]:maxV[2]
        zr = minV[3]:voxelSize[3]:maxV[3] 
        origin = minV
    elseif testCase == 4
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)

        F2, V2= geosphere(2, 20.0)
        V2 .-= Point{3, Float64}(0.0, 10.0, 20.0)
        invert_faces!(F2)
        F,V = joingeom(F,V,F2,V2)

        pointSpacing = pointspacingmean(F,V)
        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (pointSpacing/2.0, pointSpacing/2.0, pointSpacing/2.0)
        xr = minV[1]:voxelSize[2]:maxV[1]
        yr = minV[2]:voxelSize[1]:maxV[2]
        zr = minV[3]:voxelSize[3]:maxV[3] 
        origin = Point{3, Float64}(xr[1], yr[1], zr[1])
    elseif testCase == 5
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_ = mergevertices(F,V)

        F2, V2= geosphere(2, 20.0)
        V2 .-= Point{3, Float64}(0.0, 10.0, 20.0)
        invert_faces!(F2)
        F,V = joingeom(F,V,F2,V2)

        pointSpacing = pointspacingmean(F,V)
        minV = minp(V)
        maxV = maxp(V)
        voxelSize= (pointSpacing/2.0, pointSpacing/2.0, pointSpacing/2.0)
        xr = minV[1]:voxelSize[2]:maxV[1]
        yr = minV[2]:voxelSize[1]:maxV[2]
        zr = -25:voxelSize[3]:50  
        origin = Point{3, Float64}(xr[1], yr[1], zr[1])
    end

    B = mesh2bool(F, V, xr, yr, zr; tolEps = eps(Float64))

    FB, VB, _ = image2voxelmesh(B, B; meshType=:boundaryfaces, voxelSize=voxelSize, origin=origin)
    FBs, VBs = separate_vertices(FB,VB)

    fig = Figure(size=(1200,800))   
    ax1 = AxisGeom(fig[1,1])
    hp1 = meshplot!(ax1, F, V, color=(:green, 0.5), strokewidth=1.0, transparency=true)
    # normalplot(ax1, F, V)
    ax2 = AxisGeom(fig[1,2])
    hp2 = meshplot!(ax2, F, V, color=(:green, 0.5), strokewidth=0.0, transparency=true)
    hp3 = meshplot!(ax2, FBs, VBs, color=(:white, 0.5), strokewidth=1.0, transparency=true)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end