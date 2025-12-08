using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics
using Roots

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1
        nSteps = 25
        w = 2.0
        xr,yr,zr = ntuple(_->range(-w, w,nSteps),3)
        I = [sqrt(x^2 + y^2 + z^2) for x in xr, y in yr, z in zr]    
        d = 2.0
        voxelSelection = [i<=d for i in I] # Bool array
    elseif testCase == 2
        nSteps = 75
        w = 4.77
        xr,yr,zr = ntuple(_->range(-w, w,nSteps),3)
        function thisImage(x,y,z)
            phi=(1+sqrt(5))/2
            return 1.0/6.0 *(2.0 - (cos(x + phi*y) + cos(x - phi*y) + cos(y + phi*z) + cos(y - phi*z) + cos(z - phi*x) + cos(z + phi*x)))        
        end
        I = [thisImage(x,y,z) for x in xr, y in yr, z in zr]    
        B = [i<0.0 for i in I]
        voxelSelection = findall(B) # Cartesian indices
   end

    meshType = :boundaryfaces
    voxelSize = ((2.0*w)/(nSteps-1), (2.0*w)/(nSteps-1), (2.0*w)/(nSteps-1))
    origin = Point{3,Float64}(-w, -w, -w)

    F, V, C = image2voxelmesh(I, voxelSelection; meshType=meshType, voxelSize=voxelSize, origin=origin)


    k_pb = 0.1
    N = 50

    f_fun(k, λ, μ) = (1.0 - k*λ) * (1.0 - k*μ)
    function fun_check1(λ)
        μ = 1.0 /(k_pb-(1.0/λ))
        f1 = f_fun(1.0, λ, μ)
        f2 = f_fun(2.0, λ, μ)
        f = f_fun(k_pb, λ, μ)
        v = abs(f1+f2)
        if f>1.0
            v += abs(f-1.0)*1e9
        end
        return v  
    end

    λ = find_zero(fun_check1,0.5)
    μ = 1.0 /(k_pb-(1.0/λ))

    k = k_pb
    f = f_fun(k, λ, μ)
    fN = f^N

    println("λ = ", λ)
    println("μ = ", μ)
    println("f = ", f)
    println("f^N = ", fN)

    V = smoothmesh_taubin(F, V, N, λ, μ)

    # Visualization
    Fs, Vs = separate_vertices(F,V)
    Cs = simplex2vertexdata(Fs,C)

    fig = Figure(size=(800,800))   
    ax1 = AxisGeom(fig[1,1]; limits=(-w-voxelSize[2],w+voxelSize[2],-w-voxelSize[1],w+voxelSize[1],-w-voxelSize[3],w+voxelSize[3]))
    hp1 = meshplot!(ax1, Fs, Vs, color=Cs, strokewidth=1.0, colormap=Reverse(:Spectral))
    # scatter!(ax1,V, markersize=10, color=:red)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end