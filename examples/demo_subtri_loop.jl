using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0 #radius
F,V = platonicsolid(4,r)

# Fn,Vn = subTri(F,V,1; method=:loop)

################# 

fig = Figure(size=(800,800))

stepRange = 0:1:5
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 0,linewidth=30)

titleString = lift(hSlider.value) do stepIndex
    "n = " * string(stepIndex) * " refinement steps"
end

Mn = lift(hSlider.value) do stepIndex
    Fn,Vn = subtri(F,V,stepIndex; method=:Loop) 
    return GeometryBasics.Mesh(Vn,Fn)
end

ax = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
slidercontrol(hSlider,ax)

hp1=wireframe!(ax,GeometryBasics.Mesh(V,F),linewidth=3,color=:red, overdraw=false)
hp2=poly!(ax,Mn,strokewidth=1,color=:white, shading = FastShading)

# hp2=poly!(scene,M,strokewidth=1,color=:white, shading = FastShading)

fig