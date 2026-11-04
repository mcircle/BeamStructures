using BeamStructures
using ForwardDiff
using SciMLSensitivity
using Zygote
using LinearAlgebra
using Statistics
using GLMakie
using NonlinearSolve
using BenchmarkTools
using Optimisers
using Printf
import ChainRulesCore as CRC

clampmarker = Makie.BezierPath([
    MoveTo(Point2f(0,0)),
    LineTo(Point2f(-1,0)),
    LineTo(Point2f(-1,0.15)),

    LineTo(Point2f(-0.7,0.15)),
    LineTo(Point2f(-0.9,0.9)),
    LineTo(Point2f(-0.7,0.9)),
    LineTo(Point2f(-0.5,0.15)),

    LineTo(Point2f(-0.3,0.15)),
    LineTo(Point2f(-0.5,0.9)),
    LineTo(Point2f(-0.3,0.9)),
    LineTo(Point2f(-0.1,0.15)),
    
    LineTo(Point2f(0.1,0.15)),
    LineTo(Point2f(-0.1,0.9)),
    LineTo(Point2f(0.1,0.9)),
    LineTo(Point2f(0.3,0.15)),

    LineTo(Point2f(0.5,0.15)),
    LineTo(Point2f(0.3,0.9)),
    LineTo(Point2f(0.5,0.9)),
    LineTo(Point2f(0.7,0.15)),

    LineTo(Point2f(0,0.15)),
    LineTo(Point2f(1,0.15)),
    LineTo(Point2f(1,0)),
    ClosePath(),
])

slider_marker =Makie.BezierPath(Makie.Polygon(
    0.8f0 .*Point2f[(-1, -0.5), (-1, -0.1), (-1.5, -0.1), (-1.5, 0.1), (-1,0.1), (-1,0.5), (1,0.5), (1,0.1), (1.5,0.1), (1.5,-0.1), (1,-0.1), (1,-0.5), (-1,-0.5)],
    [0.8f0 .* Point2f[(-1,-0.3), (-1, -0.1), (1,-0.1), (1, -0.3)], 0.8f0 .* Point2f[(-1,0.3), (-1, 0.1), (1,0.1), (1, 0.3)]]
))
Makie.DEFAULT_MARKER_MAP[:clamp] = clampmarker
Makie.DEFAULT_MARKER_MAP[:slider] = slider_marker

w1 = 30f-1
b = 10f-1
t = 1.6f-1
λ = 0.5f0
β = π/4

w2 = w1 * (1-λ)/λ

h = Float32((w1 + w2)/(2 * tan(β)))
L = Float32(w1/(2λ * sin(β))  * 1.2)


Fx = 00f0 #5.5f0
Fy = 0f0 #-10
Mz = 0f0 #2.4
E = 1.4f+5/(1-0.4f0^2)
no1 = Clamp(zeros(Float32,6))
no2 = Clamp{Float32}(w1 ,0,zeros(Float32,4)...)
no3 = Branch{Float32}(w1,h ,zeros(Float32,4)...)
no4 = Branch{Float32}(0 ,h,zeros(Float32,4)...)
# no5 = Branch{Float32}(w2/2,h,0f0,Fx,Fy,Mz)
no5 = Clamp{Float32}(w2/2,h,0.f0,Fx,Fy,Mz)

kappa =  [0f0,-1f0,0f0,1f0,0f0] 
beam1 = CurvedBeam{Float32}( L,t,b,kappa,θs = β;E = E,θe = π/2)
beam2 = CurvedBeam{Float32}( L,t,b,-kappa,θs = π/2 +  β;E= E,θe = π/2)

# beam1 = Beam{Float32}(L,t,b,0f0,θs = β;E = E)
# beam2 = Beam{Float32}(L,t,b,0f0,θs = π/2 +  β;E= E)
# beam3 = CurvedBeam{Float32}(w2/2,10t,b,kappa,θs = deg2rad(180);E=10f+9,θe = π/2)
beam3 = Beam{Float32}(w2/2,10t,b,0f0,θs = deg2rad(180);E=10f+9)
beam4 = Beam{Float32}(w2/2,10t,b,0f0,θs = deg2rad(0);E=10f+9)

adj = [ 0 0 1 0 0;# 0 0; #1 Clamp
        0 0 0 1 0;#0 0; #2 Branch 
        1 0 0 0 1;#1 0; #3 
        0 1 0 0 1;#1 0; #4 
        0 0 1 1 0]

str = Structure(adj)

beams,nodes = BeamStructures.prepare(no1,no2,beam1,no3,beam2,beam3,no4,beam4,no5)

sols,beams,nodes_ = str(zeros(Float32,3,7),beams,nodes,false)
r = zeros(Float32,3,6)
res = BeamStructures.residuals!(r,str,sols,beams,nodes)
ini = zeros(Float32,3,6)
sols = str(ini,beams,nodes,true)

GLMakie.activate!()
fig =Figure();
ax = Makie.Axis(fig[1,1])
for (n,sol) in enumerate(sols)
    l = beams[n].l
    lines!(ax,sol[3,:].*l,sol[4,:].*l,label = "Beam $(n)")
end
axislegend(ax)
display(fig)

function snapthroughNon(str::Structure,x,beams,nodes)
    r = zero(x)
    sols,beams,nodes_ = str(x,beams,nodes)
    t = BeamStructures.residuals!(r,str,sols,beams,nodes_)
    return t
end
snapthroughNon(str,ini,(beams,nodes)) = snapthroughNon(str,ini,beams,nodes)


angles = range(-1f0,1f0,length=9)


function movetop(str::Structure,inits::NTuple{N,<:AbstractMatrix{T}},beams,nodes) where{N,T}
    l = zero(T)
    for  (ins,a) in zip(inits,angles)
        x =  h * cos(π/2f0 + a) ./2
        y =  h/2 * (1 -sin(π/2f0 + a) )
        n5 = nodes.Node_5 + (;x = x,y = -y,ϕ = a)
        nodes_ = (;nodes..., Node_5 = n5)
        r = zero(ins)
        sols,beams,nodes__ = str(ins,beams,nodes_)
        t = BeamStructures.residuals!(r,str,sols,beams,nodes__)
        l += l + sum(abs2, t)
        
    end
    return l
end


function movetopForce(str::Structure,inits::AbstractArray{T,N},beams,nodes) where{N,T}
    l = zero(T)
    for (i,a) in enumerate(angles)
        x =  h .* cos.(π/2f0 + a) ./2
        y =  h/2 .* (1 .-sin.(π/2f0 + a) )
        n5 = nodes.Node_5 + (;x = x,y = -y,ϕ = a)
        nodes_ =(;nodes..., Node_5 = n5)
        r = zeros(T,3,size(inits,2))
        sols,beams,nodes__ = str(inits[:,:,i],beams,nodes_)
        t = BeamStructures.residuals!(r,str,sols,beams,nodes__)
        n = Zygote.@ignore BeamStructures.normfactor_f(beams[end-1])
        n2 = Zygote.@ignore BeamStructures.normfactor_f(beams[end])
                                             #  -fx * sin(a) + fy * cos(a)
                                             
        l += 10 * mean(abs2, t) + abs2(-sols[6,2,3:4]' * [1/n,1/n2] * sin(a) + sols[5,2,3:4]' * [1/n,1/n2] * cos(a)) / length(angles)
    end
        
    return l
end


function movetopnl(str::Structure,inits::AbstractArray{T,N},beams,nodes) where{N,T}
    l = zero(inits)
    for (i,a) in enumerate(angles)
        x =  h * cos(π/2f0 + a) /2
        y =  h/2 * (1 -sin(π/2f0 + a) )
        n5 = nodes.Node_5 + (;x = x,y = -y,ϕ = a)
        nodes_ = (;nodes..., Node_5 = n5)
        r = zeros(T,3,size(inits,2))
        sols,beams,nodes__ = str(inits[:,:,i],beams,nodes_)
        t = BeamStructures.residuals!(r,str,sols,beams,nodes__)
        l[:,:,i] +=t
    end
        
    return l
end
movetopnl(str,inits,(beams,nodes)) = movetopnl(str,inits,beams,nodes)

inits1_ = randn(Float32,3,6)
snapthroughNon(str,inits1_,beams,nodes)

ini4move = zeros(Float32,3,6,9)
movetop(str,ini4move,beams,nodes)

using NonlinearSolve

nlfn = NonlinearFunction((x,p)->snapthroughNon(str,x,p))

nlprb = NonlinearLeastSquaresProblem(nlfn,inits1_,p = (beams,nodes))
nlsolf = solve(nlprb,TrustRegion(;autodiff = AutoForwardDiff()))
nlsolz = solve(nlprb,TrustRegion(;autodiff = AutoZygote()))

randini = 5randn(Float32,3,6)
randini[:,1:2] ./= 20
sol = str(randini,beams,nodes,true)
GLMakie.activate!()
# CairoMakie.activate!()
figdyn = Figure();
axdyn = Makie.Axis(figdyn[1,1],xlabel = L"$x$ in [mm]",ylabel = L"$y$ in [mm]")
nodesobs = Observable([Point2f(n.x,n.y) for n in nodes])
nodesobs[][3] += randini[1:2,1]
nodesobs[][4] += randini[1:2,2]

draw1 = Observable(Point2f.(sol.u[1](0:0.01:1;idxs = 3:4).u .* beams[1].l))
draw2 = Observable(Point2f.(sol.u[2](0:0.01:1;idxs = 3:4).u .* beams[2].l))
draw3 = Observable(Point2f.(sol.u[3](0:0.01:1;idxs = 3:4).u .* beams[3].l))
draw4 = Observable(Point2f.(sol.u[4](0:0.01:1;idxs = 3:4).u .* beams[4].l))
heights = [Observable(10 * b.h) for b in beams] 
lines!(axdyn,draw1,linewidth = heights[1],color = :blue, label = "Beam")
lines!(axdyn,draw2,linewidth = heights[2],label = "Beam",color = :blue,)
lines!(axdyn,draw3,linewidth = heights[3],label = "Beam",color = :blue,)
lines!(axdyn,draw4,linewidth = heights[4],label = "Beam",color = :blue,)
scatter!(axdyn,nodesobs,marker = [:clamp,:clamp,:circle,:circle,:circle],markersize = 10,color = :blue)
# ax2dyn = Makie.Axis(figdyn[1,2],xlabel = "Iteration",ylabel = "Loss")
# lossobs = Observable(Float32[])
# lines!(ax2dyn,lossobs,label = "Loss")
display(figdyn)
save("initial_solution.svg",figdyn)

beams_ = beams

nodes_ = nodes
# inits1_ = zeros(Float32,3,7)
inits = Tuple([randn(Float32,3,6) for _ in 1:9])


optpars_state = Optimisers.setup(Adam(0.02,(0.9,0.99)),inits)
optbeam_state = Optimisers.setup(Adam(0.0001,(0.9,0.99)),beams)
optnode_state = Optimisers.setup(Adam(0.0001,(0.9,0.99)),nodes)
Optimisers.adjust!(optpars_state,0.01)
Optimisers.adjust!(optbeam_state,0.01)
Optimisers.adjust!(optnode_state,0.01)

display(figdyn)
for n in 1:5000

    l,grads = withgradient((x) -> movetop(str,x,beams,nodes),inits)
    optpars_state,inits = Optimisers.update(optpars_state,inits,grads[1])
    # optbns_state,bnstr_ = Optimisers.update(optbns_state,bnstr_,grads[2])

    @printf "iteration: %i \t Loss: %.6f \n" n l
    if n % 1 == 0
        odesol = str(inits[5],beams,nodes,true)
        off_points = zeros(Float32,2,5)
        off_points[:,3:end-1] .= inits[5][1:2,1:2]
        for (i,b) in enumerate(beams)
            heights[i][] = 10 * b.h
        end
        nodesobs[] = [Point2f(n.x + op[1],n.y + op[2]) for (n,op) in zip(nodes,eachcol(off_points))]
        draw1[] = Point2f.(odesol.u[1](0:0.01:1;idxs = 3:4).u .* beams[1].l)
        draw2[] = Point2f.(odesol.u[2](0:0.01:1;idxs = 3:4).u .* beams[2].l)
        draw3[] = Point2f.(odesol.u[3](0:0.01:1;idxs = 3:4).u .* beams[3].l)
        draw4[] = Point2f.(odesol.u[4](0:0.01:1;idxs = 3:4).u .* beams[4].l)
    end 
end

inits
nlfunc = NonlinearFunction((x,p)->movetopnl(str,x,p...))
nlprob = NonlinearLeastSquaresProblem(nlfunc,cat(inits...,dims =3),p = (beams,nodes))
snl = solve(nlprob,TrustRegion(;autodiff = AutoForwardDiff()))
snl.u[:,:,5]
beams_[1]
sum(snl.u[:,end-1:end,:],dims = 2)
sum(inits[:,end-1:end,:],dims = 2)

using CairoMakie
CairoMakie.activate!()

GLMakie.activate!()
fig =Figure();
ax = Makie.Axis(fig[1,1],xlabel = L"$x$ in [mm]",ylabel = L"$y$ in [mm]")
solinit = str(zeros(Float32,3,6),beams,nodes,true)

for (n,sol) in enumerate(solinit)
    l = beams[n].l
    height = 10 * beams[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :gray,alpha = 0.3,label = "Initial")
end
display(fig)
nodes_ = [Point2f(n.x,n.y) for n in nodes]
scatter!(ax,nodes_[1:2],marker = [:clamp,:clamp],markersize = 10,color = :black)
scatter!(ax,nodes_[3:end],marker = [:circle,:circle,:circle],markersize = 10,color = :grey,alpha = 0.5)
lines!(ax,Makie.Circle(Point2f(1.5,1.5),1.5),linestyle = :dash,color = :black,label = "path")

display(fig)
inits[:,:,5]
for (ind,init) in  enumerate(inits[[1,5,9]])
    sols = str(init,beams,nodes,true)
    for (n,sol) in enumerate(sols)
        l = beams[n].l
        height = 10 * beams[n].h
        lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :blue,alpha = 0.5,label = "Adam")
    end 
end 
display(fig)

# lossobs = Observable(Float32[])
# ax2 = Makie.Axis(fig[1,2],xlabel = "Iteration",ylabel = "Loss")
# lines!(ax2,lossobs,label = "Loss")
sols = str(snl.u[:,:,5],beams,nodes,true)
for (n,sol) in enumerate(sols)
    l = beams[n].l
    height = 10 * beams[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :blue,alpha = 1,label = "Trust Region")
end 
for (ind,init) in  enumerate(eachslice(snl.u[:,:,[1,9]],dims = 3))
    sols = str(init,beams,nodes,true)
    for (n,sol) in enumerate(sols)
        l = beams[n].l
        height = 10 * beams[n].h
        lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :red,alpha = 0.5,label = "Trust Region")
    end 
end 
scatter!(ax,[1.5],[1.5],marker = :cross,markersize = 10,color = :red,label = "Pivot")
axislegend(ax, merge = true, unique = true)
display(fig)

save("initial_trust_region_curved.svg",fig)

using CairoMakie
CairoMakie.activate!()
save("adam_vs_trustregion.svg",fig)

fig = Figure();
ax = Makie.Axis(fig[1,1],xlabel = L"$x$ in [mm]",ylabel = L"$y$ in [mm]")
sols = str(snl.u[:,:,5],beams,nodes,true)
    for (n,sol) in enumerate(sols)
        l = beams[n].l
        height = 10 * beams[n].h
        lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :red,alpha = 1,label = "Trust Region")
    end 
display(fig)
save("initial_sol_curved.svg",fig)
#forces
sum(snl.u[2,end-1:end,5]) #fx
sum(snl.u[3,end-1:end,5]) #fy


beams_ = beams
nodes_ = nodes
inits_ = zeros(Float32,3,6,9)/10
inits_ = snl.u
optpars_state = Optimisers.setup(Adam(0.02,(0.9,0.99)),inits_)
optbeam_state = Optimisers.setup(Adam(0.0001,(0.9,0.99)),beams)
optnode_state = Optimisers.setup(Adam(0.0001,(0.9,0.99)),nodes)
Optimisers.adjust!(optpars_state,0.01)
Optimisers.adjust!(optnode_state,0.001)
Optimisers.adjust!(optbeam_state,0.001)
GLMakie.activate!()

ax2 = Makie.Axis(figdyn[1,2],xlabel = "angle in [°]",ylabel = "force in [N]")
f = Vector{Float32}(undef,9)
f2 = Observable(Vector{Float32}(undef,9))
n = BeamStructures.normfactor_f(beams[ end ])
n2 = BeamStructures.normfactor_f(beams_[ end ])
for (i,(init,innl)) in enumerate(zip(inits,eachslice(snl.u,dims = 3)))
    sols = str(init,beams_,nodes_,true)    
    sols2 = str(innl,beams_,nodes_,true)    
                        #  -fx * sin(a) + fy * cos(a)
    f2[][i] = (-sum(sols[5,2,3:4]) * sin(angles[i]) + sum(sols[6,2,3:4]) *cos(angles[i])) / n
    f[i] = (-sum(sols2[5,2,3:4]) * sin(angles[i]) + sum(sols2[6,2,3:4]) *cos(angles[i])) /n2
end 

lines!(ax2,rad2deg.(angles),f,label = "without TO",color = :red)
lines!(ax2,rad2deg.(angles),f2,label = "TO",color = :blue)
axislegend(ax2)

inits_ = snl.u
display(figdyn)

# stream = VideoStream(figdyn)
lossval = Vector{Float32}()
for n in 1:5000
    l,grads = withgradient((x,y,z) -> movetopForce(str,x,y,z),inits_,beams_,nodes_)
    optpars_state,inits_ = Optimisers.update(optpars_state,inits_,grads[1])
    optbeam_state,beams_ = Optimisers.update(optbeam_state,beams_,grads[2])
    optnode_state,nodes_ = Optimisers.update(optnode_state,nodes_,grads[3])

    @printf "iteration: %i \t Loss: %.6f \n" n l
    push!(lossval,l)
    # push!(lossfun,l)
    # notify(lossobs)
    if n % 1 == 0
        odesol = str(inits_[:,:,5],beams_,nodes_,true)
        off_points = zeros(Float32,2,5)
        off_points[:,3:end-1] .= inits_[1:2,1:2,5]
        for (i,b) in enumerate(beams_)
            heights[i][] = 10 * b.h
        end
        nodesobs[] = [Point2f(n.x + op[1],n.y + op[2]) for (n,op) in zip(nodes_,eachcol(off_points))]
        draw1[] = Point2f.(odesol.u[1](0:0.01:1;idxs = 3:4).u .* beams_[1].l)
        draw2[] = Point2f.(odesol.u[2](0:0.01:1;idxs = 3:4).u .* beams_[2].l)
        draw3[] = Point2f.(odesol.u[3](0:0.01:1;idxs = 3:4).u .* beams_[3].l)
        draw4[] = Point2f.(odesol.u[4](0:0.01:1;idxs = 3:4).u .* beams_[4].l)
        # n = BeamStructures.normfactor_f(beams_[ end ])
        # for (i,init) in enumerate(eachslice(inits_,dims = 3))
        #     # sols,_,_ = str(init,beams_,nodes_,false)    
        #                         #  -fx * sin(a) + fy * cos(a)
        #     f2[][i] = sum(-init[2,end-1:end] .* sin(angles[i]) .+ init[3,end-1:end] .* cos(angles[i]))
        # end 
        # notify(f2)
        autolimits!(axdyn)
        # autolimits!(ax2)
        # recordframe!(stream)
    end 
end

using JLD2

# @save "Flexible_Pivot_Data/optimisation_curved_zero_26_07_29.jld2" inits_ beams_ nodes_ lossval
@load "Flexible_Pivot_Data/optimisation_curved_zero_flexible_pivot.jld2" inits_ beams_ nodes_ lossval

fig = Figure();
ax = Makie.Axis(fig[1,1],xlabel = "Iteration",ylabel = "Loss",yscale = log10)
lines!(lossval,label = "Loss",color = :blue)
display(fig)
using CairoMakie
CairoMakie.activate!()
save("loss_curve_curved.eps",fig)
GLMakie.activate!()

# save("optimisation_curved_zero.mp4",stream)

sol = str(snl.u[:,:,5],beams,nodes,true)
draw1_ = Observable(Point2f.(sol[1](0:0.01:1;idxs = 3:4).u .* beams[1].l))
draw2_ = Observable(Point2f.(sol[2](0:0.01:1;idxs = 3:4).u .* beams[2].l))
draw3_ = Observable(Point2f.(sol[3](0:0.01:1;idxs = 3:4).u .* beams[3].l))
draw4_ = Observable(Point2f.(sol[4](0:0.01:1;idxs = 3:4).u .* beams[4].l))
heights_ = [Observable(10 * b.h) for b in beams] 
lines!(axdyn,draw1_,linewidth = heights_[1],color = :grey,alpha = 0.5, label = "initial")
lines!(axdyn,draw2_,linewidth = heights_[2],label = "initial",color = :grey,alpha = 0.5)
lines!(axdyn,draw3_,linewidth = heights_[3],label = "initial",color = :grey,alpha = 0.5)
lines!(axdyn,draw4_,linewidth = heights_[4],label = "initial",color = :grey,alpha = 0.5)
lines!(axdyn,Makie.Circle(Point2f(1.5,1.5),1.5),linestyle = :dash,color = :black,label = " Initial path")
nodes_[end]
mid = Point2f(nodes_[end].x,nodes_[end].y -1.5)
lines!(axdyn,Makie.Circle(mid,1.5),linestyle = :dash,color = :red,label = "current path")
axislegend(axdyn, merge = true, unique = true)
display(figdyn)

CairoMakie.activate!()
save("optimisation_curved_optimized.eps",figdyn)

nlprob = NonlinearLeastSquaresProblem(nlfunc,inits_,p = (beams_,nodes_))
snl_ = solve(nlprob,TrustRegion(;autodiff = AutoForwardDiff()),maxiters = 10000)
CairoMakie.activate!()

fig =Figure();
ax = Makie.Axis(fig[1,1],xlabel = L"$x$ in [mm]",ylabel = L"$y$ in [mm]")
solinit = str(zeros(Float32,3,6),beams_,nodes_,true)

for (n,sol) in enumerate(solinit.u)
    l = beams[n].l
    height = 10 * beams[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :gray,label = "Init")
end
nodeso = [Point2f(n.x,n.y) for n in nodes_]
scatter!(ax,nodeso[1:2],marker = [:clamp,:clamp],markersize = 10,color = :grey)
solinit = str(inits_[:,:,5],beams_,nodes_,true)

for (n,sol) in enumerate(solinit.u)
    l = beams_[n].l
    height = 10 * beams_[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :red,label = "Optim")
end
nodeso = [Point2f(n.x,n.y) for n in nodes_]
scatter!(ax,nodeso[1:2],marker = [:clamp,:clamp],markersize = 10,color = :red)
scatter!(ax,[1.5],[1.5],marker = :cross,markersize = 10,color = :red,)
axislegend(ax, merge = true, unique = true)
display(fig)
nodes_[5]

CairoMakie.activate!()
save("init_vs_optim_curved.eps",fig)

fig = Figure();
ax = Makie.Axis(fig[1,1],xlabel = L"$x$ in [mm]",ylabel = L"$y$ in [mm]")
sols = str(inits_[:,:,5],beams_,nodes_,true)
    
for (n,sol) in enumerate(sols)
    l = beams[n].l
    height = 10 * beams[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :grey,alpha = 0.5,label = "Ruhe")
end 

for (ind,init) in  enumerate(eachslice(inits_[:,:,[1,5,9]],dims = 3))
    sols = str(init,beams_,nodes_,true)
    c = Cycled(ind)
    for (n,sol) in enumerate(sols)
        l = beams_[n].l
        height = 10 * beams_[n].h
        lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = c,label = "Adam")
    end 
end 
display(fig)

nodes_o = [Point2f(n.x,n.y) for n in nodes_]
scatter!(ax,nodes_o[1:2],marker = [:clamp,:clamp],markersize = 10,color = :black)
scatter!(ax,nodes_o[3:end],marker = [:circle,:circle,:circle],markersize = 10,color = :grey,alpha = 1)

# inits_ = snl.u
x = zeros(Float32,9)
y = zeros(Float32,9)
for (n,init) in enumerate(eachslice(inits_,dims = 3))
    sols,_ = str(init,beams_,nodes_)
    l = beams_[ end ].l
    x[n] = sols[3,2,end] * l
    y[n] = sols[4,2,end] * l
    
end
mid = (nodes_o[end][1],nodes_o[end][2] -1.5)
x .- mid[1]
y .-mid[2]
r = hypot.(x.-mid[1],y .- mid[2])
r_ = hypot.(x[5] .- mid[1],y[5] .- mid[2] )
(r.-1.5)./1.5 .*100

c = Circle(Point2f(mid...),1.5)
lines!(ax,c,linestyle = :dash,color = :black,label = "Zielpfad")
scatter!(ax,mid,marker = :cross,markersize = 10,color = :red,label = "Pivot")
display(fig)
save("optim_solutions_curved.eps",fig)

f_ = Vector{Float32}(undef,9)
f2_ = Vector{Float32}(undef,9)
n = BeamStructures.normfactor_f(beams[ end-1 ])
n2 = BeamStructures.normfactor_f(beams_[ end ])
for (i,(init,innl)) in enumerate(zip(eachslice(inits_,dims = 3),eachslice(snl.u,dims = 3)))
    sols = str(init,beams_,nodes_,true)    
    sols2 = str(innl,beams_,nodes_,true)    
                        #  -fx * sin(a) + fy * cos(a)
    f_[i] = sum(-sols[5,2,3:4] ./n .* sin(angles[i]) .+ sols[6,2,3:4] ./n .*cos(angles[i])) 
    f2_[i] = sum(-sols2[5,2,3:4] ./n2 .* sin(angles[i]) .+sols2[6,2,3:4] ./n2 .*cos(angles[i])) 
end 
f
f2
fig = Figure();
ax = Makie.Axis(fig[1,1],xlabel = L"\text{Winkel in [°]}",ylabel = L"$F_r$ in [N]")
lines!(ax,rad2deg.(angles),f2,label = "optimized",color = :red)
lines!(ax,rad2deg.(angles),f,label = "initial",color = :blue)
axislegend(ax)
display(fig)
save("force_comparison.eps",fig)


### ________________________________Prepare for Export to CAD________________________________

using DataFrames,CSV,JLD2 

@load "Flexible_Pivot_Data/optimisation_curved_zero_26_07_29.jld2" inits_ beams_ nodes_ lossval
beams_
nodes_
inits_
sol = str(zeros(Float32,beams_,nodes_),beams_,nodes_,true)
GLMakie.activate!()
fig =Figure();
ax = Makie.Axis(fig[1,1],xlabel = L"$x$ in [mm]",ylabel = L"$y$ in [mm]")

for (n,sol) in enumerate(sol.u)
    l = beams_[n].l
    height = 10 * beams_[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :gray,alpha = 0.3,label = "Initial")
end
display(fig)

_,_,nodes2_ = str(inits_[:,:,5],beams_,nodes_)
sol_ = str(inits_[:,:,5],beams_,nodes_,true)
for (n,sol) in enumerate(sol_.u)
    l = beams_[n].l
    height = 10 * beams_[n].h
    lines!(ax,sol,idxs = ((x,y)->(x,y).*l,3,4),linewidth = height,color = :red,alpha = 0.7,label = "Initial")
end

beams_[1].θe
beams_[2].θe
nodes2_[3].x 
nodes2_[3].y
nodes2_[4].x
nodes2_[4].y

nodes[5].y

beams_[1].θe + nodes2_[3].ϕ
beams_[2].θe + nodes2_[4].ϕ

beams_[1].w
beams_[1].h
beams_[2].w
beams_[2].h

mat = Vector{Float32}()

for n in 1:size(inits_,3)
    sol_ = str(inits_[:,:,n],beams_,nodes_,true)
    b1 = (sol_.u[1](1;idxs = [2,3,4]) .- sol.u[1](1;idxs = [2,3,4])) .* [1,beams_[1].l,beams_[1].l]
    b2 = (sol_.u[2](1;idxs = [2,3,4]) .- sol.u[2](1;idxs = [2,3,4])) .* [1,beams_[2].l,beams_[2].l]
    append!(mat, b1[[2,3,1]])
    append!(mat, b2[[2,3,1]])
end 
mat = reshape(mat,6,:)

df = DataFrame(x1 = mat[1,:],y1 = mat[2,:],phi1 = mat[3,:],x2 = mat[4,:],y2 = mat[5,:],phi2 = mat[6,:])
CSV.write("Flexible_Pivot_Data/ansysdata.csv",df,delim = '\t',decimal = ',',writeheader = true)
sol_ = str(inits_[:,:,5],beams_,nodes_,true)

using DataFrames,CSV
for (n,sol) in enumerate(sol.u[1:2])
    df = DataFrame(x = sol(0:0.01:1;idxs = 3).u .* beams_[n].l,y = sol(0:0.01:1;idxs = 4).u .* beams_[n].l,z = zeros(Float32,length(0:0.01:1)))
    CSV.write("Flexible_Pivot_Data/beam_$(n)_optimized.csv",df,delim = ' ',quotechar = '"',writeheader = false)
end



### ________________________________Testing stuff like gradients________________________________
using ForwardDiff
inits = randn(Float32,3,6,9)
function testgrads_str(str::Structure,init,beams,nodes)
    sol,beams,nodes_ = str(init,beams,nodes)
    sum(sol[:,2,:])
end 

fgr = ForwardDiff.gradient((x) -> testgrads_str(str,x,beams,nodes),inits[:,:,5])
l,bgr = withgradient((x) -> testgrads_str(str,x,beams,nodes),inits[:,:,5])

isapprox(fgr, bgr[1])
fgr
bgr[1]

sol,beams,nodes_ = str(inits[:,:,5],beams,nodes,false)

function test_grads_res(sol::AbstractArray{T,N},str,beams,nodes) where{N,T}
    r = zeros(T,3,6)
    t = BeamStructures.residuals!(r,str,sol,beams,nodes)
    sum(t)
end

fgr2 = ForwardDiff.gradient((x) -> test_grads_res(x,str,beams,nodes),sol)
l,bgr2 = withgradient((x) -> test_grads_res(x,str,beams,nodes),sol)

isapprox(bgr2[1],fgr2)
fgr2
bgr2[1]


function testgradsall_str(str::Structure,init::AbstractArray{T,N},beams,nodes) where{N,T}
    sol,beams,nodes_ = str(init,beams,nodes)
    # sol[:,1,:] .= 0
    r = zeros(T,3,6)
    t = BeamStructures.residuals!(r,str,sol,beams,nodes_)
    mean(abs2,t) 
end 

fgr = ForwardDiff.gradient((x) -> testgradsall_str(str,x,beams,nodes),inits[:,:,5])
l,bgr = withgradient((x) -> testgradsall_str(str,x,beams,nodes),inits[:,:,5])

isapprox(fgr ,bgr[1])

fgr[:,1:2]
bgr[1][:,1:2]
fgr[:,3:end] 
bgr[1][:,3:end]
fgr[:,3:end] .-bgr[1][:,3:end]

using ComponentArrays
function forwardgrad(str,ca::ComponentArray{T},types) where{T}
    beamtypes,nodetypes = types
    beams_ = Vector{BeamStructures.BeamElement}()
    for (i,bt) in enumerate(beamtypes)
        s = Symbol("beam$i")
        push!(beams_,bt{T}(ca.beams[s].l,ca.beams[s].h,ca.beams[s].w,ca.beams[s].κ0;θs = ca.beams[s].θs,E = ca.beams[s].E,θe = ca.beams[s].θe))
    end

    nodes_ = Vector{BeamStructures.Boundary}()
    for (i,nt) in enumerate(nodetypes)
        s = Symbol("node$i")
        push!(nodes_,nt{T}(ca.nodes[s].x,ca.nodes[s].y,ca.nodes[s].ϕ,ca.nodes[s].fx,ca.nodes[s].fy,ca.nodes[s].mz))
    end
    inits = ca.inits
    beams,nodes = BeamStructures.prepare(nodes_...,beams_...)
    testgradsall_str(str,inits,beams,nodes)
end 


gentup(struct_T) = NamedTuple{( fieldnames(struct_T)...,), Tuple{(fieldtype(struct_T,i) for i=1:fieldcount(struct_T))...}}

@generated function to_named_tuple_generated(x)
    nt = Expr(:quote, gentup(x))
    tup = Expr(:tuple)
    for i=1:fieldcount(x)
        push!(tup.args, :(getfield(x, $i)) )
    end
    return :($nt($tup))
end

to_named_tuple_generated(beam1)

# inits = randn(Float32,3,6)
# inits[:,3:4] .= 1
sol,beams,nodes_ = str(inits[5],beams,nodes)
beams
bca = ComponentArray(beam1 = to_named_tuple_generated(beams[1]),
                     beam2 = to_named_tuple_generated(beams[2]),
                     beam3 = to_named_tuple_generated(beams[3]),
                     beam4 = to_named_tuple_generated(beams[4]))

nca = ComponentArray(node1 = to_named_tuple_generated(nodes[1]),
                     node2 = to_named_tuple_generated(nodes[2]),
                     node3 = to_named_tuple_generated(nodes[3]),
                     node4 = to_named_tuple_generated(nodes[4]),
                     node5 = to_named_tuple_generated(nodes[5]))

# inits = sol
ca = ComponentArray(beams = bca, nodes = nca,inits = inits[5])

types = (beamtypes = [BeamStructures.gettype(t) for t in beams], nodetypes = [BeamStructures.type(t) for t in nodes])

forwardgrad(str,ca,types)

testgradsall_str(str,inits[5],beams,nodes)

gradsforward = ForwardDiff.gradient((x) -> forwardgrad(str,x,types),ca)

@run l,gradsbackward = withgradient((x,y,z)-> testgradsall_str(str,x,y,z),inits[5],beams,nodes)
l,gradsbackward = withgradient((x,y,z)-> testgradsall_str(str,x,y,z),inits[5],beams,nodes)

isapprox(gradsbackward[1],collect(gradsforward.inits))
gradsbackward[1][:,3:end]
collect(gradsforward.inits)[:,3:end]
gradsbackward[1][:,1:2]
collect(gradsforward.inits)[:,1:2]

gradsbackward[1].-collect(gradsforward.inits)

gradbeams = gradsbackward[2]

gradsforward.beams.beam1
gradbeams[1]

gradsforward.beams.beam2
gradbeams[2]

gradsforward.beams.beam3
gradbeams[3]

gradsforward.beams.beam4
gradbeams[4]

gradnodes = gradsbackward[3];

gradnodes[1]
gradsforward.nodes.node1

gradnodes[2]
gradsforward.nodes.node2

gradnodes[3]
gradsforward.nodes.node3

gradnodes[4]
gradsforward.nodes.node4

gradnodes[5]
gradsforward.nodes.node5
