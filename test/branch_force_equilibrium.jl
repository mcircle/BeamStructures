using Test, LinearAlgebra, Zygote
import ChainRulesCore as CRC
import BeamStructures as BS

@testset "Branch force endpoints" begin
    beams=(b1=BS.Beam(35.0,1.0,5.0,0.01;θs=0.2),
           b2=BS.Beam(45.0,1.0,5.0,-0.01;θs=-0.1),
           b3=BS.Beam(30.0,1.0,5.0,0.005;θs=0.3))
    nodes=(n1=BS.Clamp(0.0,0.0,0.0,0.0,0.0,0.0),
           n2=BS.Branch(25.0,20.0,0.0,0.2,-0.1,0.05),
           n3=BS.Clamp(50.0,5.0,0.0,0.0,0.0,0.0))
    states=reshape(collect(range(-0.3,0.4;length=42)),7,2,3)
    adj=[0.0 0.25 0.55;0.25 0.0 0.8;0.55 0.8 0.0]
    incoming=BS.scaleforce(beams[1],states[[1,5,6],2,1])
    outgoing=BS.scaleforce(beams[3],states[[1,5,6],1,3])
    external=nodes.n2[[6,4,5]]
    for weighted in (false,true)
        a,b=weighted ? (0.25,0.8) : (1.0,1.0)
        forward(y)=BS.residuals!(zeros(3,4),adj,y,beams,nodes;dfac=weighted)
        residual=forward(states)
        @test residual[:,1] ≈ a.*incoming .- b.*outgoing .- external
        seed=zeros(3,4);seed[2,1]=1.0
        _,back=CRC.rrule(BS.residuals!,zeros(3,4),adj,states,beams,nodes,Val(weighted))
        grads=back(seed)
        @test grads[6].n2.fx ≈ -1.0
        for (component,side,edge) in ((5,2,1),(5,1,3),(5,1,1),(5,2,3))
            direction=zeros(size(states));direction[component,side,edge]=1
            f(t)=forward(states .+ t.*direction)[2,1]
            h=1e-5; fd=(f(h)-f(-h))/(2h)
            @test isapprox(grads[4][component,side,edge],fd;rtol=1e-6,atol=1e-8)
        end
    end
end
