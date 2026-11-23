using Test, BeamStructures, LinearAlgebra, Zygote
const BS=BeamStructures

@testset "Constrained tangent system" begin
    for T in (Float32,Float64), E in (100,200)
        beam=BS.Beam(T(10),T(2),T(1),T(0);E=T(E),θs=T(0),θe=T(0))
        beams=(b1=beam,)
        nodes=(n1=BS.Clamp(T(0),T(0),T(0),T(0),T(0),T(0)),
               n2=BS.Clamp(T(10),T(0),T(0),T(0),T(0),T(0)))
        adj=T[0 1;1 0];str=BS.GroundStructure()
        states=first(str(zeros(T,3,1),beams,nodes,adj))
        system=BS.admittance_matrix(states,adj,str,values(beams))
        @test system isa BS.TangentSystem
        k=BS.effective_stiffness(system,[1,2,3,4,5]=>[6])
        c=BS.effective_compliance(system,[1,2,3,4,5]=>[6])
        @test isapprox(k,T(0.008E);rtol=T(3e-3))
        @test isapprox(k,inv(c+eps(T));rtol=T(3e-3))
        @test_throws ArgumentError BS.effective_stiffness(system,[1]=>[1])
        # The current inextensible horizontal beam cannot change its axial length.
        @test_throws BS.TangentCompatibilityError BS.effective_stiffness(system,[1,2,3,4,6]=>[5])
        function objective(a)
            w=T[0 1;a 0]
            BS.effective_stiffness(BS.admittance_matrix(states,w,str,values(beams)),[1,2,3,4,5]=>[6])
        end
        a=T(0.6);h=T(1e-3)
        reverse=only(Zygote.gradient(objective,a))
        finite=(objective(a+h)-objective(a-h))/(2h)
        @test isapprox(reverse,finite;rtol=T(3e-2),atol=T(1e-3))
        function compliance_objective(weight)
            w=T[0 1;weight 0]
            BS.effective_compliance(BS.admittance_matrix(states,w,str,values(beams)),[1,2,3,4,5]=>[6])
        end
        cad=only(Zygote.gradient(compliance_objective,a))
        cfd=(compliance_objective(a+h)-compliance_objective(a-h))/(2h)
        @test isapprox(cad,cfd;rtol=T(3e-2),atol=T(1e-3))
    end
end

@testset "Clamp-Branch-Clamp tangent and adjacency gradients" begin
    positions=((0.0,0.0),(5.0,2.0),(10.0,0.0))
    function beam_at(i,j,E)
        dx=positions[j][1]-positions[i][1];dy=positions[j][2]-positions[i][2]
        angle=atan(dy,dx)
        BS.Beam(hypot(dx,dy),1.0,1.0,0.0;E=E,θs=angle,θe=angle)
    end
    nodes=(n1=BS.Clamp(0.,0.,0.,0.,0.,0.),n2=BS.Branch(5.,2.,0.,0.,0.,0.),
           n3=BS.Clamp(10.,0.,0.,0.,0.,0.))
    base=[0.0 0.7 0.4;0.7 0.0 0.8;0.4 0.8 0.0]
    dofs=[1,2,3,7,8]=>[9]
    str=BS.GroundStructure()
    for E in (100.0,200.0)
        beams=(b1=beam_at(1,2,E),b2=beam_at(1,3,E),b3=beam_at(2,3,E))
        full=[0.0 1 1;1 0 1;1 1 0]
        states=first(str(zeros(3,4),beams,nodes,full))
        k=BS.effective_stiffness(BS.admittance_matrix(states,full,str,values(beams)),dofs)
        # Reference from independent residual/reaction finite differences, with corrected Branch force ends.
        @test isapprox(k,E/100 * 0.1928511;rtol=2e-3)
        function objective(a,edge)
            adj=[0.0 0.7 0.4;edge==1 ? a : 0.7 0.0 0.8;edge==2 ? a : 0.4 edge==3 ? a : 0.8 0.0]
            BS.effective_stiffness(BS.admittance_matrix(states,adj,str,values(beams)),dofs)
        end
        for (edge,value) in ((1,0.7),(2,0.4),(3,0.8))
            h=1e-4
            ad=only(Zygote.gradient(a->objective(a,edge),value))
            fd=(objective(value+h,edge)-objective(value-h,edge))/(2h)
            @test isapprox(ad,fd;rtol=3e-3,atol=1e-5)
        end
        @test isapprox(objective(0.0,2),objective(1e-5,2);atol=2e-5)
        @test_throws BS.TangentCompatibilityError BS.effective_stiffness(
            BS.admittance_matrix(states,zeros(3,3),str,values(beams)),dofs)
    end
end

@testset "Deformed open Clamp-Branch chain" begin
    positions=((0.0,0.0),(5.0,2.0),(10.0,0.0))
    function beam_at(i,j)
        dx=positions[j][1]-positions[i][1];dy=positions[j][2]-positions[i][2]
        angle=atan(dy,dx)
        BS.Beam(hypot(dx,dy),1.0,1.0,0.0;E=100.0,θs=angle,θe=angle)
    end
    beams=(b1=beam_at(1,2),b2=beam_at(1,3),b3=beam_at(2,3))
    adj=[0.0 1 0;1 0 1;0 1 0];str=BS.GroundStructure()
    nodes(u)=(n1=BS.Clamp(0.,0.,0.,0.,0.,0.),n2=BS.Branch(5.,2.,0.,0.,0.,0.),
              n3=BS.Clamp(10.,u,0.,0.,0.,0.))
    function evaluate(q,u)
        local_states,_,moved=str(reshape(q,3,4),beams,nodes(u),adj)
        local_residual=vec(BS.residuals!(zeros(3,4),adj,local_states,beams,moved;dfac=true))
        reaction=BS.scaleforce(beams[3],local_states[[1,5,6],2,3])[3]
        local_residual,reaction,local_states
    end
    function derivatives(q,u,h)
        J=zeros(12,12);fq=zeros(12)
        for j in 1:12
            d=zeros(12);d[j]=h
            rp,fp,_=evaluate(q+d,u);rm,fm,_=evaluate(q-d,u)
            J[:,j]=(rp-rm)/(2h);fq[j]=(fp-fm)/(2h)
        end
        rp,fp,_=evaluate(q,u+h);rm,fm,_=evaluate(q,u-h)
        J,(rp-rm)/(2h),fq,(fp-fm)/(2h)
    end
    u=0.01;q=zeros(12)
    for _ in 1:8
        residual,_,_=evaluate(q,u)
        norm(residual)<1e-7 && break
        J,_,_,_=derivatives(q,u,1e-3)
        q=q-pinv(J;rtol=1e-8)*residual
    end
    residual,_,states=evaluate(q,u)
    @test norm(residual)<1e-7
    J,ru,fq,fu=derivatives(q,u,1e-3)
    reference=fu-dot(fq,pinv(J;rtol=1e-8)*ru)
    actual=BS.effective_stiffness(BS.admittance_matrix(states,adj,str,values(beams)),[1,2,3,7,8]=>[9])
    @test isapprox(actual,reference;rtol=1e-3,atol=1e-6)
end

@testset "Float32 Clamp-Branch tangent gradient" begin
    T=Float32
    pos=((T(0),T(0)),(T(5),T(2)),(T(10),T(0)))
    function beam_at(i,j)
        dx=pos[j][1]-pos[i][1];dy=pos[j][2]-pos[i][2];angle=atan(dy,dx)
        BS.Beam(hypot(dx,dy),T(1),T(1),T(0);E=T(100),θs=angle,θe=angle)
    end
    beams=(b1=beam_at(1,2),b2=beam_at(1,3),b3=beam_at(2,3))
    nodes=(n1=BS.Clamp(T(0),T(0),T(0),T(0),T(0),T(0)),
           n2=BS.Branch(T(5),T(2),T(0),T(0),T(0),T(0)),
           n3=BS.Clamp(T(10),T(0),T(0),T(0),T(0),T(0)))
    str=BS.GroundStructure();base=T[0 0.7 0.4;0.7 0 0.8;0.4 0.8 0]
    states=first(str(zeros(T,3,4),beams,nodes,base));dofs=[1,2,3,7,8]=>[9]
    f(a)=BS.effective_stiffness(BS.admittance_matrix(states,T[0 0.7 0.4;a 0 0.8;0.4 0.8 0],str,values(beams)),dofs)
    a=T(0.7);h=T(1e-3)
    ad=only(Zygote.gradient(f,a));fd=(f(a+h)-f(a-h))/(2h)
    @test isfinite(f(a))
    @test isapprox(ad,fd;rtol=T(1e-2),atol=T(1e-4))
end
