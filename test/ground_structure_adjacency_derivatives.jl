function adjacency_residual_fixture(::Type{T}=Float64) where{T}
    beams = (
        Beam_1=BS.Beam(T(35),T(1),T(5),T(0.01);θs=T(0.2)),
        Beam_2=BS.Beam(T(45),T(1),T(5),T(-0.01);θs=T(-0.1)),
        Beam_3=BS.Beam(T(30),T(1),T(5),T(0.005);θs=T(0.3)),
    )
    nodes = (
        Node_1=BS.Clamp(T(0),T(0),T(0),T(0),T(0),T(0)),
        Node_2=BS.Branch(T(25),T(25),T(0),T(0.2),T(-0.1),T(0.05)),
        Node_3=BS.Branch(T(50),T(10),T(0.1),T(-0.1),T(0.15),T(-0.02)),
    )
    states = reshape(T.(range(-0.3,0.4;length=42)),7,2,3)
    beams,nodes,states
end

function adjacency_matrix(factors)
    z = zero(eltype(factors))
    [z factors[1] factors[2];
     factors[1] z factors[3];
     factors[2] factors[3] z]
end

function adjacency_residual_loss(factors,beams,nodes,states,weights;dfac)
    residual = zeros(eltype(factors),3,5)
    adjacency = adjacency_matrix(factors)
    sum(BS.residuals!(residual,adjacency,states,beams,nodes;dfac=dfac) .* weights)
end

@testset "GroundStructure adjacency residual switch" begin
    factors = [0.25,0.55,0.8]
    beams,nodes,states = adjacency_residual_fixture()
    weights = reshape(collect(1.0:15.0),3,5)
    adjacency = adjacency_matrix(factors)
    default_residual = zeros(3,5)
    explicit_residual = zeros(3,5)
    weighted_residual = zeros(3,5)

    BS.residuals!(default_residual,adjacency,states,beams,nodes)
    BS.residuals!(explicit_residual,adjacency,states,beams,nodes;dfac=false)
    BS.residuals!(weighted_residual,adjacency,states,beams,nodes;dfac=true)

    @test default_residual == explicit_residual
    @test weighted_residual != default_residual

    reverse = only(Zygote.gradient(
        x -> adjacency_residual_loss(x,beams,nodes,states,weights;dfac=true),
        factors))
    forward = ForwardDiff.gradient(
        x -> adjacency_residual_loss(x,beams,nodes,states,weights;dfac=true),
        factors)
    @test reverse ≈ forward rtol=1e-8 atol=1e-10
    @test any(!iszero,reverse)

    disabled = only(Zygote.gradient(
        x -> adjacency_residual_loss(x,beams,nodes,states,weights;dfac=false),
        factors))
    @test disabled ≈ zeros(3) atol=1e-12
end
