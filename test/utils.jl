@testset "Utility naming and preparation" begin
    nodes = (
        BS.Clamp(0.0,0.0,0.0,0.0,0.0,0.0),
        BS.Branch(1.0,0.0,0.0,0.0,0.0,0.0),
    )
    beams = (
        BS.Beam{Float64}(1,1,1,0),
        BS.Beam{Float32}(1,1,1,0),
    )

    named_nodes = BS.getnames(nodes...)
    named_beams = BS.getnames(beams...)
    @test keys(named_nodes) == (:Node_1,:Node_2)
    @test keys(named_beams) == (:Beam_1,:Beam_2)
    @test values(named_nodes) == nodes
    @test values(named_beams) == beams

    prepared_beams,prepared_nodes = BS.prepare(
        beams[1],nodes[1],beams[2],nodes[2])
    @test prepared_beams == named_beams
    @test prepared_nodes == named_nodes
end
