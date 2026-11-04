using Test
using LinearAlgebra
using Random
using ForwardDiff
using Zygote
import ChainRulesCore as CRC
import BeamStructures as BS

include("effective_properties.jl")
include("utils.jl")

include("beam_derivatives.jl")
include("ground_structure_adjacency_derivatives.jl")
