using BehavioralSystems
using ControlSystemsBase
using LinearAlgebra
using Test

@testset "BehavioralSystems.jl" begin
    include("matrices.jl")
    include("utils.jl")
    include("behaviors.jl")
end
