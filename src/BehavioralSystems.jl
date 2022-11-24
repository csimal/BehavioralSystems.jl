module BehavioralSystems

using ControlSystemsBase
using Optimization
using LinearAlgebra
using Random

include("utils.jl")
export drss, random_trajectory, sizes

include("matrices.jl")
export hankel_matrix, hankel_projection, antidiagonal
export multiplication_matrix, multiplication_projection

include("behaviors.jl")

end
