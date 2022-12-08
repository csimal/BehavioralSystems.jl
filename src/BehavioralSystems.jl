module BehavioralSystems

using ControlSystemsBase
using Optimization
using LinearAlgebra
using Random
using Lasso

include("utils.jl")
export drss, random_trajectory, sizes, canonical_permutation

include("matrices.jl")
export hankel_matrix, hankel_projection, antidiagonal
export multiplication_matrix, multiplication_projection
export range_basis, kernel_basis, compare_spans

include("behaviors.jl")
export ss2BT_datadriven, ss2BT_hankel, ss2BT_modelbased
export ss2r_modelbased, ss2r_datadriven

include("complexity.jl")
export lag_modelbased, lag_datadriven
export complexity_modelbased, complexity_datadriven

include("mpum.jl")
export complexity_mpum, most_powerful_unfalsified_model

include("interpolation.jl")
export data_interpolation, data_simulation
export impulse_response, step_response

end
