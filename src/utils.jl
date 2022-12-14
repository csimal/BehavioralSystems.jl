
# Shamelessy ported from python-control
"""
    drss(states, inputs, outputs, [rng])

Generate a random discrete time state space LTI system

## Example
```julia
sys = drss(3,4,2)
u = rand(4, 10) # Time must be on second dimension
x₀ = rand(3)
sol = lsim(sys, u, 1:10, x₀)
```
"""
function drss(states, inputs, outputs, rng=Random.GLOBAL_RNG)
    poles = zeros(ComplexF64, states)
    p_repeat = 0.05 # Probability of repeating a pole
    p_real = 0.6 # Probability of generating a real pole
    p_bc_mask = 0.8 # Probability of masking an element of B,C
    p_d_mask = 0.3 # Probability of masking an element of D
    p_d_zero = 0.5

    i = 1
    while i < states
        if rand(rng) < p_repeat && i > 1 && i < states
            # Copy pole
            if imag(poles[i-1]) == 0.0
                poles[i] = poles[i-1]
                i += 1
            else
                poles[i:i+1] .= poles[i-2:i-1]
                i += 2
            end
        elseif rand(rng) < p_real || i == states
            # Add real pole
            poles[i] = 2.0 * rand(rng) - 1.0
            i += 1
        else # Add Complex pole
            r = rand(rng)
            ϕ = 2.0 * π * rand(rng)
            poles[i] = r * (cos(ϕ) + sin(ϕ)*im)
            poles[i+1] = conj(poles[i])
            i += 2
        end
    end
    A = zeros(states, states)
    i = 1
    while i <= states
        if imag(poles[i]) == 0
            A[i,i] = real(poles[i])
            i += 1
        else
            a, b = real(poles[i]), imag(poles[i])
            A[i:i+1,i:i+1] .= [a b; -b a]
            i += 2
        end
    end
    T = rand(rng, states, states)
    A = (T \ A) * T
    B = randn(rng, states, inputs)
    C = randn(rng, outputs, states)
    D = randn(rng, outputs, inputs)

    B_mask = rand(rng, states, inputs) .< p_bc_mask
    C_mask = rand(rng, outputs, states) .< p_bc_mask
    D_mask = rand(rng, outputs, inputs) .< p_d_mask
    B .*= B_mask 
    C .*= C_mask 
    D .*= D_mask 
    ss(A, B, C, D, 1.0)
end

"""
    random_trajectory(sys, T)

Compute a random trajectory of the discrete time LTI system `sys` up to time `T`.

The output is a matrix `w = [u; y]` containing the input and the output in its rows and succesive time steps in its columns.
"""
function random_trajectory(A,B,C,D,T)
    random_trajectory(ss(A,B,C,D,1.0), T)
end

function random_trajectory(sys, T)
    n,m,_ = sizes(sys)
    u = rand(m,T)
    x₀ = rand(n)
    y = lsim(sys, u, 0:T-1; x0=x₀)[1]
    return [u; y]
end

function random_trajectory!(w, sys, T)
    n,m,p = sizes(sys)
    if size(w) != (m+p,T)
        error("w needs to be of size ($(m+p),$T)")
    end
    u = rand(m,T)
    w[1:m,:] .= u
    x₀ = rand(n)
    y = lsim(sys, u, 0:T-1; x0=x₀)[1]
    w[m+1:end,:] .= y
    return w
end

"""
    sizes(sys)

Compute the number of states, inputs and output of an LTI system.

## Examples
```jldoctest
julia> sys = drss(3,2,4)
ControlSystemsBase.StateSpace{ControlSystemsBase.Discrete{Float64}, Float64}
...
julia> sizes(sys)
(3, 2, 4)
```
"""
function sizes(sys)
    n = size(sys.A,1)
    m = size(sys.B,2)
    p = size(sys.C,1)
    return (n,m,p)
end

"""
    canonical_permutation(m,p,T)

Return the canonical permutation Π that maps a vector `w = [u; y]` with `u` of length `m*T` and `y` of size `p*T` to a vector `[u_1; y_1; ...; u_T; y_T]`.

`u_t` and `y_t` are the successive slices of size `m` and `p` of `u` and `y` respectively.

The result is a vector of size `(m+p)*T` representing the permutation for use with [`permute!`](@ref) or [`permutedims`](@ref)
"""
function canonical_permutation(m, p, T)
    q = m+p
    Π = Vector{Int}(undef, (m+p)*T)
    for t in 1:T
        Π[((t-1)*q) .+ (1:m)] .= ((t-1)*m) .+ (1:m)
        Π[((t-1)*q) .+ ((m+1):q)] .= (T*m+(t-1)*p)
        Π[((t-1)*q) .+ ((m+1):q)] .+= (1:p)
    end
    return Π
end

trajectory_length(w::AbstractVector) = length(w)
trajectory_length(w::AbstractMatrix) = size(w,2)