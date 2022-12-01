
"""
    BehaviorRepresentation

The abstract for representations of the behavior of an LTI system.
"""
abstract type BehaviorRepresentation end

"""
    StateSpaceRepresentation{T,P} <: BehaviorRepresentation

A representation of the behavior using state space parameters A, B, C, D and a permutation that distinguishes inputs from outputs.

`w = [u, y]` is trajectory of the system if there exists `x` such that
- ``σx = Ax + Bu``
- ``y = Cx + Du``
"""
struct StateSpaceRepresentation{T<:StateSpace,P} <: BehaviorRepresentation 
    sys::T
    Π::P
end

"""
    DataDrivenRepresentation{T} <: BehaviorRepresentation

A representation of the behavior using the Hankel matrix ℋ of a trajectory.

`w` is a trajectory of the system if ``w \\in \\range (ℋ)``
"""
struct DataDrivenRepresentation{T<:AbstractMatrix} <: BehaviorRepresentation
    H::T
end

basis(ddr::DataDrivenRepresentation) = range_basis(ddr.H)

"""
    KernelRepresentation{T} <: BehaviorRepresentation

A representation of the behavior as the kernel of an operator `R`.

`w` is a trajectory of the system if ``R(σ)w = 0``
"""
struct KernelRepresentation{T<:AbstractMatrix} <: BehaviorRepresentation
    R::T
end

"""
ss2BT_modelbased(sys,T)

Compute a basis of the behaviour ℬ_T using the analytic solution.
"""
function ss2BT_modelbased(sys::StateSpace,Π,T)
    A,B,C,D = sys.A,sys.B,sys.C,sys.D
    n = size(A,1)
    m = size(B,2)
    p = size(C,1)
    M = [zeros(m*T,n) I; zeros(p*T, n+m*T)]
    O = copy(C) 
    M[(m*T+1):(m*T+p), 1:n] .= C
    for t in 1:T
        rows = m*T + (t-1)*p .+ (1:p)
        cols = n + (t-1)*m .+ (1:m)
        M[rows,cols] .= D
    end
    for t in 2:T
        rows = m*T + (t-1)*p .+ (1:p)
        M[rows,1:n] .= O
        for s in t:T
            rows_ = m*T + (s-1)*p .+ (1:p)
            cols_ = n + (s-t)*m .+ (1:m)
            M[rows_,cols_] .= O*B
        end
        O *= A
    end
    return M
end

"""
ss2BT_datadriven(sys, T)

Compute an orthonormal basis of ℬ_T by sampling q*T trajectories.
"""
function ss2BT_datadriven(sys::StateSpace,T)
    _,m,p = sizes(sys)
    q = m + p
    W = zeros(q*T, q*T)
    for i in 1:q*T
        w = random_trajectory(sys,T)
        W[:,i] .= vec(w)
    end
    return range_basis(W)
end


"""
    ss2BT_hankel(w, L)

Compute an orthonormal basis of ℬ_L using the Hankel matrix of a trajectory `w`.
"""
function ss2BT_hankel(w,L)
    ℋ = hankel_matrix(w,L)
    return range_basis(ℋ)
end

function ss2BT_hankel(sys::StateSpace, L)
    n,m,_ = sizes(sys)
    w = random_trajectory(sys,(m+1)*L+n)
    return ss2BT_hankel(w, L)
end

function ss2r_modelbased(sys::StateSpace)
    n,m,p = sizes(sys)
    q = m + p
    l = lag_modelbased(sys)
    M = ss2BT_modelbased(sys, nothing, l+1)
    O = M[m*(l+1)+1:end-p, 1:n]
    C = M[m*(l+1)+1:end,n+1:end]
    R_y = sys.C * sys.A^(l) * pinv(O)
    R_u = R_y * C[1:end-p,:] + C[end-p+1:end,:]

    return [R_u -R_y -I]
end

"""
    ss2r_datadriven(w,L)

Compute a Kernel Representation of the behavior of a state space system based on one of its trajectories.
"""
function ss2r_datadriven(w, L)
    kernel_basis(hankel_matrix(w,L))'
end

function ss2r_datadriven(sys::StateSpace)
    n,m,_ = sizes(sys)
    l = lag_datadriven(sys) + 1
    w = random_trajectory(sys, (m+1)*l+n)
    return ss2r_datadriven(w, l)
end

"""
    r2BT(R,T)

Compute a basis for the restricted behavior ℬ_T from the parameter `R` of a kernel representation of ℬ
"""
function r2BT(R,T)
    return nullspace(multiplication_matrix(R,T))
end

function r2BT(ℬ::KernelRepresentation, T)
    return r2BT(ℬ.R, T)
end