
abstract type BehaviorRepresentation end

struct StateSpaceRepresentation{T<:StateSpace,P} <: BehaviorRepresentation 
    sys::T
    Π::P
end

struct DataDrivenRepresentation{T<:AbstractMatrix} <: BehaviorRepresentation
    H::T
end

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
function ss2BT_hankel(w::AbstractMatrix,L)
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
    l = lag_model_based(sys)
    M = ss2BT_modelbased(sys, nothing, l)
    rows = m*l+1:(m+p)*l
    O = M[rows, 1:n]
    C = M[rows,n+1:end]
    tmp = O*pinv(O)
    
    return [tmp*C -(tmp+I)]
end

"""
    ss2r_datadriven(w,L)

Compute a Kernel Representation of 
"""
function ss2r_datadriven(w, L)
    kernel_basis(hankel_matrix(w,L))
end

function ss2r_datadriven(sys::StateSpace)

end

"""
    r2BT(R,T)
"""
function r2BT(R,T)
    return kernel_basis(multiplication_matrix(R,T))
end