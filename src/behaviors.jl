
abstract type Behaviour end

struct StateSpaceBehaviour <: Behaviour end

struct HankelBehaviour <: Behaviour end

struct KernelBehaviour{T} <: Behaviour
    R::T
end

"""
ss2BT_modelbased(sys,T)

Compute a basis of the behaviour ℬ_T using the analytic solution.
"""
function ss2BT_modelbased(sys,Π,T)
A,B,C,D = sys.A,sys.B,sys.C,sys.D
n = size(A,1)
m = size(B,2)
p = size(C,1)
q = m + p
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
function ss2BT_datadriven(sys,T)
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

function ss2BT_hankel(sys, L)
n,m,_ = sizes(sys)
w = random_trajectory(sys,(m+1)*L+n)
return ss2BT_hankel(w, L)
end

function ss2r_modelbased(sys)
    n,m,_ = sizes(sys)
    return kernel_basis(ss2BT_modelbased(sys, nothing, (m+1)*L+n))
end

function ss2r_datadriven(w, L)
    kernel_basis(hankel_matrix(w,L))
end

"""
    lag_model_based(sys)

Compute the lag of a dynamical system using the analytical solution.
"""
function lag_model_based(sys)
n,_,p = sizes(sys)
A, C = sys.A, sys.C
O = zeros((n-1)*p, n)
l = 1
M = C
rows = axes(C,1)
O[rows .+ (l-1)*p,:] .= M
while rank(view(O,1:l*p,1:n)) < n && l < n
    l += 1
    M  = M*A
    O[rows .+ (l-1)*p,:] .= M
end
return l
end

"""
    lag_datadriven(sys)

Compute the lag of a dynamical system by studing the dimension of its restricted behaviors.
"""
function lag_datadriven(sys)
n,m,_ = sizes(sys)
w = random_trajectory(sys, (m+1)*n +n)
l = n
while rank(ss2BT_hankel(w,l)) == m*l + n
    l -= 1
end
return l+1
end

function lag_datadriven(w::AbstractMatrix)
    
end

function r2BT(sys)
    
end