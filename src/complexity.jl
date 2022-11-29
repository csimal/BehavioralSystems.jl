
"""
    lag_model_based(sys)

Compute the lag of a dynamical system using the analytical solution.

This is done by constructing observability matrices of increasing order until they reach full rank.
"""
function lag_modelbased(sys::StateSpace)
    n,_,p = sizes(sys)
    A, C = sys.A, sys.C
    O = zeros(n*p, n)
    l = 1
    M = C
    rows = axes(C,1)
    O[rows .+ (l-1)*p,:] .= M
    while rank(view(O,1:l*p,1:n)) < n && l < n
        l += 1
        M = M*A
        O[rows .+ (l-1)*p,:] .= M
    end
    return l
end

"""
    lag_datadriven(sys::StateSpace)
    lag_datadriven(w::AbstractMatrix)
    lag_datadriven(w, m, n)

Compute the lag of a dynamical system by studing the dimension of its restricted behaviors.

- If the input is a `Statespace`, the lag is computed by generating a random trajectory of sufficient length.
- If the input is a data matrix representing a trajectory, the lag is estimated by looking at the ranks of hankel matrices of the data. If the provided trajectory is too short, the result will be incorrect.
- If the number of states `n` and input `m` are known, these can provided alongside the data matrix.
"""
function lag_datadriven(sys::StateSpace)
    n,m,_ = sizes(sys)
    w = random_trajectory(sys, (m+1)*n +n)
    return lag_datadriven(w, m, n)
end

function lag_datadriven(w::AbstractMatrix, m, n)
    l_min = 1
    l_max = n
    while l_min + 1 < l_max
        l_mid = div(l_min + l_max, 2)
        if rank(ss2BT_hankel(w, l_mid)) == m*l_mid + n
            l_max = l_mid
        else
            l_min = l_mid
        end
    end
    return l_max
end

function lag_datadriven(w::AbstractMatrix)
    T = size(w,2)
    l_max = _findmax_rank(w) # after that point, we don't have enough data to correctly estimate the rank
    r_max = rank(ss2BT_hankel(w,l_max))
    m = r_max - rank(ss2BT_hankel(w,l_max-1)) # get the best estimate we can for m
    n = r_max - m*l_max
    return lag_datadriven(w, m, n)
end

"""
    _findmax_rank(w)

Find `L` such that the rank of `hankel_matrix(w,L)` is maximal.
"""
function _findmax_rank(w::AbstractMatrix)
    T = size(w,2)
    t_min = 1
    t_max = T
    while t_min + 1 < t_max
        t_mid = div(t_min+t_max, 2)
        if rank(ss2BT_hankel(w, t_mid)) < rank(ss2BT_hankel(w,t_mid+1))
            t_min = t_mid
        else
            t_max = t_mid
        end
    end
    return t_max
end

"""
    complexity_modelbased(sys::StateSpace)

Compute the complexity of a discrete-time LTI system.

The complexity is given by a triplet `(n,m,l)` where `n` is the number of states, `m` is the number of inputs and `l` is the lag.
"""
function complexity_modelbased(sys::StateSpace)
    n,m,_ = sizes(sys)
    l = lag_modelbased(sys)
    return m,n,l
end

"""
    complexity_datadriven(w::AbstractMatrix)

Estimate the complexity of a LTI system based on one of its trajectories.


"""
function complexity_datadriven(w::AbstractMatrix)
    l_max = _findmax_rank(w) # after that point, we don't have enough data to correctly estimate the rank
    r_max = rank(ss2BT_hankel(w,l_max))
    m = r_max - rank(ss2BT_hankel(w,l_max-1)) # get the best estimate we can for m
    n = r_max - m*l_max
    return m, n, lag_datadriven(w, m, n)
end