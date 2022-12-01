
function most_powerful_unfalsified_model(w::AbstractMatrix)
    _, _, l = complexity_mpum(w)
    return ss2r_datadriven(w,l+1)
end

"""
    complexity_mpum(w)

Compute the complexity of the most powerful unfalsified model that generates `w`.
"""
function complexity_mpum(w::AbstractMatrix)
    q,T = size(w)
    l_max = floor((T+1)/(q+1))
    l_mpum = 1
    m_mpum, n_mpum = estimate_complexity(w, l_mpum)
    for l in 2:l_max
        m, n = estimate_complexity(w,l)
        if m < m_mpum || n < n_mpum
            l_mpum = l
            m_mpum = m
            n_mpum = n
        end
    end
    return m_mpum, n_mpum, l_mpum
end

"""
    estimate_complexity(w, l)

Estimate the number of inputs `m` and states `n` of a LTI system based on one of its trajectories and assuming lag `l`.
"""
function estimate_complexity(w::AbstractMatrix, l)
    r_1 = rank(hankel_matrix(w,l+1))
    r_2 = rank(hankel_matrix(w,l+2))
    m = r_2 - r_1
    n = r_1 - (l+1)*m
    return m, n
end