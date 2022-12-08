
"""
    data_interpolation(w_d, w[, λ]; kw...)

Solve the data interpolation problem for `w` based on the exact trajectory `w_d`.
"""
function data_interpolation(
    w_d, w;
    S = ones(size(w)),
    m = nothing,
    n = nothing,
    L = nothing,
    ℓ = nothing,
    kernel=false
)
    T = size(w,2)
    B = ss2BT_hankel(w_d, T) # returns an orthonormal basis by default
    D = Diagonal(vec(S))
    if !isnothing(m) && !isnothing(n)
        B = low_rank_approximation(B, m*T+n)
    end
    w_h = B * ((D*B) \ (D * vec(w)))
    if kernel
        N = nullspace(D * B)
        return reshape(w_h, size(w)), N
    else
        return reshape(w_h, size(w))
    end
end

function data_interpolation!(w_h, A, w)
    ldiv!(vec(w_h), A, vec(w))
end

function data_interpolation(
    w_d, w, λ;
    S = ones(length(vec(w))), 
    m = nothing, 
    n = nothing, 
    L = nothing, 
    ℓ = nothing,
    kernel = false
)
    T = size(w,2)
    B = ss2BT_hankel(w_d, T)
    if !isnothing(m) && !isnothing(n)
        B = low_rank_approximation(B, m*T+n)
    end
    model = fit(LassoModel, B, vec(w), λ=λ, wts=vec(S), cd_tol=1e-10)
    w_h = predict(model)
    if kernel
        N = nullspace(B)
        return reshape(w_h, size(w)), N
    else
        return reshape(w_h, size(w))
    end
end

function data_interpolation!(w_h, A, w, λ)
    
end

function weights(w)
    f = x -> (ismissing(x) || isnan(x)) ? 0.0 : 1.0
    return [f(x) for x in w]
end

function obfuscate_output(w, p)
    q = size(w,1)
    w_ = Array{Union{Missing,eltype(w)},2}(undef, size(w))
    w_ .= w
    w_[q-p+1:q,:] .= missing
    return w_
end

"""
    data_simulation(w_d, w_ini, w_input)


"""
function data_simulation(w_d, w_ini, w_input)
    w = [w_ini w_input]
    return data_interpolation(w_d, coalesce.(w, zero(eltype(w))); S = weights(w))
end

function data_simulation(w_d, w_ini, w_input, p)
    return data_simulation(w_d, w_ini, obfuscate_output(w_input, p))
end

function impulse_response(sys::StateSpace, T)
    l = lag_datadriven(sys)
    n,m,_ = sizes(sys)
    imp = impulse(sys, l+1)
    w_d = random_trajectory(sys, (m+1)*T + n)
    return impulse_response(w_d, imp.y, m, T)
end

function impulse_response(w_d, y_ini, m, T)
    q = size(w_d, 1)
    y = zeros(eltype(w_d), q-m, T, m) 
    for i in 1:m
        w = impulse_response(w_d, y_ini[:,:,i], i, m, T)
        y[:,:,i] .= w[m+1:end,:]
    end
    return y
end

function impulse_response(w_d, y_ini, i, m, T)
    q = size(w_d,1)
    l = trajectory_length(y_ini)
    w_ini = zeros(q,l)
    w_ini[i,1] = 1.0
    w_ini[m+1:end,:] .= y_ini
    w_input = zeros(q,T-l)
    return data_simulation(w_d, w_ini, w_input, q-m)
end

function step_response(sys::StateSpace, T)
    l = lag_datadriven(sys)
    n,m,_ = sizes(sys)
    imp = step(sys, l+1)
    w_d = random_trajectory(sys, (m+1)*T + n)
    return step_response(w_d, imp.y, m, T)
end

function step_response(w_d, y_ini, m, T)
    q = size(w_d, 1)
    y = zeros(eltype(w_d), q-m, T, m) 
    for i in 1:m
        w = step_response(w_d, y_ini[:,:,i], i, m, T)
        y[:,:,i] .= w[m+1:end,:]
    end
    return y
end

function step_response(w_d, y_ini, i, m, T)
    q = size(w_d,1)
    l = trajectory_length(y_ini)
    w_ini = zeros(q,l)
    w_ini[i,:] .= 1.0
    w_ini[m+1:end,:] .= y_ini
    w_input = zeros(q,T-l)
    w_input[i,:] .= 1.0
    return data_simulation(w_d, w_ini, w_input, q-m)
end