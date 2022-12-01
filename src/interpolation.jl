
function data_interpolation(
    w_d, w;
    S = ones(length(vec(w))), 
    λ = [0.0], 
    m = nothing, 
    n = nothing, 
    L = nothing, 
    ℓ = nothing
)
    T = size(w,2)
    B = ss2BT_hankel(w_d, T)
    m = fit(LassoModel, B, vec(w), λ=λ, wts=vec(S))
    return reshape(predict(m), size(w))
end