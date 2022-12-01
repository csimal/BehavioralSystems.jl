import Base.getindex

struct HankelMatrix{M,T} <: AbstractMatrix{T} 
    w::M
    L::Int
end

function Base.size(H::HankelMatrix)
    q,T = size(H.w)
    return (q*(H.L), T-(H.L)+1)
end

function Base.getindex(H::HankelMatrix, i, j)
    q = size(H.w, 1)
    if rem(i, q) == 0

    else

    end
    t = div(i, q) + 1
    α = rem(i, q)
    β = j+t-1
    return H.w[α,β]
end