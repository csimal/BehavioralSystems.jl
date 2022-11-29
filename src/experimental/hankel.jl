
struct HankelMatrix{M<:AbstractVecOrMat{T},T} <: AbstractMatrix{T} where {T}
    w::M
    L::Int
end