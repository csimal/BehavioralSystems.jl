

"""
    hankel_matrix(W, L)

Compute the block Hankel matrix of `W` by taking slices of length `L`.

The resulting matrix depends on the type of `W`.
* If `W` is a vector of scalars, the result has `L` lines with each column being successive slices of length `L` of `W`.
* If `W` is a matrix, the columns of `W` are taken as blocks, and the result has `L` block lines, constructed as in the scalar case.
* If `W` is a vector of vectors or matrices, the Hankel matries of the elements of `W` are concatenated along their columns.

## Examples
```jldoctest
julia> hankel_matrix(1:5, 2)
2×4 Matrix{Int64}:
 1  2  3  4
 2  3  4  5

 julia> x = [1:5;;1:5]'
 2×5 adjoint(::Matrix{Int64}) with eltype Int64:
  1  2  3  4  5
  1  2  3  4  5

julia> hankel_matrix(x, 2)
4×4 Matrix{Int64}:
 1  2  3  4
 1  2  3  4
 2  3  4  5
 2  3  4  5

julia> hankel_matrix([1:5,1:3], 2)
2×6 Matrix{Int64}:
 1  2  3  4  1  2
 2  3  4  5  2  3
```
"""
function hankel_matrix(W::AbstractMatrix, L)
    if L < 1
        throw(DomainError("L should be at least 1. Got $L"))
    end
    q, T = size(W)
    if L > T
        throw(DomainError("Not enough samples. T should be larger than L. Got T=$T and L=$L"))
    end
    ℋ = zeros(eltype(W), L*q, T-L+1)
    for j in 1:T-L+1, i in 1:q, t in 1:L
        ℋ[(t-1)*q+i,j] = (W[i,(j-1)+t])
    end
    return ℋ
end

function hankel_matrix(W::AbstractVector{T}, L) where {T<:Union{AbstractVector,AbstractMatrix}}
    return hcat(hankel_matrix.(W,L)...)
end

function hankel_matrix(W::AbstractVector, L)
    return hankel_matrix(W', L)
end

"""
    antidiagonal(t,m,n)

Return an iterator over the indices of the `t`-th antidiagonal of a matrix with `m` rows and `n` columns.

The `t`-th antidiagonal is the set of entries `(i,j)` that satisfy `i+j = t+1`.
"""
function antidiagonal(t,m,n)
    if m < 0 || n < 0
        throw(DomainError("m and n should be non-negative. Got m=$m and n=$n"))
    end
    return Iterators.filter(Iterators.map(i -> (i,t+1-i), 1:m)) do (i,j)
        checkindex(Bool, 1:m, i) && checkindex(Bool, 1:n, j)
    end
end

"""
    hankel_projection(D,L)

Compute the orthogonal projection of `D` on the space of `L`-Hankel matrices.
"""
function hankel_projection(D,L)
    if L < 1
        throw(DomainError("L should be at least 1. Got $L"))
    end
    q = size(D,1) ÷ L
    T = size(D,2) + L - 1
    w = zeros(q,T)
    for t in axes(w,2)
        count = 0
        for (i,j) in antidiagonal(t, L, size(D,2))
            w[:,t] .+= D[(i-1)*q .+ (1:q), j]
            count += 1
        end
        w[:,t] ./= count
    end
    return hankel_matrix(w,L)
end

function hankel_projection(D,L,n)
    idx = cumsum(n)
    pushfirst!(idx, 0)
    return hcat([hankel_projection(
        view(D,axes(D,1),idx[i]+1:idx[i+1]),
        L
    ) for i in 1:length(idx)-1]...)
end

function multiplication_matrix(r::AbstractMatrix, T)
    
end

function multiplication_matrix(r::AbstractVector, T)
    l = length(r)
    if l > T
        throw(DomainError(""))
    end
    ℳ = zeros(eltype(r), T-l+1, T)
    for i in axes(ℳ,1), j in eachindex(r)
        ℳ[i, (i-1) + j] = r[j]
    end
    return ℳ
end

function multiplication_projection()
    
end

"""
    range_basis(A)

Find an orthonormal basis of the columns of `A`.

See also [`kernel_basis`](@ref) and [`compare_spans`](@ref)
"""
function range_basis(A::AbstractMatrix)
    F = svd(A)
    r = rank(A) # NB. This is horribly inefficient, since we can get the rank from the QR decomposition we just computed
    return F.U[:,1:r]
end

"""
    kernel_basis(A)

Find an orthonormal basis of the kernel of `A`.

See also [`range_basis`](@ref) and [`compare_spans`](@ref)
"""
function kernel_basis(A::AbstractMatrix)
    # Why is there a function to get an orthonormal basis of the kernel of a matrix, but not for its range? ¯\_(ツ)_/¯
    return nullspace(A')
end

"""
    compare_spans(P,Q)

Check of the images of two matrices with orthonormal columns are equal.
"""
function compare_spans(P, Q)
    if size(P) != size(Q)
        return false
    end
    C = P' * Q
    return !any(all(x -> isapprox(x,0.0), C[i,j] for j in axes(C,2)) for i in axes(C,1))
end