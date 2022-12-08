

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
    hankel_matrix!(ℋ, W, L)
    return ℋ
end

function hankel_matrix(W::AbstractVector{T}, L) where {T<:Union{AbstractVector,AbstractMatrix}}
    return hcat(hankel_matrix.(W,L)...)
end

function hankel_matrix(W::AbstractVector, L)
    return hankel_matrix(W', L)
end

function hankel_matrix!(ℋ, W, L)
    q, T = size(W)
    for j in 1:T-L+1, i in 1:q, t in 1:L
        ℋ[(t-1)*q+i,j] = (W[i,(j-1)+t])
    end
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

"""
    multiplication_matrix(r, T, q=1)

Compute the multiplication matrix ℳ_T associated to `r`.

The result depends on the shape of `r`.
- If `r` is a row or column vector, the result is a `T-length(r) × T` matrix whose consists of shifted copies of `r`. If `q` is greater than 1, `r` is shifted by `q` elements each time.
- If `r` is a matrix, the multiplication matrices of each row are concatenated along their columns.

## Examples
```jldoctest
julia> multiplication_matrix(1:5, 5)
1×5 Matrix{Int64}:
 1  2  3  4  5

julia> multiplication_matrix(1:5, 6)
2×6 Matrix{Int64}:
 1  2  3  4  5  0
 0  1  2  3  4  5

julia> multiplication_matrix(1:4, 4, 2)
3×8 Matrix{Int64}:
 1  2  3  4  0  0  0  0
 0  0  1  2  3  4  0  0
 0  0  0  0  1  2  3  4

julia> multiplication_matrix([1 2 3 4 5; 6 7 8 9 10], 6)
4×6 Matrix{Int64}:
 1  2  3  4   5   0
 0  1  2  3   4   5
 6  7  8  9  10   0
 0  6  7  8   9  10
```
"""
function multiplication_matrix(r::AbstractMatrix, T, q=1)
    g = size(r,1)
    l = div(size(r,2), q) -1
    if l >= T
        throw(DomainError("T must be larger than l. Got T=$T and l=$l."))
    end
    ℳ = zeros(eltype(r), g*(T-l), q*T)
    for i in axes(r,1)
        multiplication_matrix!(view(ℳ,((i-1)*(T-l)) .+ (1:T-l), :), r[i,:], T, q, l)
    end
    return ℳ
end

function multiplication_matrix(r::AbstractVector, T, q=1)
    l = div(length(r), q) -1
    if l >= T
        throw(DomainError("T must be larger than l. Got T=$T and l=$l."))
    end
    ℳ = zeros(eltype(r), T-l, q*T)
    multiplication_matrix!(ℳ, r, T, q, l)
    return ℳ
end

"""
    multiplication_matrix!(ℳ, r, T, q, l)

Compute the multiplication matrix of length `T` associated with vector `r` and write it to `ℳ`

For `q>1`, `r` is assumed to be of length `(l+1)*q` and the matrix is constructed with blocks of length `q`.
"""
function multiplication_matrix!(ℳ, r::AbstractVector, T, q, l)
    ℳ .= zero(eltype(r))
    for i in axes(ℳ,1), j in 1:l+1
        ℳ[i, ((i+j-2)*q) .+ (1:q)] = r[((j-1)*q) .+ (1:q)]
    end
end

"""
    multiplication_projection(D, p, q)

Compute the orthogonal projection of `D` on the space of multiplication matrices, assuming blocks of size `p × q`.
"""
function multiplication_projection(D::AbstractMatrix, p, q)
    T = div(size(D,2), q)
    @assert size(D,2) == q*T
    l = T - div(size(D,1), p)
    @assert size(D,1) == p*(T-l)
    R = zeros(p, q*(l+1))
    for i in 1:p
        multiplication_projection!(
            view(R, i, :), 
            view(D, (i-1)*(T-l) .+ (1:T-l), :), 
            q, 
            T, 
            l
        )
    end
    return multiplication_matrix(R, T, q)
end

function multiplication_projection(D::AbstractMatrix, q)
    T = div(size(D,2), q)
    @assert size(D,2) == q*T
    l = T - size(D,1)
    @assert l < T
    R = zeros(q*(l+1))
    multiplication_projection!(R, D, q, T, l)
    return multiplication_matrix(R, T, q)
end

"""
    multiplication_projection!(R, D, q, T, l)

Compute the orthogonal projection of `D` on the space of multiplication matrices and write it to `R`.

`D` is assumed to be of size `T-l × q*T` and `R` is assumed to be size `(l+1)*q`.
"""
function multiplication_projection!(R, D, q, T, l)
    cols = 1:(l+1)*q
    for t in 1:T-l
        R .+= D[t, (t-1)*q .+ cols]
    end
    R ./= (T-l)
end

"""
    range_basis(A)

Find an orthonormal basis of the columns of `A`.

See also [`kernel_basis`](@ref) and [`compare_spans`](@ref)
"""
function range_basis(A::AbstractMatrix)
    F = svd(A)
    r = _rank(F) 
    return F.U[:,1:r]
end

"""
    kernel_basis(A)

Find an orthonormal basis of the left kernel of `A`.

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

# Evil Type piracy
function _rank(F::SVD; atol=1e-10, rtol=1e-10)
    r = 0
    ϵ = max(atol, rtol*minimum(size(F)) * F.S[1])
    for σ in F.S
        if σ < ϵ
            return r
        else
            r += 1
        end
    end
    return r
end

function low_rank_approximation(A, r)
    F = svd(A)
    if _rank(F) > r
        return F.U[:,1:r] * Diagonal(F.S[1:r]) * F.Vt[1:r,:]
    else
        return A
    end
end