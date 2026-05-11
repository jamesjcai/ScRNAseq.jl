module GeneRegulatoryNetwork
export pcregress, affprop, scode

using LinearAlgebra, Statistics, Arpack, Clustering
# using CommunityDetection   # https://github.com/afternone/CommunityDetection.jl

function pcregress(X::AbstractMatrix{T}; nev::Int=3) where T<:Real
    n = size(X, 2)
    A = 1.0 .- Matrix(I, n, n)
    for k in 1:n
        y  = X[:, k]
        𝒳  = X[:, 1:end .≠ k]
        _, ϕ = Arpack.eigs(𝒳'𝒳, nev=nev, which=:LM)
        s  = 𝒳 * ϕ
        s ./= (norm.(s[:,i] for i = 1:size(s,2)) .^ 2)'
        b  = sum(y .* s, dims=1)
        𝒷  = ϕ * b'
        A[k, A[k,:] .== 1.0] = 𝒷
    end
    return A
end


function affprop(A::AbstractMatrix{T}; pref=nothing) where T<:AbstractFloat
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("A must be a square matrix ($(size(A)) given)."))
    r = isnothing(pref) ? Clustering.affprop(A) : Clustering.affprop(A; pref=pref)
    return r
end

include("SCODE.jl")

end
