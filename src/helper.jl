# Basic types for scRNA-seq data

using Statistics: mean, std
using LinearAlgebra: norm

"""
ScRNAseqData with metadata
Fields:
* `X`: gene-by-cell expression matrix
* `g`: gene list
* `s`: coordinates of embedding of cells
* `c`: class ID of cells
"""
struct ScRNAseqData
    X::AbstractMatrix
    g::Vector{String}
    s::Matrix{Float64}
    c::Vector{<:Integer}
    function ScRNAseqData(X::AbstractMatrix, g::Vector{String}, s::Matrix{Float64}, c::Vector{<:Integer})
        @assert size(X,1) == length(g)
        @assert size(X,2) == size(s,1)
        @assert size(X,2) == length(c)
        new(X, g, s, c)
    end
end

rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())

vecnorm(x::AbstractMatrix) = norm.(x[:,i] for i in 1:size(x,2))
function normc!(x)
    for i in 1:size(x,2)
        x[:,i] = x[:,i] ./ norm(x[:,i])
    end
end
