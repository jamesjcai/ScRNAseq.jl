# Basic types for scRNA-seq data

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
    g::Array{String,1}
    s::Array{Float64,2}
    c::Array{Integer,1}
end


# Constructs a ScRNAseqData 
function ScRNAseqData(X::AbstractMatrix, g::String, s::Array{Float64,2}, c::Array{Integer,1})
    # a=ScRNAseqData(rand(3,4),["g1","g2","g3"],[4.8 3.2 1.2 2.3; 1.0 2.0 3.0 4.0],[3, 1, 2, 1])    
    @assert size(X,1)==length(g)
    @assert size(X,2)==size(s,1)
    @assert size(X,2)==length(c)
    new(X, g, s, c)
end

# using TableReader

rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())


# vecnorm(x) = x./norm.(x[:,i] for i in 1:size(x,2))'
vecnorm(x::AbstractMatrix) = norm.(x[:,i] for i in 1:size(x,2))
function normc!(x)
    for i in 1:size(x,2)
        x[:,i]=x[:,i]./norm(x[:,i])
    end
end

#=
function [A,G,H] = pythagoreanMeans(list)
    A = mean(list);           % arithmetic mean
    G = exp(mean(log(list))); % geometric mean
    H = 1./mean(1./list);     % harmonic mean
end
=#
