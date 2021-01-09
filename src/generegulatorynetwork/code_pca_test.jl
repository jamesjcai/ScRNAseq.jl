using LinearAlgebra, Statistics, MAT, Arpack

# cd("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\julia_test\\pcnet")
# file=matopen("testdata.mat")
# X=read(file,"X")
# A0=read(file,"A")
# close(file)


# https://discourse.julialang.org/t/how-to-get-the-principal-components-and-variances-from-multivariatestats/15843/4

"""
    pca(data)
perform PCA using SVD
inputs:
    - data: M x N matrix of input data. (M dimensions, N trials)
outputs:
    - PC: each column is a principle component
    - V: M x 1 matrix of variances
"""
function pca(data::Array{T,2}) where T
    X = data .- mean(data, dims=2)    
    Y = X' ./ sqrt(T(size(X,2)-1))
    U,S,PC = svd(Y)
    S = diagm(0=>S)
    V = S .* S
    # find the least variance vector
    indexList = sortperm(diag(V); rev=true)
    # PCs = map(x->PC[:,x], indexList)
    return PC, diag(V)[indexList]
end
X=randn(Float64, (20,4))
pc0,d0=pca(collect(X'))


using MultivariateStats
p=fit(PCA,X')
pc1=p.proj
d1=p.prinvars

function pca2(X; k::Int=3)
    X=X.-mean(X,dims=1)
    Σ = X'X./(size(X,1)-1)              # Covariance natrix
    # D,V = eigen(Σ,sortby=-)                 # Factorise into Σ = U * diagm(S) * V'
    D,V = eigen(Σ,sortby=x -> -abs(x))                 # Factorise into Σ = U * diagm(S) * V'
    # sortby = x -> -abs(x)
    Xrot = X*V                      # Rotate onto the basis defined by U
    # pvar = sum(D[1:k]) / sum(D)  # Percentage of variance retained with top k vectors
    # X̃ = Xrot[:,1:k]                 # Keep top k vectors
    return V, D
end
pc2,d2=pca2(X)

function pca3(X, k::Int=3)
    X=X.-mean(X,dims=1)
    X=X./sqrt(size(X,1)-1)
    # var(X,dims=1)
    F = svd(X)                   # Factorise into Σ = U * diagm(S) * V'
    # Xrot = X*F.V               # Rotate onto the basis defined by U
    # pvar = sum(F.S[1:k]) / sum(F.S)  # Percentage of variance retained with top k vectors
    # X̃ = Xrot[:,1:k]                 # Keep top k vectors
    return F.V, F.S.*F.S
end
pc3,d3=pca3(X)

function pca4(X, k::Int=3)
    X=X.-mean(X,dims=1)
    X=X./sqrt(size(X,1)-1)
    # var(X,dims=1)
    F = svds(X;nsv=3)[1]                   # Factorise into Σ = U * diagm(S) * V'    
    # Xrot = X*F.V               # Rotate onto the basis defined by U
    # pvar = sum(F.S[1:k]) / sum(F.S)  # Percentage of variance retained with top k vectors
    # X̃ = Xrot[:,1:k]                 # Keep top k vectors
    return F.V, F.S.*F.S
end
pc4,d4=pca4(X)

[d0 d1 d2 d3 d4] 

```
pc0
pc1
pc2
pc3
```
