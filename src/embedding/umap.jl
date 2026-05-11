import UMAP

function umap(X)
    X2 = convert(Array{Float64,2}, X)
    Y = UMAP.umap(X2, 3)
    return Y'
end
