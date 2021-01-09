import TSne, Plots

function tsne(X)
    X = convert(Array{Float64,2}, X);
    Y = TSne.tsne(X, 3);
 
    return Y
end


