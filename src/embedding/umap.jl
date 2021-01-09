import UMAP, Plots

function umap(X)
    X2 = convert(Array{Float64,2}, X)
    Y = UMAP.umap(X2, 3)  # ,2;n_neighbors=5);
    f2 = Plots.plot(Y[1, :], Y[2, :], seriestype = :scatter)

    Y = Y';
    # theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)))  # , color=Int.(allabels[1:size(Y,1)]))
    #=
    f = scatter3d(
        Y[:, 1],
        Y[:, 2],
        Y[:, 3],
        marker = (2, 2, :auto, stroke(0)),
        color = Int.(sum(X, dims = 1)),
    )
    =#
    return Y
end


