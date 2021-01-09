# Pkg.add("TSne")
using Statistics, DelimitedFiles

a = readdlm("GSM3204305_P_N_Expr_999cells.csv", ',');
X = a[2:end, 2:end]
# X=a[2:500,2:500]
rescale(A; dims = 1) =
    (A .- mean(A, dims = dims)) ./ max.(std(A, dims = dims), eps())

# using TSne
# Y = tsne(X', 2, 50, 1000, 20.0);    # samples in row

using UMAP, Plots
X2 = convert(Array{Float64,2}, X)
Y = umap(X2, 3)  # ,2;n_neighbors=5);
f2 = plot(Y[1, :], Y[2, :], seriestype = :scatter)

Y = Y';
# theplot = scatter(Y[:,1], Y[:,2], marker=(2,2,:auto,stroke(0)))  # , color=Int.(allabels[1:size(Y,1)]))
f = scatter3d(
    Y[:, 1],
    Y[:, 2],
    Y[:, 3],
    marker = (2, 2, :auto, stroke(0)),
    color = Int.(sum(X, dims = 1)),
)
# Plots.pdf(f, "myplot.pdf")
