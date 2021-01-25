module Clustering

using Clustering, ParallelKMeans

export kmeans,
       kmeanspar

function kmeans(X, k)
    Clustering.kmeans(X, k);
end

function kmeanspar(X, k)
    # https://pydatablog.github.io/ParallelKMeans.jl/stable/#How-To-Use
    ParallelKMeans.kmeans(X, k);
end

end