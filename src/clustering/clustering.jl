module Clustering

using Clustering

export kmeans,
       affprop

function kmeans(X, k)
    Clustering.kmeans(X, k);
end

function affprop(X, k)
    Clustering.affprop(X, k);
end

end