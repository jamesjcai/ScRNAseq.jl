module QualityControl

using Statistics, SparseArrays

export selectg, 
       scstats, 
       emptyrate

function ismtgene(genelist)
    startswith.(uppercase.(genelist),"MT-")
end

function selectg(X,genelist)
    ng=size(X,1);
    i=vec(sum(!iszero,X,dims=2)./ng.>0.05)
    X=X[i,:];
    genelist=genelist[i];
    return X,genelist
end

function scstats(X)
    logmean    = log10.(max.(mean(X, dims=2), eps()))
    logvar     = log10.(max.(var(X, dims=2), eps()))
    dropoutrate = mean(X .== 0, dims=2)
    return logmean, logvar, dropoutrate
end

function emptyrate(X)
    count(!iszero, X) / length(X)
end

end