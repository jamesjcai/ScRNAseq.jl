
module DifferentialExpression
export 
    de_mannwhitney,
    DESeq2


using GLM: glm, coeftable, LogLink, NegativeBinomial
using HypothesisTests:MannWhitneyUTest
using MultipleTesting
using Statistics: median
using StatsBase: geomean

# scRNASeq Differential Gene Expression with two-sample
# MannWhitneyUTest

function de_mannwhitney(X,Y)
    assert(size(X,1)==size(Y,1))
    n=size(X,1)
    p=zeros(n,1)
    for k in 1:n
        p[k]=pvalue(MannWhitneyUTest(X[k,:],Y[k,:]))
    end
    q=adjust(vec(p), BenjaminiHochberg())    
    return p,q
end


add_constant(X) = hcat(ones(size(X, 1)), X)
loge_to_log2(β) = β |> exp |> log2
median_of_ratios_normalization(X) = normalize_geomean(X) |> scale_size_factors
normalize_geomean(X) = X ./ mapslices(geomean, X, dims=1)
scale_size_factors(X) = X .* mapslices(median, X, dims=2)


struct DESeq2_results
    β
    SE
    z
    p
end


function get_DESeq2_results(models)
    # Needs to be fixed to use Log2 in the future
    β  = [coeftable(model).cols[1][2] for model in models]
    SE = [coeftable(model).cols[2][2] for model in models]
    z  = [coeftable(model).cols[3][2] for model in models]
    p  = [coeftable(model).cols[4][2] for model in models]
    DESeq2_results(β, SE, z, p)
end


function DESeq2(X, y)
    X_norm = median_of_ratios_normalization(X)
    models = mapslices(X -> glm(add_constant(X), y, NegativeBinomial(), LogLink()), X_norm, dims=1)
    get_DESeq2_results(models)
end


"""
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8#Sec2

Work in progress
----------------
Some low hanging fruit
1. Convert analysis to log2
2. Offer extra analysis/plots that DESeq2 provides

Challenging work
----------------
1. Ensuring this module behaves exactly like DESeq2
2. Going through the rest of the paper for details
"""


end
