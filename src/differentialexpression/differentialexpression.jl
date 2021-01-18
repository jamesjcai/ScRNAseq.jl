
module DifferentialExpression
export 
de_mannwhitney

using HypothesisTests:MannWhitneyUTest
using MultipleTesting

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

end
