
module Normalization
export 
norm_libsize,
norm_deseq

using ScRNAseq.DifferentialExpression

function norm_libsize(X)
    lbsz=sum(X,dims=1)
    X=(X./lbsz)*1e4;
    return X
end

function norm_deseq(X)
    X = DifferentialExpression.median_of_ratios_normalization(X);
    return X
end
end
