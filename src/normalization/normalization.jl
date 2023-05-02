module Normalization
using Statistics

export 
norm_libsize
# norm_deseq


function norm_libsize(X)
    lbsz=sum(X,dims=1)
    γ=mean(lbsz)
    X=(X./lbsz)*γ
    return X
end

# Not sure how to include another function in the same package -- pls help!
#
# include("../differentialexpression/differentialexpression.jl")
# using .DifferentialExpression

# function norm_deseq(X)
#    X = DifferentialExpression.median_of_ratios_normalization(X);
#    return X
# end

end
