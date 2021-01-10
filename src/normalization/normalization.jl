
module Normalization
export 
norm_libsize

# using StatsBase

function norm_libsize(X)
    lbsz=sum(X,dims=1)
    X=(X./lbsz)*1e4;
    return X
end

end
