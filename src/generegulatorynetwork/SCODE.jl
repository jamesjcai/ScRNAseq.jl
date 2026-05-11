using DelimitedFiles

"""
    scode(fdata, ftime; tfnum, pnum, maxite, maxB, minB) -> (A, W, B)

Infer a gene regulatory network via the SCODE algorithm.

- `fdata`: path to tab-delimited gene-by-cell expression file
- `ftime`: path to tab-delimited pseudotime file (pseudotime in column 2)
- Returns adjacency matrix `A`, weight matrix `W`, and diagonal `B`.

Reference: Matsumoto et al., Bioinformatics 2017.
"""
function scode(fdata::AbstractString, ftime::AbstractString;
               tfnum::Int=100, pnum::Int=4, maxite::Int=10,
               maxB::Float64=2.0, minB::Float64=-10.0)
    X          = readdlm(fdata, '\t')
    pseudotime = readdlm(ftime, '\t')[:, 2]
    pseudotime = pseudotime ./ maximum(pseudotime)

    cnum  = length(pseudotime)
    W     = zeros(tfnum, pnum)
    Z     = zeros(pnum, cnum)
    new_B = rand(pnum) .* (maxB - minB) .+ minB
    old_B = copy(new_B)
    RSS   = Inf

    function sample_Z!()
        for i in 1:pnum, j in 1:cnum
            Z[i, j] = exp(new_B[i] * pseudotime[j]) + (rand() * 0.002 - 0.001)
        end
    end

    for ite in 1:maxite
        target = rand(1:pnum)
        new_B[target] = rand() * (maxB - minB) + minB
        if ite == maxite
            new_B .= old_B
        end
        sample_Z!()
        ZZt = Z * Z'
        for i in 1:tfnum
            W[i, :] = ZZt \ (Z * X[i, :])
        end
        WZ      = W * Z
        tmp_RSS = sum((X .- WZ) .^ 2)
        if tmp_RSS < RSS
            RSS           = tmp_RSS
            old_B[target] = new_B[target]
        else
            new_B[target] = old_B[target]
        end
    end

    B    = Diagonal(new_B)
    invW = pinv(W)
    A    = W * B * invW
    return A, W, B
end
