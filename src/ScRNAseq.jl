### ScRNAseq.jl
###
### A julia package for scRNA-seq exploratory data analysis
##
### This file is NOT YET a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

module ScRNAseq

export
    FileIO,
    QualityControl,
    Transformation,
    Normalization,
    Embedding,
    GeneRegulatoryNetwork

    include("helper.jl")
    include("fileio/fileio.jl")
    include("qualitycontrol/qualitycontrol.jl")
    include("transformation/transformation.jl")
    include("normalization/normalization.jl")
    include("embedding/embedding.jl")
    include("differentialexpression/differentialexpression.jl")
    include("generegulatorynetwork/generegulatorynetwork.jl")

using .FileIO
using .QualityControl
using .Transformation
using .Normalization
using .Embedding
using .DifferentialExpression
using .GeneRegulatoryNetwork

end  # module ScRNAseq