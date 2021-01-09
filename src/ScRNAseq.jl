### ScRNAseq.jl
###
### A julia package for scRNA-seq exploratory data analysis
##
### This file is NOT YET a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

module ScRNAseq

export
    Fileio,
    Qualitycontrol,
    Transformation,
    Embedding

    include("helper.jl")
    include("fileio/fileio.jl")
    include("qualitycontrol/qualitycontrol.jl")
    include("transformation/transformation.jl")
    include("embedding/embedding.jl")

using .Fileio
using .Qualitycontrol
using .Transformation
using .Embedding

end  # module ScRNAseq