using Documenter
using ScRNAseq

makedocs(
    sitename = "ScRNAseq",
    format = Documenter.HTML(),
    modules = [ScRNAseq]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
