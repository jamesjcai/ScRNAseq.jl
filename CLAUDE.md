# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

ScRNAseq.jl is a Julia package for exploratory single-cell RNA-seq analysis. It wraps external Julia libraries (UMAP.jl, TSne.jl, Clustering.jl, etc.) into a cohesive pipeline with a single data container type.

## Commands

```julia
# Install dependencies (first time)
using Pkg; Pkg.instantiate()

# Run tests
include("test/runtests.jl")

# Run a single test file
include("test/code_test1.jl")        # full pipeline demo (requires data files in test/)

# Build docs
include("docs/make.jl")
```

## Core Data Type

All analysis operates on `ScRNAseqData` (defined in [`src/helper.jl`](src/helper.jl)):

```julia
struct ScRNAseqData
    X::AbstractMatrix  # gene-by-cell expression matrix (sparse)
    g::Vector          # gene names
    s::Matrix          # 2D/3D cell embedding coordinates
    c::Vector          # cell cluster/class IDs
end
```

The constructor enforces that `size(X, 1) == length(g)`.

## Module Architecture

The main module ([`src/ScRNAseq.jl`](src/ScRNAseq.jl)) includes eight submodules, each in its own subdirectory:

| Submodule | File | Key exports |
|---|---|---|
| FileIO | `src/fileio/fileio.jl` | `readtxt`, `readmtx`, `readmat`, `readhdf`, `read10xh5`, `readgenelist` |
| QualityControl | `src/qualitycontrol/qualitycontrol.jl` | `selectg`, `scstats`, `emptyrate` |
| Transformation | `src/transformation/pearsonresiduals.jl` | `pearsonresiduals` |
| Normalization | `src/normalization/normalization.jl` | `norm_libsize` |
| Embedding | `src/embedding/` | `umap`, `tsne` |
| Clustering | `src/clustering/clustering.jl` | `kmeans`, `kmeanspar` |
| DifferentialExpression | `src/differentialexpression/differentialexpression.jl` | `de_mannwhitney`, `DESeq2` |
| GeneRegulatoryNetwork | `src/generegulatorynetwork/` | `pcregress`, `affprop`, `SCODE` |

## Conventions

- Functions follow a pipeline style: take `ScRNAseqData`, return modified data or a result matrix.
- Embedding and clustering submodules are thin wrappers — the real logic lives in UMAP.jl, TSne.jl, Clustering.jl.
- Pearson residuals are capped at ±√n (n = number of cells) for numerical stability.
- DE testing uses Benjamini-Hochberg correction via MultipleTesting.jl.
- `MAST` (in DifferentialExpression) and `SCODE` (in GeneRegulatoryNetwork) are partially implemented stubs — do not rely on them.
- HDF5 reading (`readhdf`) is stubbed out with a TODO.
- Sparse matrices (SparseArrays) are used throughout; avoid dense operations on `X`.
