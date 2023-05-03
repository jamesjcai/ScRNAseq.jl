struct SingleCellExperiment0
    X::Matrix{Int16}
    g::Vector{String}
    SingleCellExperiment0(X,g) = size(X,1)!=length(g) ? error("size unmatch") : new(X,g)
end


struct SingleCellExperiment{T<:AbstractMatrix}
    X::T
    g::Vector{String}
    s::Matrix    
    SingleCellExperiment(X,g,s) = size(X,1)!=length(g) ? error("size unmatch") : new{AbstractMatrix}(X,g,s)
end


struct SCE{T<:AbstractMatrix}
    X::T
    g::Vector{String}
    s::T
    function SCE(X::S, g::Vector{String}, s::S) where S<:AbstractMatrix
        size(X,1)==length(g) || throw(DimensionMismatch())
        new{S}(X, g, s)
    end
end


struct SCE5{T<:Any, M<:AbstractMatrix{T}}
    X::M
    g::Vector{String}
    s::M
#    function SCE2(X::S, g::Vector{String}, s::S) where S<:AbstractMatrix
#        size(X,1)==length(g) || throw(DimensionMismatch())
#        new{S}(X, g, s)
#    end
end


function SingleCellExperiment(X::AbstractMatrix,g::Vector{String})
    s=rand(size(X,2),3)
    SingleCellExperiment(X, g, s)
end


# SingleCellExperiment(X,g) = size(X,1)!=length(g) ? error("size unmatch") : new{AbstractMatrix}(X,g)



using ScRNAseq
cd("C:\\Users\\jcai\\Desktop\\refwrap_test\\forjulia")


using HDF5
f=HDF5.h5open("hugedata.mat","r")
X=read(f,"X")
close(f)
X=convert(Array{Int16,2}, X)
X=convert(Matrix{Int16},X)


f=HDF5.h5open("smalldata_Int16.mat","r")
X=read(f,"X")
close(f)

using SparseArrays
X=sparse(X)


using MAT
f=MAT.matopen("smalldata_sparse.mat")
X=read(f,"X")    # out of memory error for large full matrix
g=read(f,"g")    # not support MATLAB string, but cellstr is okay
g=convert(Matrix{String},g)
close(f)

X=convert(Matrix{Int16},X)
# g=convert(Array{String,1},g)
g=vec(g)
s=rand(3,3)
sce = SingleCellExperiment(X, g, s)

SCE2(X,vec(g),s)

#X=readdlm("X.txt",',',Int16)
#genelist=vec(readdlm("genelist.txt",String))

#using Statistics
X,g=ScRNAseq.QualityControl.selectg(X,g)
X1=ScRNAseq.Transformation.pearsonresiduals(X)
X2=ScRNAseq.Normalization.norm_libsize(X)
Y=ScRNAseq.Embedding.umap(X2)


# Y2=ScRNAseq.Embedding.tsne(X)
using Plots
scatter3d(Y[:,1],Y[:,2],Y[:,3])