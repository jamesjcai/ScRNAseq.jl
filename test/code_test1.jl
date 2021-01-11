using ScRNAseq
using DelimitedFiles
cd(dirname(@__FILE__))


X=readdlm("X.txt",',',Int16)
genelist=vec(readdlm("genelist.txt",String))

X,genelist=ScRNAseq.QualityControl.selectg(X,genelist)
X1=ScRNAseq.Transformation.pearsonresiduals(X)
X2=ScRNAseq.Normalization.norm_libsize(X)
Y=ScRNAseq.Embedding.umap(X1)


# Y2=ScRNAseq.Embedding.tsne(X)
using Plots
scatter3d(Y[:,1],Y[:,2],Y[:,3])
