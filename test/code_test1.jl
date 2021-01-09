using ScRNAseq                          # obviously the tests use the ScTenifoldNet module... 
using DelimitedFiles
cd(dirname(@__FILE__))
X=readdlm("X.txt",',',Int16)
genelist=vec(readdlm("genelist.txt",String))

X,genelist=ScRNAseq.Qualitycontrol.selectg(X,genelist)
X1=ScRNAseq.Transformation.pearsonresiduals(X)
Y=ScRNAseq.Embedding.umap(X1)

# Y2=ScRNAseq.Embedding.tsne(X)
using Plots
scatter(Y[:,1],Y[:,2])
