#using Pkg
#Pkg.add("DelimitedFiles")
#Pkg.add("MAT")
#Pkg.add("Arpack")

using Statistics, DelimitedFiles, LinearAlgebra, Statistics, MAT, Arpack, Pkg

cd("E:\\GitHub\\julia_test\\scrnaseq_code")
file=matopen("s1131_cr.mat")
X=read(file,"X");
s=read(file,"t_sne");
close(file)
cd("..")
X = convert(Array{Float64,2}, X)

# Pkg.activate("E:\\GitHub\\julia_test\\pcnet\\pcnet.jl")
using pcrnet

#@time A1=pcrnet.pcnetwork1(rand(200,300));
#@time A3=pcrnet.pcnetwork3(rand(200,300));
@time A=pcrnet.pcnetwork3(X');
