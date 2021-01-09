
# using MarketMatrix
# https://imgur.com/a/y3C0Vd2
println(pwd())
cd(dirname(@__FILE__))
println(pwd()) 

using MatrixMarket
A=mmread("Ydf_matrix.mtx")
using UnicodePlots
spy(A)