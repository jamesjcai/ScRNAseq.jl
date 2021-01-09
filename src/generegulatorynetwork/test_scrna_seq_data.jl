println("hello")
using DelimitedFiles
a=readdlm("GSM3204305_P_N_Expr_999cells.csv",',');
b=a[2:end,2:end]




using CSV, DataFrames, Statistics
# CSV.File("GSM3204305_P_N_Expr_999cells.csv"; datarow=2)
df = CSV.File("GSM3204305_P_N_Expr_999cells.csv"; datarow=2) |> DataFrame!
X=convert(Matrix, df[:,2:end])
typeof(X)

# df = DataFrame(x = rand(3),w=rand(3))
# dv = @data([NA, 3, 2, 5, 4])
# mean(dv)
# b=readtable("GSM3204305_P_N_Expr_999cells.csv")

libsize=sum(X,dims=1)
libsize[libsize.>20000]
filter(x -> x > 20000, libsize)

a = [1 2; 3 4]
a[a .== 1]
a[[false true; false true]]

X[[X[:,1].>=0,libsize.>20000]]

