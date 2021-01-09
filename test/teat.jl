using DelimitedFiles
cd(dirname(@__FILE__))
X=readdlm("X.txt",',',Int16)
genelist=vec(readdlm("genelist.txt",String))

using Statistics
glogmean=log10.(mean(X,dims=2));
glogvar=log10.(var(X,dims=2));
dropoutrate=mean(X.==0,dims=2);
using SparseArrays
nnz(sparse(X))
count(!iszero,X)


using Plots
Plots.scatter(glogmean,glogvar)
Plots.scatter(glogmean,1 .-dropoutrate)

using Distributions
# fit_mle(NegativeBinomial,X[3,:].+0.0)

fit_mle(Normal,X[3,:].+0.0)
# fit(Normal,X[3,:].+0.0)
# using GLM

# http://naobioml.blogspot.com/2017/02/how-to-fit-count-data-with-negative.html

data=X[3,:].+0.0;
function f2(x::Vector)
    sum = 0
    for i in data
        sum += log(pdf(NegativeBinomial(x[1], x[2]), i))
    end
    return -1*sum
end

using Optim, Distributions

lower = [0.0, 0.0]
upper = [Inf, 1]
initial_x = [0.5, 0.5]
# x2 = optimize(DifferentiableFunction(f2), initial_x, lower, upper, Fminbox(), optimizer = GradientDescent)
x2 = optimize(f2, [0.0, 0.0], [Inf, 1],[0.5, 0.5])


function f1(x)
    sum = 0
    for i in data
        sum += log(pdf(Poisson(x), i))
    end
    -1*sum
end
x1 = optimize(f1, 0.0, 7)



#=
 ------------

using Distributions, Optim, StatsPlots
julia> data = vcat([0 for i in 1:70],
           [1 for i in 1:38],
           [2 for i in 1:17],
           [3 for i in 1:10],
           [4 for i in 1:9],
           [5 for i in 1:3],
           [6 for i in 1:2],
           [7 for i in 1:1]);
julia> function f2(x)
           sum = 0.0
           for i in data
               sum += log(pdf(NegativeBinomial(x[1], x[2]), i))
           end
           return -sum
       end
f2 (generic function with 1 method)
julia> opt_result = optimize(f2, [0.0, 0.0], [Inf, 1],[0.5, 0.5])
 * Status: success
 * Candidate solution
    Minimizer: [1.02e+00, 4.72e-01]
    Minimum:   2.224372e+02
 * Found with
    Algorithm:     Fminbox with L-BFGS
    Initial Point: [5.00e-01, 5.00e-01]
 * Convergence measures
    |x - x'|               = 0.00e+00 ≤ 0.0e+00
    |x - x'|/|x'|          = 0.00e+00 ≤ 0.0e+00
    |f(x) - f(x')|         = 0.00e+00 ≤ 0.0e+00
    |f(x) - f(x')|/|f(x')| = 0.00e+00 ≤ 0.0e+00
    |g(x)|                 = 7.79e-08 ≰ 1.0e-08
 * Work counters
    Seconds run:   0  (vs limit Inf)
    Iterations:    4
    f(x) calls:    226
    ∇f(x) calls:   226
julia> histogram(data, normalize = true, label = "Data", alpha = 0.5, linecolor = "white"); plot!(NegativeBinomial(opt_result.minimizer[1], opt_result.minimizer[2]), label = "Negative Binomial fit")

=#