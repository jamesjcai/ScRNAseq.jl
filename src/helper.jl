rescale(A; dims=1) = (A .- mean(A, dims=dims)) ./ max.(std(A, dims=dims), eps())


# vecnorm(x) = x./norm.(x[:,i] for i in 1:size(x,2))'
vecnorm(x::AbstractMatrix) = norm.(x[:,i] for i in 1:size(x,2))
function normc!(x)
    for i in 1:size(x,2)
        x[:,i]=x[:,i]./norm(x[:,i])
    end
end

#=
function [A,G,H] = pythagoreanMeans(list)
    A = mean(list);           % arithmetic mean
    G = exp(mean(log(list))); % geometric mean
    H = 1./mean(1./list);     % harmonic mean
end
=#
