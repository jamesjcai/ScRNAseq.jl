module GeneRegulatoryNetwork
export pcregress


using LinearAlgebra, Statistics, Arpack

function pcregress(X)
    n=size(X,2)
    A=1.0 .-Matrix(I,n,n)
    for k in 1:n   
        y=X[:,k]
        𝒳=X[:,1:end.≠k]
        _,ϕ=Arpack.eigs(𝒳'𝒳,nev=3,which=:LM)
        s=𝒳*ϕ
        s ./=(norm.(s[:,i] for i=1:size(s,2)).^2)'
        b=sum(y.*s,dims=1)
        𝒷=ϕ*b'
        A[k,A[k,:].==1.0]=𝒷
    end
    return A
  end  

end
