
using LinearAlgebra, Statistics, MAT, Arpack, MultivariateStats, Random

# p=fit(PCA,X')

# cd("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\julia_test\\pcnet")
# file=matopen("testdata.mat")
# X=read(file,"X")
# A0=read(file,"A")
# close(file)
rng = MersenneTwister(1234);
X=randn(rng,Float64,(200,300));

function pcnetwork1(X)
  n=size(X,2)
  A=1.0 .-Matrix(I,n,n)
  for k in 1:n   
      y=X[:,k]
      𝒳=X[:,1:end.≠k]
      F=svds(𝒳,nsv=3)[1]
      ϕ=F.V
      s=𝒳*ϕ
      s ./= (norm.(s[:,i] for i=1:size(s,2)).^2)'
      b=sum(y.*s,dims=1)    
      𝒷=ϕ*b'
      A[k,A[k,:].==1.0]=𝒷
  end
  return A
end


function pcnetwork2(X)
  n=size(X,2)
  A=1.0 .-Matrix(I,n,n)
  for k in 1:n   
      y=X[:,k]
      𝒳=X[:,1:end.≠k]
      F=svd(𝒳)
      ϕ=F.V[:,1:3]
      s=𝒳*ϕ
      s ./=(norm.(s[:,i] for i=1:size(s,2)).^2)'
      b=sum(y.*s,dims=1)    
      𝒷=ϕ*b'
      A[k,A[k,:].==1.0]=𝒷
  end
  return A
end


function pcnetwork3(X)
  n=size(X,2)
  A=1.0 .-Matrix(I,n,n)
  for k in 1:n   
      y=X[:,k]
      𝒳=X[:,1:end.≠k]
      # _,v=eigen(𝒳'𝒳,sortby=-)
      # v=eigvecs(𝒳'𝒳,sortby=-)
      # v=v[:,1:3]
      _,ϕ=eigs(𝒳'𝒳,nev=3,which=:LM)
      s=𝒳*ϕ
      s ./=(norm.(s[:,i] for i=1:size(s,2)).^2)'
      b=sum(y.*s,dims=1)
      𝒷=ϕ*b'
      A[k,A[k,:].==1.0]=𝒷
  end
  return A
end


function pcnetwork4(X)
  n=size(X,2)
  A=1.0 .-Matrix(I,n,n)
  for k in 1:n   
      y=X[:,k]
      𝒳=X[:,1:end.≠k]
      p=fit(PCA,𝒳')
      v=p.proj
      v=v[:,1:3]
      s=𝒳*v
      s=s./(norm.(s[:,i] for i=1:size(s,2)).^2)' 
      b=sum(y.*s,dims=1)    
      𝒷=v*b'
      A[k,A[k,:].==1.0]=𝒷
  end
  return A
end


function pcnetwork5(X)
  # http://hua-zhou.github.io/teaching/biostatm280-2017spring/slides/16-eigsvd/eigsvd.html
  n=size(X,2)
  A=1.0 .-Matrix(I,n,n)
  for k in 1:n   
      y=X[:,k]
      𝒳=X[:,1:end.≠k]
      # _,v=eigen(𝒳'𝒳,sortby=-)
      # v=eigvecs(𝒳'𝒳,sortby=-)
      # v=v[:,1:3]
      U,S,V=svd(𝒳)
      ϕ=V[:,1:3]
      b=V*inv(diagm(S))*U'*y
      𝒷=ϕ*b'
      A[k,A[k,:].==1.0]=𝒷
  end
  return A
end

@time A1=pcnetwork1(X);
@time A2=pcnetwork2(X);
@time A3=pcnetwork3(X);
# @time A5=pcnetwork5(X);
# @time A4=pcnetwork4(X);
A1≈A2≈A3
# A1≈A5