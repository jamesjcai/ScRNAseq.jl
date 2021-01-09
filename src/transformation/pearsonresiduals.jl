function pearsonresiduals(X)
    u=(sum(X,dims=2)*sum(X,dims=1))./sum(X);
    s=sqrt.(u+(u.^2)./100);
    X=(X-u)./s;
    n=size(X,2);
    sn=sqrt(n);
    X[X.>sn].=sn;
    X[X.<-sn].=-sn;
    return X
end