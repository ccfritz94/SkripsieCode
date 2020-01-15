function pca_loocv(X)

%// loop over data points 
for n=1:size(X,1)
    Xtrain = X([1:n-1 n+1:end],:);
    mu = mean(Xtrain);
    Xtrain = bsxfun(@minus, Xtrain, mu);
    [~,~,V] = svd(Xtrain, 'econ');
    Xtest = X(n,:);
    Xtest = bsxfun(@minus, Xtest, mu);

    %// loop over the number of PCs
    for j=1:min(size(V,2),25)
        P = V(:,1:j)*V(:,1:j)';        %//'
        err1 = Xtest * (eye(size(P)) - P);
        err2 = Xtest * (eye(size(P)) - P + diag(diag(P)));
        for k=1:size(Xtest,2)
            proj = Xtest(:,[1:k-1 k+1:end])*pinv(V([1:k-1 k+1:end],1:j))'*V(:,1:j)'; 
            err3(k) = Xtest(k) - proj(k);
        end

        error1(n,j) = sum(err1(:).^2);
        error2(n,j) = sum(err2(:).^2);
        error3(n,j) = sum(err3(:).^2);
    end    
end

error1 = sum(error1);
error2 = sum(error2);
error3 = sum(error3);