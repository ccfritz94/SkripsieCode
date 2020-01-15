function mse = plscv(X,Y,ncomp,cvp,mcreps,ParOptions)

[n,dx] = size(X);

% Return error for as many components as asked for; some columns may be NaN
% if ncomp is too large for CV.
mse = NaN(2,ncomp+1);

% The CV training sets are smaller than the full data; may not be able to fit as
% many PLS components.  Do the best we can.
if isa(cvp,'cvpartition')
    cvpType = 'partition';
    maxncomp = min(min(cvp.TrainSize)-1,dx);
    nTest = sum(cvp.TestSize);
else
    cvpType = 'Kfold';
%    maxncomp = min(min( floor((n*(cvp-1)/cvp)-1), dx));
    maxncomp = min( floor((n*(cvp-1)/cvp)-1), dx);
    nTest = n;
end
if ncomp > maxncomp
    warning(message('stats:plsregress:MaxComponentsCV', maxncomp));
    ncomp = maxncomp;
end

% Cross-validate sum of squared errors for models with 1:ncomp components,
% simultaneously.  Sum the SSEs over CV sets, and compute the mean squared
% error
CVfun = @(Xtr,Ytr,Xtst,Ytst) sseCV(Xtr,Ytr,Xtst,Ytst,ncomp);
sumsqerr = crossval(CVfun,X,Y,cvpType,cvp,'mcreps',mcreps,'options',ParOptions);
mse(:,1:ncomp+1) = reshape(sum(sumsqerr,1)/(nTest*mcreps), [2,ncomp+1]);
