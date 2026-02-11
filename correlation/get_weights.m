function [weights,Sxx] = get_weights(Sxx,lambda,variance,nIter)
if nargin < 4
    nIter = 3;
end

%%
nDims = ndims(Sxx);
sizeSxx = size(Sxx);
nTapers = sizeSxx(end);
if nTapers < 2
    fprintf("only 1 taper, doing nothing\n");
    return;
end
defaultW = 1/nTapers;
weights = defaultW*ones(sizeSxx);

%%
indices = repmat({':'},1,nDims-1);
Sest_xx = mean(Sxx(indices{:},1:2),nDims);
iter = 0;
while iter < nIter
    for k = 1:nTapers
        lambda_k = lambda(k);
        sqrt_lambda = sqrt(lambda_k);
        lambda_k_S = sqrt_lambda.*Sest_xx;
        weights_ = lambda_k_S./(sqrt_lambda*lambda_k_S + (1-lambda_k)*variance);
        %weights_(~isfinite(weights_) | weights_ <= 0) = defaultW;
        weights(indices{:},k) = weights_;
    end
    normWeights = 1./sum(weights,nDims);
    weights = normWeights.*weights;
    Sest_xx = sum(weights.*Sxx,nDims);
    iter = iter + 1;
end
Sxx = Sest_xx;