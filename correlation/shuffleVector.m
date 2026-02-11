function [numShifts,ncol,excess,shuffleVec,flipVec] = shuffleVector(n)
n1 = n - 1;
ncol = n*n1*0.5;
numShifts = ceil(n1*0.5);
excess = n*numShifts - ncol;
shuffleVec = [];
flipVec = 1;

%%
if n < 3
    return;
end

%%
if n < 4
    shuffleVec = [1; 3; 2];
    flipVec = 2;
    return;
end

%%
shuffleVec = ones(n,numShifts);
shuffleVec(1,:) = 1:numShifts;

%%
normalBlockN = n-numShifts;
nextInc = n1;
for k = 2:normalBlockN
    shuffleVec(k,:) = shuffleVec(k-1,:)+nextInc;
    nextInc = nextInc - 1;
end

ns1 = numShifts-1;
for k = normalBlockN+1:n-1
    shuffleVec(k,ns1+1) = n-ns1-1;
    shuffleVec(k,1:ns1) = shuffleVec(k-1,ns1+1)+(1:ns1);
    ns1 = ns1 - 1;
end

cumdiffsum = cumsum((n1:-1:n1-numShifts)');
for k = 1:numShifts
    x_ = [n1; n1+cumdiffsum];
    shuffleVec(end-k+1:end,k) = x_(1:k);
    n1 = n1 - 1;
end
shuffleVec = shuffleVec(1:ncol)';

%%
n1 = n - 1;
if n < 5
    flipVec = 3;
elseif n < 6
    flipVec = [3; 4; 7];
elseif n < 7
    flipVec = [4; 5; 9];
else
    lastIndex = cumsum((n1:-1:1)');
    flipWidth = numShifts;
    if excess
        flipWidth = flipWidth - 1;
    end
    nFlips = sum(1:flipWidth);
    flipVec = zeros(nFlips,1);
    ei = cumsum(flipWidth:-1:1);
    si = 1;
    nn = 1;
    for i = flipWidth:-1:1
        flipVec(si:ei(nn)) = (lastIndex(nn)-i+1:lastIndex(nn))';
        si = ei(nn)+1;
        nn = nn + 1;
    end
end
