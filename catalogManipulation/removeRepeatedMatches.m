function [t,cc,removeIndices,keepIndices] = removeRepeatedMatches(t,cc,timeThresh)
% assumes t is sorted
% t must also be either datenum or datetime format

%%
[t,sI] = sort(t);
cc = cc(sI);

%%
difft = diff(t);
if isduration(difft)
    dTI = difft <= seconds(timeThresh);
else
    dTI = difft <= (timeThresh/86400);
end

loopMax = 100; %up to 100 familes in a single dataset?
removeIndices = cell(loopMax,1);
keepIndices = removeIndices;
n = 0;
while sum(dTI) && n < loopMax
    n = n+1;

    %%
    ldTI = sum(dTI);
    findI = find(dTI);
    removeIndex = NaN(ldTI,1);
    keepIndex = removeIndex;
    for i = 1:ldTI
        o1 = findI(i);
        o2 = o1+1;
        if cc(o1) >= cc(o2)
            removeIndex(i) = o2;
            keepIndex(i) = o1;
        else
            removeIndex(i) = o1;
            keepIndex(i) = o2;
        end
    end
    removeIndices{n} = removeIndex;
    keepIndices{n} = keepIndex;
    t(removeIndex) = [];
    cc(removeIndex) = [];

    % do it all over again
    difft = diff(t);
    if isduration(difft)
        dTI = difft <= duration(00,00,timeThresh);
    else
        dTI = difft <= (timeThresh/86400);
    end
end

%%
removeIndices = removeIndices(1:n);
keepIndices = keepIndices(1:n);
