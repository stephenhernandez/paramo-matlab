function [xu,yu] = removeDuplicates(x,y)
xu = unique(x);
maxxu = max(xu);
xu(end+1) = maxxu+1;
yu = y;
N = histcounts(x,xu);
NI = N > 1;             %find which elements have more than one instance (dups)
duplicated = xu(NI);    %catalog those dups
xu = xu(1:end-1);

ld = length(duplicated);
if ld
    discardIndex = cell(ld,1);
    for k = 1:ld
        disp(k)
        dupIndex = sort(find(x == duplicated(k)));
        keepIndex = dupIndex(1); %sample to keep
        discardIndex_ = dupIndex(2:end);
        discardIndex{k} = discardIndex_; %samples to discard (at least one will exist!)
        yu(keepIndex,:) = mean(y(dupIndex,:));
    end
    dI = cat(1,discardIndex{:});
    yu(dI,:) = [];
end