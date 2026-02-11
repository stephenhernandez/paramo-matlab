function [t1_commonI,t2_commonI,just_t1I,just_t2I] = ...
    synchronizeCatalog(t1,t2,thresh,strictFlag,plotFlag)
if nargin < 3
    thresh = 45;
end

if nargin < 4
    strictFlag = true;
end

if nargin < 5
    plotFlag = false;
end

%%
tAll = [t1; t2];
iAll = [zeros(size(t1)); ones(size(t2))];
[tAll,sortI] = sort(tAll);
iAll = iAll(sortI);

if plotFlag
    plotData(tAll,thresh);
    return;
end

difft = seconds(diff(tAll));
diffi = diff(iAll);
timeThreshI = difft <= thresh;
goodI1 = timeThreshI & diffi == 1;
goodIm1 = timeThreshI & diffi == -1;

sum_1 = sum(goodI1);
sum_m1 = sum(goodIm1);

if ~strictFlag
    if sum_1 >= sum_m1
        t2I = sort([find(goodI1)+1; find(goodIm1)]);
        t2_common = tAll(t2I);
        t1I = sort([find(goodI1); find(goodIm1)+1]);
        t1_common = tAll(t1I);
        disp("t2 generally after t1");
    else
        t1I = sort([find(goodIm1)+1; find(goodI1)]);
        t1_common = tAll(t1I);
        t2I = sort([find(goodIm1); find(goodI1)+1]);
        t2_common = tAll(t2I);
        disp("t1 generally after t2");
    end
else
    if sum_1 >= sum_m1
        t2I = sort(find(goodI1)+1);
        t2_common = tAll(t2I);
        t1I = sort(find(goodI1));
        t1_common = tAll(t1I);
        disp("t2 strictly after t1");
    else
        t1I = sort(find(goodIm1)+1);
        t1_common = tAll(t1I);
        t2I = sort(find(goodIm1));
        t2_common = tAll(t2I);
        disp("t1 strictly after t2");
    end
end

%%
lia = ismember(t1,t1_common);
t1_commonI = find(lia);
just_t1I = find(~lia);

lia = ismember(t2,t2_common);
t2_commonI = find(lia);
just_t2I = find(~lia);

function plotData(tAll, thresh)
figure();
semilogy(tAll(2:end),seconds(diff(tAll)),'.');
zoom on; grid on;
hold on;
plot([min(tAll) max(tAll)],[thresh thresh],'--','linewidth',2);
