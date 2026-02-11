function [shifted_data,maxccp_,G,plags_,raw_shifts,meancc,medcc] = ...
    apply_vdcc(indiv_events,ccData,weightedFlag,plotFlag,verboseFlag,CANUSEGPU)
if nargin < 6
    CANUSEGPU = false;
end

if nargin < 5
    verboseFlag = false;
end

if nargin < 4
    plotFlag = false;
end

if nargin < 3
    weightedFlag = false;
end

if nargin < 2
    ccData = [];
end

n = size(indiv_events,2);
if isempty(ccData) || size(ccData,2) < 3
    if verboseFlag
        fprintf("Performing Cross-Correlations\n");
    end

    indiv_events = normalizeWaveforms(indiv_events);
    if size(ccData,2) < 2
        if CANUSEGPU
            indiv_events = gpuArray(single(indiv_events));
        end
        [maxccp_,plags_,~] = doccFreqCircShiftPolarities(indiv_events,verboseFlag);
    else
        maxccp_ = ccData(:,1);
        plags_ = ccData(:,2);
    end
    [~,~,flipper] = svd(indiv_events,"econ");
    flipper = sign(flipper(:,1))';
else
    maxccp_ = ccData(:,1);
    plags_ = ccData(:,2);
    pol_ = ccData(:,3);
    pol = squareform(pol_);
    pol(pol==0) = 1;
    flipper = sign(sum(pol));
    sumFlipper = sum(flipper);

    if sumFlipper == 0 || sumFlipper == -n
        flipper = pol(1,:);
    end

    if sum(flipper == 0)
        flipper = pol(1,:);
    end

    if any(flipper < 0)
        fprintf("some elements were flipped\n");
        disp(flipper')
        find(flipper < 0)
    end
end

%%
indiv_events = double(indiv_events);
indiv_events = indiv_events.*flipper;

%%
plags_ = double([plags_; 0]);

%%
G = Gvdcc(n);
if weightedFlag
    n2 = size(G,1);
    W = [maxccp_; 0.1];
    %W = W/sum(W);
    W = sparse(W').*speye(n2);
    %raw_shifts = ((G'*W*G)^(-1))*G'*W*plags_;
    raw_shifts = lscov(G,plags_,W);
else
    W = (1/n)*speye(n);
    raw_shifts = W*G'*plags_;
end

%%
raw_shifts = -round(raw_shifts);

%%
if verboseFlag
    fprintf("Applying Shifts\n");
end

%%
shifted_data = apply_shifts(indiv_events,raw_shifts);

%%
if plotFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(raw_shifts,'.'); zoom on; grid on;
end

%%
meancc = mean(maxccp_);
medcc = median(maxccp_);