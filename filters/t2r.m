function [rate,varargout] = t2r(t,timeWidth,amp,causalFlag)
if nargin < 4
    causalFlag = true;
end

if nargin < 3
    amp = [];
end

if nargin < 2
    timeWidth = days(1);
end

if isduration(timeWidth)
    if causalFlag
        tAll = [t; t + timeWidth - timeWidth/1000];
    else
        tAll = [t; t - timeWidth + timeWidth/1000];
    end

    %%
    iAll = [ones(size(t)); zeros(size(t))];
    [~,sI] = sort(tAll);
    iAll = iAll(sI);

    %
    posI = iAll == 1;

    %
    cum1 = cumsum(iAll);
    cum2 = cumsum(~iAll);
    newSum2 = cum1 - cum2;
    if causalFlag
        rate = newSum2(posI);
    else
        rate = -newSum2(posI)+1;
    end
else
    N = floor(timeWidth);
    fprintf("Using fixed N: %d\n",N);
    if N < 2
        fprintf(2,"Error. N must be greater than 1\n");
        return;
    end
    lt = length(t);
    if lt < N
        fprintf(2,"Length of time vector is too small\n");
        return;
    end
    fprintf("rate is expressed as number of events per day\n");
    A = (1:N)';
    B = N*ones(lt-N,1);
    [~,si,ei] = cutWindows(t,N,N-1,false);
    dur_ = days(t(ei)-t(si));
    if causalFlag
        rate = [A; B];
        dur = [days(t(1:N-1)-t(1)); dur_];
        %rate = rate./dur;
        %rate(1) = 0;
    else
        rate = [B; flipud(A)];
        dur = [dur_; days(t(end)-t(end-N+2:end))];
        %rate = rate./dur;
        %rate(end) = 0;
    end
end

%%
if isempty(amp)
    if ~isduration(timeWidth)
        if causalFlag
            rate = rate./dur;
            rate(1) = 0;
        else
            rate = rate./dur;
            rate(end) = 0;
        end
    end
    return;
end

uniqueRates = unique(rate);
lUniqueRates = length(uniqueRates);
meanMagsFixedTimeWin = NaN(length(t),1);
medianMagsFixedTimeWin = NaN(length(t),1);
sumEnergyFixedTimeWin = NaN(length(t),1);

%%
for i = 1:lUniqueRates
    rate_ = uniqueRates(i);
    indexThisUniqueRate = find(rate == rate_);
    nThisUniqueRate = length(indexThisUniqueRate);
    theseMags = NaN(rate_,nThisUniqueRate);
    for j = 1:nThisUniqueRate
        if rate_ == 1
            ampI = indexThisUniqueRate(j);
        else
            if causalFlag
                ampI = indexThisUniqueRate(j)-rate_+1:indexThisUniqueRate(j);
            else
                ampI = indexThisUniqueRate(j):indexThisUniqueRate(j)+rate_-1;
            end
        end
        theseMags(:,j) = amp(ampI);
    end

    meanMagsFixedTimeWin(indexThisUniqueRate) = mean(theseMags,1)';
    medianMagsFixedTimeWin(indexThisUniqueRate) = median(theseMags,1)';
    %sumEnergyFixedTimeWin(indexThisUniqueRate) = sum(theseMags.^2,1); %
    sumEnergyFixedTimeWin(indexThisUniqueRate) = sum(10.^(1.5*theseMags + 4.8),1)';
end

%%
varargout(1) = {meanMagsFixedTimeWin};
varargout(2) = {medianMagsFixedTimeWin};
varargout(3) = {sumEnergyFixedTimeWin};

if ~isduration(timeWidth)
    if causalFlag
        rate = rate./dur;
        rate(1) = 0;
    else
        rate = rate./dur;
        rate(end) = 0;
    end
end