function [maxccp,plags,maxccn,nlags] = doccFreqCircShift(data,verboseFlag,refBlock,maxN,maxLag)
%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
%

%%
if nargin < 2
    verboseFlag = true;
end

if nargin < 3
    refBlock = [];
end

if nargin < 4
    maxN = [];
end

%%
data = data./rssq(data);
[winlen,n] = size(data);
if nargin < 5
    maxLag = winlen-1;
end

%% n == 1
if n < 2
    if verboseFlag
        disp('single vector. doing nothing.');
    end
    maxccp = NaN(1);
    maxccn = maxccp;
    plags = maxccp;
    nlags = maxccp;
    return;
end

%%
if verboseFlag
    fprintf("Window Length: %d\n",winlen);
    fprintf("Number of Events: %d\n",n);
end

%% get proper fft length
if mod(winlen,2)
    nfft = 2*(winlen+1); %if odd, add 1, multiply by two (is now EVEN)
else
    nfft = 2*winlen;
end

%%
data = fft(data,nfft);
nfft2 = nfft/2;
zeroIndex = nfft2+1;

%% get the job done
if n < 3
    %% n == 2
    if verboseFlag
        disp('Number of Combos: 1');
    end
    data(:,2) = conj(data(:,2));
    data = prod(data,2);
    data = fftshift(ifft(data,[],1,'symmetric'),1);

    [maxccp,plags] = max(data);
    [maxccn,nlags] = min(data);
    maxccn = -maxccn;

    %% adjust lag indices
    plags = plags-zeroIndex;
    nlags = nlags-zeroIndex;
    return;
end

%% if no user defined refBlock, make our own [default behaviour, experimental otherwise]
if isempty(refBlock)
    if isempty(maxN)
        refBlock = conj(data);
    end
else
    %experimental section, populate here
end
refBlock = circshift(refBlock,-1,2);

%% n is 3 or greater (n >= 3)
[numShifts,ncol,excess,shuffleVec,flipVec] = shuffleVector(n);
if verboseFlag
    disp(['Number of Combos: ',num2str(ncol)]);
end
maxccp = NaN(ncol,1);
maxccn = maxccp;
plags = maxccp;
nlags = maxccp;

%%
if numShifts > 1
    %% n >= 4 (numShifts >= 2)
    ei = n*(1:numShifts)';
    si = 1 + n*(0:numShifts-1)';

    for i = 1:numShifts-1
        if verboseFlag
            fprintf("%d %d\n",i,numShifts);
        end
        ei_ = ei(i);
        si_ = si(i);
        [maxccp_,plags_,maxccn_,nlags_] = docc_(refBlock,data,zeroIndex,maxLag);
        refBlock = circshift(refBlock,-1,2);
        maxccp(si_:ei_) = maxccp_';
        plags(si_:ei_) = plags_';
        maxccn(si_:ei_) = -maxccn_';
        nlags(si_:ei_) = nlags_';
    end

    %%
    si_ = si(numShifts);
    ei_ = ei(numShifts);
    if excess
        ei_ = ei_ - excess;
        [maxccp_,plags_,maxccn_,nlags_] = docc_(refBlock(:,1:end-excess),data(:,1:end-excess),zeroIndex,maxLag);
    else
        [maxccp_,plags_,maxccn_,nlags_] = docc_(refBlock,data,zeroIndex,maxLag);
    end
    maxccp(si_:ei_) = maxccp_';
    plags(si_:ei_) = plags_';
    maxccn(si_:ei_) = -maxccn_';
    nlags(si_:ei_) = nlags_';

    %% adjust lag indices
    maxccp(shuffleVec) = maxccp;
    maxccn(shuffleVec) = maxccn;
    plags(shuffleVec) = plags;
    nlags(shuffleVec) = nlags;
    plags(flipVec) = -plags(flipVec);
    nlags(flipVec) = -nlags(flipVec);
else
    %% n == 3 (numShifts < 2)
    [maxccp,plags,maxccn,nlags] = docc_(refBlock,data,zeroIndex,maxLag);

    %% adjust lag indices
    maxccp(shuffleVec) = maxccp;
    maxccn(shuffleVec) = -maxccn;
    plags(shuffleVec) = plags;
    nlags(shuffleVec) = nlags;
    plags(flipVec) = -plags(flipVec);
    nlags(flipVec) = -nlags(flipVec);

    %%
    maxccp = maxccp';
    maxccn = maxccn';
    plags = plags';
    nlags = nlags';
end
end

%%
function [maxccp,plags,maxccn,nlags] = docc_(refBlock,data,zeroIndex,maxLag)
Ctmp = refBlock.*data;
Ctmp = fftshift(ifft(Ctmp,[],1,'symmetric'),1);

si = zeroIndex - maxLag;
ei = zeroIndex + maxLag;
Ctmp = Ctmp(si:ei,:);

%% positive and negative cc coeffs. and lags
zeroIndex = maxLag+1;
[maxccp,plags] = max(Ctmp);
plags = plags - zeroIndex;
[maxccn,nlags] = min(Ctmp);
nlags = nlags - zeroIndex;
end
