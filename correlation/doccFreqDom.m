function [maxccp,plags,maxccn,nlags] = doccFreqDom(data,verboseFlag,maxN)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    verboseFlag = false;
end

if nargin < 3
    maxN = [];
end

%%
data = normalizeWaveforms(data);
[winlen,n] = size(data);
n1 = n-1;
ncol = n*n1*0.5;
if verboseFlag
    fprintf("Window Length: %d\n",winlen);
    fprintf("Number of Events: %d\n",n);
    fprintf("Number of Combos: %d\n",ncol);
end

%% pre-allocate data, get first ffts
maxccp = NaN(ncol,1);
maxccn = maxccp;
plags = maxccp;
nlags = maxccp;

%%
if mod(winlen,2)
    nfft = 2*(winlen+1);
else
    nfft = 2*winlen;
end

data = fft(data,nfft);
refBlock = conj(data);

nfft2 = nfft/2;
ci = nfft2+1;

%% start and end indices
ei = fliplr(1:n1);
ei = cumsum(ei);
si = ei+1;
si = [1 si];
si(end) = [];

indexVec = 1:n1;
if ~isempty(maxN)
    if maxN >= n
        fprintf("Oops... Youre maxN is larger than the input matrix. \n");
        fprintf("Reducing size and moving on.\n");
    else
        indexVec = 1:maxN;
    end
end

%% get the job done
for i = indexVec
    if verboseFlag
        fprintf("%d\n",i);
    end
    ni = n-i;
    ei_ = ei(i);
    si_ = si(i);

    % prep data
    refBlock = refBlock(:,2:end); % the size of refBlock is reduced by 1 column after each iteration (slow!)
    Ctmp = refBlock.*repmat(data(:,i),1,ni); % duplicate ni-times the current trace in data (slow!)

    % revert to time domain
    Ctmp = fftshift(ifft(Ctmp,[],1,"symmetric"),1);

    % positive and negative cc coeffs. and lags
    [maxccp_,plags_] = max(Ctmp);
    [maxccn_,nlags_] = min(Ctmp);

    % gather final products
    maxccp(si_:ei_) = maxccp_';
    plags(si_:ei_) = plags_';
    maxccn(si_:ei_) = -maxccn_';
    nlags(si_:ei_) = nlags_';
end

%% adjust lag indices
plags = plags-ci;
nlags = nlags-ci;

%%
if isempty(maxN)
    return;
end

nanI = isnan(maxccp);
maxccp(nanI) = [];
plags(nanI) = [];
maxccn(nanI) = [];
nlags(nanI) = [];