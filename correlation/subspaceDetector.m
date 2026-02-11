function [indiv_events,tabs,energyRatio,z2p,Neff,p2rms,kurt,eratio,t] = ...
    subspaceDetector(S,threshold,basisFunctionFileName,maxBasisFunctions,...
    recordLength,maxN,mpd,diffFlag,linearccnorm,plotFlag,verboseFlag,smoothFlag)
tic;
warning off signal:findpeaks:largeMinPeakHeight

if nargin < 2; threshold = 0.2; end
if nargin < 3; basisFunctionFileName = '~/research/now/reventador/casc_svd_basis_functions'; end
if nargin < 4; maxBasisFunctions = 5; end
if nargin < 5; recordLength = 120; end
if nargin < 6; maxN = 5e2; end
if nargin < 7; mpd = round(recordLength/4); end
if nargin < 8; diffFlag = []; end
if nargin < 9; linearccnorm = true; end
if nargin < 10; plotFlag = false; end
if nargin < 11; verboseFlag = false; end
if nargin < 12; smoothFlag = false; end

%%
if isempty(diffFlag)
    diffFlag = false;
end

%%
variableInfo = who('-file',basisFunctionFileName);

%%
basisFunctions = load(basisFunctionFileName);
basisFunction = basisFunctions.U;
basisFunction = double(normalizeWaveforms(basisFunction));
Fs = basisFunctions.newFs;
lfc = basisFunctions.lfc;
hfc = basisFunctions.hfc;
kstnm = basisFunctions.kstnm;
chan = basisFunctions.chan;
ntwk = basisFunctions.ntwk;

%%
if ismember('snr', variableInfo) % returns true
    snr = basisFunctions.snr;
else
    snr = [];
end

%%
if ismember('locID', variableInfo) % returns true
    locID = basisFunctions.locID;
else
    locID = [];
end

%%
clear basisFunctions;

%%
lntwk = length(ntwk);
MAXPOSSIBLEDATACHANNELS = length(kstnm);

%%
if ~exist('locID','var') || isempty(locID)
    locID = repmat("",MAXPOSSIBLEDATACHANNELS,1);
end

%%
if lntwk < MAXPOSSIBLEDATACHANNELS
    ntwk = repmat(ntwk(1),MAXPOSSIBLEDATACHANNELS,1);
end

%%
nsamples = size(basisFunction,1);
winlen = recordLength*Fs;

%%
% %snr = ones(MAXPOSSIBLEDATACHANNELS,1);
% if ~exist('snr','var') || isempty(snr)
%     snr = ones(MAXPOSSIBLEDATACHANNELS,1);
% end

%% sync waveform structure tp template list here
localSNCL = string([char(kstnm) char(chan) char(ntwk) char(locID)]);
goodI = false(MAXPOSSIBLEDATACHANNELS,1);
waveformKstnms = pull(S,'kstnm');
waveformKcmpnm = pull(S,'kcmpnm');
waveformKnetwk = pull(S,'knetwk');
waveformKhole = pull(S,'khole');
waveformSNCLs = strip(string([char(waveformKstnms(:)) char(waveformKcmpnm(:)) ...
    char(waveformKnetwk(:)) char(waveformKhole(:))]));
Ssort = S;

nn = 0; %counter for total number of sensors
for i = 1:MAXPOSSIBLEDATACHANNELS
    sncl_ = localSNCL(i);
    [lia,locb] = ismember(sncl_,waveformSNCLs);
    goodI(i) = lia;
    if lia
        nn = nn + 1;
        Ssort(nn) = S(locb);
    end
end

%%
if winlen <= nsamples
    if MAXPOSSIBLEDATACHANNELS > 1
        basisFunction = basisFunction(1:winlen,:,1:maxBasisFunctions);
    else
        basisFunction = basisFunction(1:winlen,1:maxBasisFunctions);
    end
else
    if MAXPOSSIBLEDATACHANNELS > 1
        basisFunction = basisFunction(1:end,:,1:maxBasisFunctions);
    else
        basisFunction = basisFunction(1:end,1:maxBasisFunctions);
    end
end

%%
if ~exist('snr','var') || isempty(snr)
    snr = ones(1,MAXPOSSIBLEDATACHANNELS);
end

%%
if MAXPOSSIBLEDATACHANNELS > 1
    nBF = size(basisFunction,3);
    if nBF > 1 %pages exist
        main = basisFunction;
        for pp = 1:nBF
            main_ = main(:,:,pp);
            main_ = normalizeWaveforms(flipud(detrend(main_)));
            main(:,:,pp) = main_;
        end
    else % no pages present
        main = normalizeWaveforms(flipud(detrend(basisFunction))); %detrend, flip, and normalize
    end
else
    nBF = maxBasisFunctions;
    main = normalizeWaveforms(flipud(detrend(basisFunction))); %detrend, flip, and normalize
end

%% do some preallocation
box = ones(winlen,1);
energyRatio = NaN(maxN,1);
Neff = energyRatio;
tabs = NaT(maxN,1);

%%
if MAXPOSSIBLEDATACHANNELS == 1
    indiv_events = NaN(winlen,maxN);
    z2p = energyRatio;
    p2rms = z2p;
    kurt = z2p;
    eratio = [];
    t = [];
else
    indiv_events = [];
    z2p = NaN(maxN,MAXPOSSIBLEDATACHANNELS);
    p2rms = z2p;
    kurt = z2p;
    eratio = [];
    t = [];
end

%%
noise = 0;
noiseWin = noise*Fs;
npoles = 4;
tw = 0.0002;
minThresh = 0.05;

%%
if verboseFlag
    fprintf('----------\n');
    fprintf('done with the pre-processing\n');
    fprintf('----------\n');
end

%%
nAll = 0;
tic;

%%
if ~nn %if at least 1 sensor...
    fprintf(2,'no waveform data, try again\n');
    return;
end

S = Ssort(1:nn);
clear waveformKstnms Ssort;

%%
snr_ = snr(goodI);
snr_ = snr_(:);
snr_ = snr_';

%%
tw = 200;
experimentalFlag = false;
S = detrendWaveforms(S);

if experimentalFlag
    S = interpolateWaveforms(S);
    S = filterWaveforms(S,lfc,hfc,npoles);
    S = nanGapWaveforms(S);
else
    S = resampleWaveforms(detrendWaveforms(S),Fs);
    S = syncWaveforms(S,false,true,true);
    S = filterWaveforms(S,lfc,hfc);
    S = nanGapWaveforms(S,0);
    S = detrendWaveforms(S);
    fprintf('done filtering\n')
end

%%
data = double(pull(S));
badI = ~isfinite(data);
data(badI) = 0;

ldf2 = size(data,1);
npts_ = pull(S,'npts');
t = getTimeVec(S(find(npts_ == ldf2,1))); %<-- hack
clear npts_;

%% chack that data meets a minimum length (at least as long as template)
if ldf2 < winlen
    fprintf(2,'waveforms not long enough, try again\n');
    return;
end

%%
fprintf('----------\n');
fprintf('processing day: %s\n',datestr(dateshift(t(1),'start','day')));
fprintf('----------\n');

%%
mpd = min([winlen round(mpd*Fs)]);
maxIndex = ldf2 - winlen + 1;
data2 = data.^2; %toc; disp('got data squared');

%%
norms = fftfilt(box,data2); %sum of squares for each time window
snorms = norms > 1;

%%
eratio = NaN(size(data));
if MAXPOSSIBLEDATACHANNELS > 1
    U = squeeze(main(:,goodI,1:nBF));
else
    U = main(:,i:nBF); %toc; disp('done with cc...');
end

for i = 1:nn %loop through each channel
    if MAXPOSSIBLEDATACHANNELS > 1
        disp(i)
        U_ = squeeze(U(:,i,:));
    else
        U_ = U;
    end
    cc_ = fftfilt(U_,data(:,i));
    cc_ = sum(cc_.^2,2); %sum energy controbutions across all nBF basis functions (this channel only)
    snorms_ = snorms(:,i);
    norms_ = norms(:,i);
    eratio_ = snorms_.*cc_./norms_;

    %%
    badI = ~isfinite(eratio_) | ~snorms_;
    eratio_(badI) = 0;
    eratio_(badI) = NaN;                      % toc; disp('done with removing nans');
    eratio(:,i) = eratio_;
end

%% compress ccnorm
notbadI = isfinite(eratio) & snorms;
nEffective = sum(notbadI,2);
sumSnr = sum(snr_.*notbadI,2);
weights = snr_./sumSnr; %weights are normalized

%%
if MAXPOSSIBLEDATACHANNELS > 1
    % ccnorm is a vector after this step
    eratio(~snorms) = NaN;
    if linearccnorm
        %fprintf("sum weights: %f\n",sum(weights(isfinite(weights)),"all"));
        eratio = sum(weights.*eratio,2,"omitnan");      %weighted sum, weights should (hopefully) sum to 1
    else
        eratio = median(eratio,2,"omitnan");
    end
elseif plotFlag
    figure(1);
    hold on;
    plot(eratio,'.-');
    zoom on;
end

%%
eratio(~isfinite(eratio)) = 0; % fix NaNs
if smoothFlag
    eratio = zpkFilter(eratio,-inf,1/mpd,1,1,1); %delete me!
end

%% subplot 2
if plotFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    aa_(1) = subplot(2,1,1);
    plot(aa_(1),t(1:end-winlen+1),eratio(winlen:end));
    zoom on; grid on;
    hold on;
end

%%
if isempty(threshold)
    threshold = median(eratio,"omitnan") + 8*mad(eratio,1);
end
[CC,locs_] = findpeaks(eratio,'MINPEAKHEIGHT',threshold,'MINPEAKDISTANCE',mpd); %,'Threshold',1e-4)
locs_ = locs_-winlen+1;

%%
lI = locs_ <= maxIndex & isfinite(CC) & locs_ >= noiseWin & locs_ > 0;
locs_ = locs_(lI);
CC = CC(lI);
tabs_ = t(locs_);
ll = length(locs_);
neff_ = nEffective(locs_+winlen-1);

if ~ll
    fprintf('----\n');
    fprintf('no events found\n');
    fprintf('----\n');
else
    if plotFlag
        plot(aa_(1),t(locs_),CC,'p');
    end

    %%
    fprintf('----\n');
    fprintf('found: %d event(s)\n',ll);
    fprintf('----\n');

    %%
    for j = 1:ll
        nAll = nAll + 1;
        tabs(nAll) = tabs_(j);
        Neff(nAll) = neff_(j);
        energyRatio(nAll) = CC(j);

        if MAXPOSSIBLEDATACHANNELS == 1
            tmp = data(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin);
            indiv_events(:,nAll) = tmp;
        else
            tmp = data(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin,:); %filtered data
            z2p(nAll,goodI) = max(abs(tmp),[],1);
            p2rms(nAll,goodI) = peak2rms(tmp);
            kurt(nAll,goodI) = kurtosis(tmp,0);
        end
    end
end

clear CC locs_ maxIndex minIndex lI tabs_ cc_ aI dOrig norms data2 nanI maxI snorms eratio_
if ~plotFlag
    clear data S
end
disp(' ');

%%
energyRatio = energyRatio(1:nAll);
Neff = Neff(1:nAll);
tabs = tabs(1:nAll);

%%
if plotFlag
    kk = (1:3)'; %(1:3:12)';
    if nn < length(kk)
        kk = 1;
    end

    %%
    %kk
    disp(pull(S(kk),'kstnm'));
    aa_(2) = subplot(2,1,2);
    plot(aa_(2),t,data(:,kk));

    %S2 = loadWaveforms(datetime(2022,07,31),2,["SAGA"],["HHZ"],"EC","",true,true);
    %S2 = detrendWaveforms(S2);
    %S2 = nanGapWaveforms(S2,0);
    %aa_(3) = subplot(3,1,3);
    %plot(aa_(3),getTimeVec(S2),zpkFilter(S2.d,0.25,2,100));
    title(aa_(end),S(kk(1)).kstnm);

    linkaxes(aa_,'x');
    set(aa_,'Box','off');
    set(aa_(1:end-1),'XTick',[]);
    linkaxes(aa_,'x');
    axis tight;
end

%%
if MAXPOSSIBLEDATACHANNELS == 1
    indiv_events = indiv_events(:,1:nAll);
    try
        indiv_events = detrend(indiv_events);
    catch ME
        warning(ME.message);
    end
else
    indiv_events = [];
end

%%
if MAXPOSSIBLEDATACHANNELS == 1
    try
        disp('filtering data for ancillary statistics...')
        z2p = max(abs(indiv_events),[],1)';
        p2rms = peak2rms(indiv_events)';
        kurt = kurtosis(indiv_events,0)';
    catch ME
        disp('not enough memory, opting NOT to filter...');
        z2p = max(abs(indiv_events),[],1)';
        p2rms = peak2rms(indiv_events)';
        kurt = kurtosis(indiv_events,0)';
        warning(ME.message);
    end
else
    z2p = z2p(1:nAll,:);
    p2rms = p2rms(1:nAll,:);
    kurt = kurt(1:nAll,:);
end
