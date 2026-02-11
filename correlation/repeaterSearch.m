function [indiv_events,tabs,scaledCC,templateIndex,z2p,NCC,kurt,Neff,...
    ccnorm,t,ccMax] = repeaterSearch(S,thresh,templateFileName,NeventsTotal,NeventStart,...
    recordLength,maxN,diffFlag,mpd,subspaceFlag,linearccnorm,plotFlag,verboseFlag)

%%
warning off signal:findpeaks:largeMinPeakHeight

%%
if nargin < 2; thresh = 0.6; end
if nargin < 3; templateFileName = '~/regions/pedernales/ispt_280templates_1Hz4Hz_50sps'; end
if nargin < 4; NeventsTotal = 100; end
if nargin < 5; NeventStart = 1; end
if nargin < 6; recordLength = 10; end
if nargin < 7; maxN = 15e3; end
if nargin < 8; diffFlag = []; end
if nargin < 9; mpd = round(recordLength/5); end
if nargin < 10; subspaceFlag = false; end
if nargin < 11; linearccnorm = false; end
if nargin < 12; plotFlag = false; end
if nargin < 13; verboseFlag = false; end


%%
if isempty(diffFlag)
    diffFlag = false;
end

%%
if strcmp(templateFileName,'~/research/now/sangay/sangay_svd_basis_functions')
    subspaceFlag = true;
end

%%
templatesData = load(templateFileName);
if subspaceFlag
    templateWaveforms = templatesData.U;
else
    templateWaveforms = templatesData.data;
end
templateWaveforms = double(normalizeWaveforms(templateWaveforms));
Fs = templatesData.newFs;
lfc = templatesData.lfc;
hfc = templatesData.hfc;
kstnm = templatesData.kstnm;
snr = templatesData.snr;
clear templatesData;

%%
debuggingMode = 0;

%%
Nsensors = length(kstnm);
NeventStop = NeventStart+NeventsTotal-1;
nsamples = size(templateWaveforms,1);
winlen = recordLength*Fs;

%%
if ~exist('snr','var') || isempty(snr)
    snr = ones(1,Nsensors);
end

%%
if winlen <= nsamples
    if Nsensors > 1
        templateWaveforms = templateWaveforms(1:winlen,:,NeventStart:NeventStop);
    else
        templateWaveforms = templateWaveforms(1:winlen,NeventStart:NeventStop);
    end
else
    if Nsensors > 1
        templateWaveforms = templateWaveforms(1:end,:,NeventStart:NeventStop);
    else
        templateWaveforms = templateWaveforms(1:end,NeventStart:NeventStop);
    end
end

%%
if Nsensors > 1
    lP = size(templateWaveforms,3);
    if lP > 1 %pages exist
        master = templateWaveforms; %preallocate
        for pp = 1:lP
            master_ = master(:,:,pp);
            master_ = normalizeWaveforms(flipud(detrend(master_)));
            master(:,:,pp) = master_;
        end
    else % no pages present
        master = normalizeWaveforms(flipud(detrend(templateWaveforms))); %detrend, flip, and normalize
    end
else
    lP = NeventsTotal;
    master = normalizeWaveforms(flipud(detrend(templateWaveforms))); %detrend, flip, and normalize
end

%%
box = ones(winlen,1);
scaledCC = NaN(maxN,1);
templateIndex = NaN(maxN,1);
Neff = templateIndex;
tabs = NaT(maxN,1);

%%
if Nsensors == 1
    indiv_events = NaN(winlen,maxN);
else
    indiv_events = [];
    z2p = NaN(maxN,Nsensors);
    kurt = z2p;
end

%%
noise = 0;
noiseWin = noise*Fs;
n = zeros(lP,1);
npoles = 6;
NCC = scaledCC;

%%
if verboseFlag
    disp('----------');
    disp('done with the pre-processing');
    disp('----------');
end

%%
freqDomFlag = false; %when true, doesnt seem to provide much speed up vs. "time-domain" fftfilt
if freqDomFlag
    ldf = 86400*Fs;
    master = fft(master,ldf);
    box = fft(box,ldf);
end

%%
nAll = 0;

%% synchwaveform structure tp template list here
goodI = false(Nsensors,1);
waveformKstnms = pull(S,'kstnm');
Ssort = S;
nn = 0;
for i = 1:Nsensors
    kstnm_ = kstnm(i);
    [lia,locb] = ismember(kstnm_,waveformKstnms);
    goodI(i) = lia;
    if lia
        nn = nn + 1;
        Ssort(nn) = S(locb);
    end
end

%%
ccMax = [];
if nn
    if nn ~= sum(goodI)
        disp('something went wrong, the two numbers dont match, fix.');
        disp([nn sum(goodI)]);
        return;
    end
    
    %%
    S = Ssort(1:nn);
    clear waveformKstnms Ssort;
    
    %%
    snr_ = snr(goodI);
    sumSnr = sum(snr_);
    weights = snr_/sumSnr;
    weights = weights(:)';  % row vector
    
    %%
    S = detrendWaveforms(S);
    dOrig = double(pull(resampleWaveforms(S,Fs)));
    S = filterWaveforms(S,lfc,hfc,npoles);
    
    %%
    if diffFlag
        S = differentiateWaveforms(S);
    end
    S = resampleWaveforms(S,Fs);
    S = padWaveforms(S);
    
    %%
    data = double(pull(S));
    ldf2 = size(data,1);
    npts_ = pull(S,'npts');
    t = getTimeVec(S(find(npts_ == ldf2,1))); %<-- hack
    clear npts_;
    
    %% check that data meets a minimun length (at least as long as template)
    if ldf2 > winlen
        %%
        if verboseFlag
            disp('----------');
            disp(['processing day: ',datestr(dateshift(t(1),'start','day'))]);
            disp('----------');
        end
        
        %%
        mpd = min([winlen/2 round(mpd*Fs)]);
        maxIndex = ldf2 - winlen + 1;
        data2 = data.^2;
        
        %%
        if freqDomFlag
            data = fft(data,ldf);
            data2 = fft(data2,ldf);
            norms = ifft(box.*data2,'symmetric');
        else
            norms = fftfilt(box,data2);
        end
        norms = sqrt(abs(norms));
        snorms = norms > 1;
        
        %%
        if plotFlag
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            aa_ = getCustomAxesPos(lP+1,1,0.12);
        end
        
        %% loop through each desired plantilla (template) or basis function (for svd option)
        for i = 1:lP
            if freqDomFlag
                if Nsensors > 1
                    ccnorm = ifft(master(:,goodI,i).*data,'symmetric');
                else
                    ccnorm = ifft(master(:,i).*data,'symmetric');
                end
            else
                if Nsensors > 1
                    ccnorm = fftfilt(master(:,goodI,i),data);
                else
                    ccnorm = fftfilt(master(:,i),data);
                end
            end
            
            %% normalize the raw cross-correlation
            ccnorm = snorms.*ccnorm./norms;
            
            %% find wird samples, set to NaN
            badI = ~isfinite(ccnorm) & ~snorms;
            ccnorm(badI) = NaN;
            
            %% subspace experiment
            if subspaceFlag
                ccnorm = ccnorm.^2;
            end
            
            %% compress ccnorm
            if nn > 1
                % ccnorm is a vector after this step
                ccnorm(~snorms) = NaN;
                ccnorm = weights.*ccnorm;
                if linearccnorm
                    ccnorm = nanmean(ccnorm,2);
                else
                    ccnorm = nanmedian(ccnorm,2);
                end
            end
            
            %% get new threshes
            nEffective = sum(snorms,2);
            
            %% smooth the discriminator
            ccnorm(~isfinite(ccnorm)) = 0; % fix NaNs
            ccnorm = abs(hilbert(ccnorm));
            ccnorm = zpkFilter(ccnorm,-inf,1/round(mpd)/1,1,1,1);
            
            %%
            if debuggingMode
                figure();
                tmpax(1) = subplot(311);
                plot(t(1:end-winlen+1),ccnorm(winlen:end)); zoom on; grid on;
                tmpax(2) = subplot(312);
                plot(t(1:end-winlen+1),nEffective(winlen:end)); grid on;
            end
            
            %%
            if nn > 0
                w_ = zeros(nn,1);
                stdtmp = w_;
                for mm = 1:nn
                    wI = nEffective == mm;
                    if sum(wI)
                        if linearccnorm
                            stdtmp_ = std(ccnorm(wI));
                        else
                            stdtmp_ = mad(ccnorm(wI),1);
                        end
                        stdtmp(mm) = stdtmp_;
                        w_(mm) = sum(wI)/ldf2;
                        ccnorm(wI) = (ccnorm(wI) - nanmedian(ccnorm(wI)))/stdtmp(mm);
                    end
                end
                wSTD = nansum(w_.*stdtmp);% weighted sum of all `nn' stds
            else
                ccnorm = ccnorm - nanmedian(ccnorm);
                if linearccnorm
                    wSTD = std(ccnorm);
                else
                    wSTD = mad(ccnorm,1);
                end
                ccnorm = ccnorm/wSTD;
            end
            
            %% new scaled threshold
            newThresh = thresh;
            if thresh < 1
                newThresh = newThresh/wSTD;
            end
            
            if debuggingMode
                tmpax(3) = subplot(313);
                plot(t(1:end-winlen+1),ccnorm(winlen:end)); zoom on; grid on;
                linkaxes(tmpax,'x');
            end
            
            %%
            ccMax = [ccMax ccnorm];
            [CC,locs_] = findpeaks(ccnorm,'MINPEAKHEIGHT',newThresh,'MINPEAKDISTANCE',mpd); %,'Threshold',1e-4);
            locs_ = locs_-winlen+1;
            
            %%
            %disp(['Max Index: ',num2str(maxIndex)]);
            lI = locs_ <= maxIndex & isfinite(CC) & locs_ >= noiseWin & locs_ > 0;
            locs_ = locs_(lI);
            CC = CC(lI);
            tabs_ = t(locs_);
            ll = length(locs_);
            neff_ = nEffective(locs_+winlen-1);
            
            %%
            if ll
                %%
                if debuggingMode
                    hold(tmpax(3),'on');
                    plot(tmpax(3),t(locs_),CC,'p'); grid on;
                end
                
                %%
                if plotFlag
                    aa_(i) = subplot(lP+1,1,i);
                    plot(t(1:end-winlen+1),ccnorm(winlen:end)); zoom on; grid on;
                    hold on;
                    plot(t(locs_),CC,'p');
                end
                
                %%
                disp('------------------------------------');
                disp(['template ',num2str(i),': ',num2str(ll),' event(s)']);
                disp('------------------------------------');
                
                %%
                for j = 1:ll
                    n(i) = n(i) + 1;
                    nAll = nAll + 1;
                    
                    %%
                    tabs(nAll) = tabs_(j);
                    scaledCC(nAll) = CC(j);
                    NCC(nAll) = CC(j)*stdtmp(neff_(j));
                    Neff(nAll) = neff_(j);
                    templateIndex(nAll) = i;
                    
                    %%
                    if Nsensors == 1
                        tmp = dOrig(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin); % unfiltered data
                        %tmp = data(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin,:); % filtered data
                        indiv_events(:,nAll) = detrend(tmp);
                    else
                        tmp = data(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin,:); % filtered data
                        z2p(nAll,goodI) = max(abs(detrend(tmp)),[],1); %RHS should be a row vector
                        kurt(nAll,goodI) = kurtosis(tmp,0); %RHS should be a row vector
                    end
                end
            else
                if verboseFlag
                    disp('------------------------------------')
                    disp(['template ',num2str(i),': no events found']);
                    disp('------------------------------------');
                end
                
                %%
                if plotFlag
                    aa_(i) = subplot(lP+1,1,i);
                    plot(t(1:end-winlen+1),ccnorm(winlen:end)); zoom on;
                    hold on;
                end
            end
        end
        clear pks_ locs_ maxIndex minIndex lI tabs_ aI dOrig norms data2 nanI maxI snorms
        if ~plotFlag
            clear data
        end
        disp(' ');
    end
else
    disp('no data');
end

%%
if plotFlag
    kk = 1;
    if nn < kk
        kk = 1;
    end
    
    %%
    aa_(lP+1) = subplot(lP+1,1,lP+1);
    plot(t,data(:,kk));
    zoom on;
    title(aa_(end),S(kk).kstnm);
    linkaxes(aa_,'x');
    set(aa_,'Box','off');
    set(aa_(1:end-1),'XTick',[]);
    linkaxes(aa_,'x');
end

%%
scaledCC = scaledCC(1:nAll);
NCC = NCC(1:nAll);
Neff = Neff(1:nAll);
tabs = tabs(1:nAll);
templateIndex = templateIndex(1:nAll);

%%
if Nsensors == 1
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
if Nsensors == 1
    try
        ieF = detrend(cumsum(zpkFilter(taper(detrend(diff(indiv_events)),0.05),lfc,hfc,Fs,npoles)));
        z2p = max(abs(ieF))';
        kurt = kurtosis(ieF,0)';
        indiv_events = ieF;
        clear ieF
    catch ME
        z2p = max(abs(ieF))';
        warning(ME.message);
    end
else
    z2p = z2p(1:nAll,:);
    kurt = kurt(1:nAll,:);
end


%%
[tabs,sI] = sort(tabs);
z2p = z2p(sI,:);
Neff = Neff(sI);
scaledCC = scaledCC(sI);
templateIndex = templateIndex(sI);
NCC = NCC(sI);
kurt = kurt(sI,:);
