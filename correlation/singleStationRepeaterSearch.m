function [indiv_events,tabs,pks,templateIndex,maxAmpRMS,pksOrig,saveStds,Neff] = ...
    singleStationRepeaterSearch(thresh,templateFileName,NeventsTotal,NeventStart,recordLength,dayStart,dayStop,maxN,diffFlag,saveFlag,saveName,rawDataDir)
tic;
warning off signal:findpeaks:largeMinPeakHeight

%%
if nargin < 1; thresh = 0.6; end
if nargin < 2; templateFileName = '~/regions/pedernales/ispt_280templates_1Hz4Hz_50sps'; end
if nargin < 3; NeventsTotal = 100; end
if nargin < 4; NeventStart = 1; end
if nargin < 5; recordLength = 10; end
if nargin < 6; dayStart = datetime(2008,07,01); end
if nargin < 7; dayStop = datetime(2017,07,01); end
if nargin < 8; maxN = 15e3; end
if nargin < 9; diffFlag = []; end
if nargin < 10; saveFlag = []; end
if nargin < 11; saveName = []; end
if nargin < 12; rawDataDir = []; end

%%
if isempty(diffFlag)
    diffFlag = false;
end

if isempty(saveFlag)
    saveFlag = false;
end

if isempty(saveName)
    saveName = '~/regions/pedernales/ispt_repeater_search_2008_2016_100TemplatesV1_50sps_forFutureStacking';
end

if isempty(rawDataDir)
    rawDataDir = '~/rawdata/';
end
%'~/regions/pichincha/ggpc_repeater_search_2008_2016_100TemplatesV1_50sps_forFutureStacking'; end

%%
% 01-Jan-2020 06:20:41.195
% 01-Jan-2020 07:28:51.495
% 01-Jan-2020 09:47:55.395
% 01-Jan-2020 11:44:24.515
% 01-Jan-2020 20:15:44.055

%%
templatesData = load(templateFileName);
data = templatesData.data;
data = double(normalizeWaveforms(data));
Fs = templatesData.newFs;
lfc = templatesData.lfc;
hfc = templatesData.hfc;
kstnm = templatesData.kstnm;
chan = templatesData.chan;
ntwk = templatesData.ntwk;
snr = templatesData.snr;
locID = "";
clear templatesData;
linearccnorm = true;
plotFlag = false;

%%
Nsensors = length(kstnm);
NeventStop = NeventStart+NeventsTotal-1;
[nsamples,lP] = size(data);
winlen = recordLength*Fs;

if winlen <= nsamples
    if Nsensors > 1
        data = data(1:winlen,:,NeventStart:NeventStop);
    else
        data = data(1:winlen,NeventStart:NeventStop);
    end
else
    if Nsensors > 1
        data = data(1:end,:,NeventStart:NeventStop);
    else
        data = data(1:end,NeventStart:NeventStop);
    end
end

%%
%%
if ~exist('snr','var')
    snr = ones(1,Nsensors);
end

%%
if Nsensors > 1
    lP = size(data,3);
    if lP > 1 %pages exist
        master = data;
        for pp = 1:lP
            master_ = master(:,:,pp);
            master_ = normalizeWaveforms(flipud(detrend(master_)));
            master(:,:,pp) = master_;
        end
    else % no pages present
        master = normalizeWaveforms(flipud(detrend(data))); %detrend, flip, and normalize
    end
else
    master = normalizeWaveforms(flipud(detrend(data))); %detrend, flip, and normalize
end

%
% figure(); plot(master); zoom on;
%

box = ones(winlen,1);
pks = NaN(maxN,1);
templateIndex = NaN(maxN,1);
Neff = templateIndex;
tabs = NaT(maxN,1);
dayUsed = NaT(maxN,1);
if Nsensors == 1
    indiv_events = NaN(winlen,maxN);
else
    indiv_events = [];
    maxAmpRMS = NaN(maxN,Nsensors);
end
noise = 0;
noiseWin = noise*Fs;
n = zeros(lP,1);
days = dayStart:dayStop;
lDays = length(days);
dayNum = 0;
npoles = 4;
Hd = zpkOperator(lfc,hfc,Fs,npoles);

%%
saveStds = NaN(lDays,lP);
pksOrig = pks;
minThresh = 0.35;
%maxThresh = 0.95;

disp(' ');
disp('done with the pre-processing');
disp(' ');

%%
derivedRate = false;
verboseFlag = false;

freqDomFlag = false; %when true, doesnt seem to provide much speed up vs. "time-domain" fftfilt
if freqDomFlag
    ldf = 86400*Fs;
    master = fft(master,ldf);
    box = fft(box,ldf);
end

%%
nAll = 0;
for i = 1:lDays
    tic;
    disp(datestr(days(i)));
    
    %% load day data here
    nn = 0;
    S = populateWaveforms([Nsensors 1]);
    goodI = false(Nsensors,1);
    for jj = 1:Nsensors
        S_ = loadWaveforms(days(i),1,kstnm(jj),chan(jj),ntwk(jj),locID,derivedRate,verboseFlag,rawDataDir);
        if ~isnat(S_.ref)
            nn = nn + 1;
            S(nn,1) = resampleWaveforms(S_,Fs);
            goodI(jj) = true;
        end
    end
    
    %%
    snr_ = snr(goodI);
    sumSnr = sum(snr_);
    weights = snr_/sumSnr;
    
    %%
    if nn
        S = S(1:nn);
        S = differentiateWaveforms(S);
        S = detrendWaveforms(S);
        if ~diffFlag
            %S = differentiateWaveforms(S);
            S = taperWaveforms(S,0.004);
            S = intWaveforms(S);
        end
        S = padWaveforms(S);
        
        dayNum = dayNum + 1;
        dayUsed(dayNum) = days(i);
        t = getTimeVec(S(1));       %toc; disp('got time vector');
        data = double(pull(S));     %toc; disp('data extracted');
        dOrig = data;
        data = filter(Hd,data);     %toc; disp('done filtering');
        ldf2 = size(data,1);
        
        if ldf2 > winlen
            mpd = min([winlen ldf2])/2;
            maxIndex = ldf2 - winlen + 1;
            data2 = data.^2; %toc; disp('got data squared');
            
            %% experimental
            %if mod(ldf,2)
            %    ldf = ldf + 1;
            %end
            
            %%
            if freqDomFlag
                data = fft(data,ldf);
                data2 = fft(data2,ldf);
                norms = ifft(box.*data2,'symmetric');
            else
                norms = fftfilt(box,data2); %toc; disp('done with sum of squares...');
            end
            norms = sqrt(abs(norms)); %toc; disp('done with square root of sum of squares...');
            snorms = norms > 1; %toc; disp('done with mask...');
            
            %%
            for jj = 1:lP
                if freqDomFlag
                    if Nsensors > 1
                        ccnorm = ifft(master(:,goodI,jj).*data,'symmetric');
                    else
                        ccnorm = ifft(master(:,jj).*data,'symmetric');
                    end
                else
                    if Nsensors > 1
                        ccnorm = fftfilt(master(:,goodI,jj),data);
                    else
                        ccnorm = fftfilt(master(:,jj),data); %toc; disp('done with cc...');
                    end
                end
                ccnorm = snorms.*ccnorm./norms;     % toc; disp('done with normalized cc...');
                nanI = ~isfinite(ccnorm);
                ccnorm(nanI) = 0;                   % toc; disp('done with removing nans');
                maxI = ccnorm >= 1;
                ccnorm(maxI) = 0;                   % toc; disp('done with removing nans');
                
                %
                if nn > 1
                    % ccnorm is a vector after this step
                    ccnorm(~snorms) = NaN;
                    ccnorm = weights'.*ccnorm;
                    if linearccnorm
                        ccnorm = nansum(ccnorm,2);
                    else
                        ccnorm = nanmedian(ccnorm,2);
                    end
                    ccnorm(~isfinite(ccnorm)) = 0;  % wherever ccnorm is not finite, set to 0
                end
                
                % get new threshes
                nEffective = sum(snorms,2);
                if Nsensors > 1
                    w_ = NaN(nn,1);
                    stdtmp = NaN(nn,1);
                    for mm = 1:nn
                        wI = nEffective == mm;
                        if linearccnorm
                            stdtmp(mm) = std(ccnorm(wI));
                        else
                            stdtmp(mm) = mad(ccnorm(wI),1);
                        end
                        w_(mm) = sum(wI)/ldf2;
                        ccnorm(wI) = ccnorm(wI)/stdtmp(mm);
                    end
                    stds_ = sum(w_.*stdtmp);
                else
                    stds_ = std(ccnorm);
                end
                saveStds(dayNum,jj) = stds_;
                
                newThresh_ = thresh;
                %if thresh < 1
                %    newThresh_ = thresh;
                %else
                %    newThresh_ = thresh*stds_;
                %end
                
                %fprintf('original threshold is: %f\n',newThresh_);
                if Nsensors == 1
                    newThresh_ = max([minThresh newThresh_]);
                end
                %newThresh_ = min([maxThresh newThresh_]);
                %fprintf('new threshold is: %f\n',newThresh_);
                
                [pks_,locs_] = findpeaks(ccnorm,'MINPEAKHEIGHT',newThresh_,'MINPEAKDISTANCE',mpd,'Threshold',1e-4);
                locs_ = locs_-winlen+1;
                
                lI = locs_ <= maxIndex & isfinite(pks_) & locs_ >= noiseWin & locs_ > 0;
                locs_ = locs_(lI);
                pks_ = pks_(lI);
                tabs_ = t(locs_);
                ll = length(locs_);
                neff_ = nEffective(locs_);
                
                if pks_
                    if plotFlag
                        kstnm_ = kstnm(goodI);
                        for kk = 6%1:nn
                            figure('units','normalized','outerposition',[0 0 1 1]);
                            aa_(1) = subplot(211);
                            plot(t(1:end-winlen+1),ccnorm(winlen:end)); zoom on;
                            hold on;
                            plot(t(locs_),pks_,'p');
                            
                            aa_(2) = subplot(212);
                            plot(t,data(:,kk)); zoom on; title(kstnm_(kk))
                            linkaxes(aa_,'x');
                            suptitle(['Template Number: ',num2str(jj)])
                        end
                    end
                    
                    disp('----');
                    disp(['template ',num2str(jj),': ',num2str(ll),' event(s)']);
                    disp('----');
                    
                    for j = 1:ll
                        n(jj) = n(jj) + 1;
                        nAll = nAll + 1;
                        tabs(nAll) = tabs_(j);
                        pks(nAll) = pks_(j)/stds_;
                        pksOrig(nAll) = pks_(j);
                        if Nsensors == 1
                            tmp = dOrig(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin);
                            indiv_events(:,nAll) = tmp;
                        else
                            tmp = data(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin,:);
                            maxAmpRMS(nAll,goodI) = rssq(tmp)';
                        end
                        templateIndex(nAll) = jj;
                        Neff(nAll) = neff_(j);
                    end
                else
                    disp('----')
                    disp(['template ',num2str(jj),': no events found']);
                    disp('----');
                end
            end
            clear pks_ locs_ maxIndex minIndex lI tabs_ t ccnorm ccnorm_ aI dOrig S norms data2 nanI maxI snorms data
            disp(' ');
        end
    else
        disp('day is empty, skipping');
    end
end

%%
dayUsed = dayUsed(1:dayNum);
saveStds = saveStds(1:dayNum,:);
pks = pks(1:nAll);
pksOrig = pksOrig(1:nAll);
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
        ieF = zpkFilter(taper(detrend(diff(indiv_events)),0.05),lfc,hfc,Fs,6);
        maxAmpRMS = rms(ieF)';
        indiv_events = ieF;
        clear ieF
    catch ME
        maxAmpRMS = rms(indiv_events)';
        warning(ME.message);
    end
else
    maxAmpRMS = maxAmpRMS(1:nAll,:);
end

%%
if saveFlag
    disp('saving mat file');
    tic;
    save(saveName,'-v7.3');
    toc;
end
