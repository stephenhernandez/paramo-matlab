function [indiv_events,tabs,maxAmp,indiv_events_orig,maxAmpRMS,snr,...
    meanFreq,medFreq,bandwidth,family,linear_stack,flipped_indices,maxccp,maxccn,...
    plags,nlags,optsToUse] = waveformSimilarity(varargin)

% Input Name - variable type,default value
% stnm - s,'RETU'
% chan - s,'SHZ'
% dayStart
% dayEnd
% lfc - f,-inf
% hfc - f,4
% thresh - f,0.65
% sta - f,5
% lta - f,10
% mph - f,1.5
% secDur - f,10
% noiseWin - f,1
% volcano - s,'tungurahua'
% diffFlag - l,false
% npoles - i,4
% hFlag - l,true
% zeroPhase - l,false
% N - i,smooth parameter (odd)
% vlines - f,time marks
% locID - s, ''
% net - s, 'EC'
% maxSNR - f,5e3, [maximum SNR allowed,test variable,could be useful in screening for glitches...]

%% function defaults, deal variables
nVarargin = length(varargin);
functionDefaults = {"RETU","SHZ",datetime(2015,04,10),datetime(2015,04,10),0.5,8,0.65,15,15,1.5,12,2,...
    'tungurahua',false,4,false,false,51,{datetime(2015,04,10,00,21,47),datetime(2015,04,10,3,40,00),...
    datetime(2015,04,10,5,15,00),datetime(2015,04,10,17,00,00),datetime(2015,04,11,13,15,00),...
    datetime(2015,04,12,12,00,00)},"","EC",5e5};
optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[stnm,chan,dayStart,dayEnd,lfc,hfc,thresh,sta,lta,mph,secDur,noiseWin,volcano,...
    diffFlag,npoles,hFlag,zeroPhase,N,vlines,locID,net,maxSNR] = deal(optsToUse{:});

%% get various plotting directories and strings
pFormat = '-djpeg';
days = dayStart:dayEnd;
[yyyy,month,dayOne] = datevec(dayStart);

dayTmp = day(days,'doy');
doy1= dayTmp(1);
doyStr = padCalendarDay(doy1);
preFixStr0 = strcat('~/recent_plots/waveform_similarity/',volcano,'/',char(stnm),...
    '/',num2str(yyyy),'/',doyStr);
unix(['\rm -rf ',preFixStr0]); %remove base directories
preFixStr1 = strcat(preFixStr0,'/uniq_fams');
preFixStr2 = strcat(preFixStr0,'/uniq_fams/time');
preFixStr3 = strcat(preFixStr0,'/summary');

%re-create base directories
mkdir(preFixStr1);
mkdir(preFixStr2);
mkdir(preFixStr3);

%% preallocate
nEventsfound = 0;
maxEventNumber = 2e4;
snr = NaN(maxEventNumber,1);
tabs = dayStart; %datetime(yyyy,month,daysOrig(1));
tabs = repmat(tabs,maxEventNumber,1);
refTime = dayStart; %datetime(yyyy,month,daysOrig(1));
interpFlag = true;

%% get basic info
S = loadWaveforms(datetime(yyyy,month,dayOne),1,stnm,chan,net,locID);
Fs = round(1/S.delta);
if Fs ~= 100
    Fs = 100;
end

if diffFlag
    winlen = secDur*Fs + 1;
else
    winlen = secDur*Fs;
end
noiseWin = noiseWin*Fs;
indiv_events = zeros(winlen,maxEventNumber);
indiv_events_orig = zeros(winlen,maxEventNumber);

%% find events
disp('finding events');
for i = 1:length(days)
    disp(days(i));
    [yyyy_,mm_,dd_] = datevec(days(i));
    %[S,gap_stops] = loadSacData(stnm,chan,yyyy_,mm_,dd_,net,locID,interpFlag);
    %[S,gap_stops] = loadSacData(datetime(yyyy_,mm_,dd_),1,stnm,chan,net,locID,interpFlag);
    S = loadWaveforms(datetime(yyyy_,mm_,dd_),1,stnm,chan,net,locID,interpFlag);

    if ~isnat(S.ref)
        FsNow = round(1/S(1).delta);
        if FsNow ~= 100
            disp('resampling sac data')
            S = resampleWaveforms(S,100);
        end
        if diffFlag
            S = differentiateWaveforms(S);
        end
        d = S.d;
        
        ldf = length(d);
        df = detrend(d);
        
        
        if any(isfinite([lfc hfc]))
            Sf = filterWaveforms(S,lfc,hfc,npoles,[],zeroPhase);
            df = Sf.d;
        end
        
        t = getTimeVec(S);
        
        tic;
        [locs_,snr_] = stalta(Sf,sta,lta,mph,hFlag,0,0,1);
        maxIndex = ldf - winlen + 1;
        minIndex = noiseWin+1;
        lI = locs_ <= maxIndex & locs_ > minIndex & snr_ < maxSNR;
        locs_ = locs_(lI);
        snr_ = snr_(lI);
        toc;
        
        tabs_ = t(locs_);
        if ~isempty(locs_)
            disp(length(locs_))
            for j = 1:length(tabs_)
                % populate snr and time vectors
                nEventsfound = nEventsfound + 1;
                snr(nEventsfound) = snr_(j);
                tabs(nEventsfound) = tabs_(j);
                
                % cut events using unfiltered waveforms
                tmp = d(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin);
                tmp2 = df(locs_(j)-noiseWin:locs_(j)+winlen-1-noiseWin);
                tmp = demean(tmp);
                tmp = detrend(tmp);
                tmp = demean(tmp);
                tmp2 = demean(tmp2);
                tmp2 = detrend(tmp2);
                tmp2 = demean(tmp2);
                indiv_events_orig(:,nEventsfound) = tmp;
                indiv_events(:,nEventsfound) = tmp2;
            end
        end
    else
        disp(['skipping day: ',datestr(datetime(yyyy_,mm_,dd_))]);
        continue
    end
end

%% clear temporary variables and truncate oversized variables
clear i j d tmp ldf snr_ locs_ tabs_ ref ref_ maxIndex minIndex t df lI;
%nfft = nextpow2(winlen);
nfft = winlen; %2^nfft;
indiv_events = indiv_events(:,1:nEventsfound);
indiv_events_orig = indiv_events_orig(:,1:nEventsfound);
tabs = tabs(1:nEventsfound);
snr = snr(1:nEventsfound);
indiv_events2 = indiv_events;
indiv_events = zeros(nfft,nEventsfound);
indiv_events(1:winlen,:) = indiv_events2;
indiv_events2 = indiv_events_orig;
indiv_events_orig = zeros(nfft,nEventsfound);
indiv_events_orig(1:winlen,:) = indiv_events2;
clear indiv_events2;

%% set aside unfiltered copy, get unfiltered amplitudes
winlen = size(indiv_events,1);
maxAmp = max(abs(indiv_events))';   %maxAmp = max(abs(indiv_events_orig))';
maxAmpRMS = rms(indiv_events)';
plotMinAmp = 0;                     %floor(min(maxAmpRMS)*0.5);
plotMaxAmp = ceil(max(maxAmpRMS)*2);

%% get spectral bounds
nfft = 2^(nextpow2(winlen)+1);
disp('getting spectra')
ltabs = length(tabs);
level05 = NaN(ltabs,1);
level95 = level05;
medFreq = level05;

tmp_events = zpkFilter(indiv_events_orig,0.5,-inf,Fs);
tmp_events = detrend(tmp_events);
[pxx,fxx] = pwelch(tmp_events,[],[],nfft,Fs);
nr = size(pxx,1);
pxx = pxx./repmat(sum(pxx),nr,1);
meanFreq = sum(repmat(fxx,1,ltabs).*pxx)';
pxx = cumsum(pxx)./repmat(sum(pxx),nr,1);

disp('removing duplicates');
for i = 1:ltabs
    disp(i)
    pxx_ = pxx(:,i);
    fxx_ = fxx;
    [pxx_,fxx_] = removeDuplicates(pxx_,fxx_);
    level05(i) = interp1(pxx_,fxx_,0.05);
    level95(i) = interp1(pxx_,fxx_,0.95);
    medFreq(i) = interp1(pxx_,fxx_,0.5);
end
bandwidth = (level95 - level05);
clear tmp_events fxx_ pxx_

%% figure 1, plot SNR
close all;
ttmp = datenum(tabs); % convert to serial number
diff_tabs = diff(ttmp);
diff_tabs = diff_tabs*86400; %convert
tFilt = medfiltSH(ttmp,N);
tFilt = dn2dt(tFilt);
markersize = 15;

% plot snr
figh(1) = figure('units','normalized','outerposition',[0 0 1 1]);
h = plot(tabs,snr,'.','markersize',markersize,'linewidth',1);
xlims = h.Parent.XLim;
ax = gca;
plot_vertical_lines(ax,vlines,10^floor(log10(min(snr))),10^ceil(log10(max(snr)))); %plot vertical lines if requested
ax.XTickLabelRotation = 45;
ax.YScale = 'log';
ylabel('SNR');
title('Signal-to-Noise Ratio')
xlim(xlims);

%% figure 2, plot interevent times
figh(2) = figure('units','normalized','outerposition',[0 0 1 1]);
ha(1) = subplot(211);
hold(ha(1),'on');
h = plot(tabs(1:end-1),diff_tabs,'.','markersize',markersize,'linewidth',1);
xlims = h.Parent.XLim;
plot_vertical_lines(ha(1),vlines,10^floor(log10(min(diff_tabs))),10^ceil(log10(max(diff_tabs)))); %plot vertical lines if requested
ietFilt = medfiltSH(diff_tabs,N);
htmp = plot(tFilt(1:end-1),ietFilt,'linewidth',2);
title('Interevent Time vs. Time')
ylabel('[sec.]')
xlim(xlims);
ha(1).XTickLabelRotation = 45;
ha(1).YScale = 'log';
legend(htmp,[num2str(N),' pt. median filter'],'location','northwest');

% plot cumulative
ha(2) = subplot(212);
cumN = 1:ltabs;
h = plot(tabs,cumN,'linewidth',2);
ylims = h.Parent.YLim;
xlims = h.Parent.XLim;
plot_vertical_lines(ha(2),vlines,plotMinAmp,2*ltabs); %plot vertical lines if requested
title('Cumulative No. of Events');
ylim(ylims);
xlim(xlims);
ha(2).XTickLabelRotation = 45;
linkaxes(ha,'x');

%% figure 3, plot events per minute
figh(3) = figure('units','normalized','outerposition',[0 0 1 1]);
ha(1) = subplot(211);
tw = N/ltabs;
trate = medfiltSH(ttmp(1:end-1),N);
[P,tplot] = pratio(dn2dt(ttmp),N);
%tplot = dn2dt(tplot);
trate = dn2dt(trate);
mint = min([min(tplot) min(trate)]);
maxt = max([max(tplot) max(trate)]);
rate = 60./medfiltSH(diff_tabs,N);
rate = taper(rate,tw);
h = plot(trate,rate,'linewidth',1); hold on;
ylims = h.Parent.YLim;
plot_vertical_lines(ha(1),vlines,plotMinAmp,plotMaxAmp); %plot vertical lines if requested
ha(1).XTickLabelRotation = 25;
title('Events Per Minute')
ylim(ylims);
xlim([mint maxt]);

% plot p ratio
ha(2) = subplot(212);
plot([mint maxt],[1 1],'k--'); hold on;
plot(tplot,P,'linewidth',1); hold on;
plot_vertical_lines(ha(2),vlines,0,ceil(max(P))); %plot vertical lines if requested
title({'\makebox[4in][c]{p-Ratio}','\makebox[4in][c]{mean(iet)/std(iet)}'})
ha(2).XTickLabelRotation = 25;
xlim([mint maxt]);
linkaxes(ha,'x');

%% figure 4, plot mean frequency vs. time
figh(4) = figure('units','normalized','outerposition',[0 0 1 1]);
ha(1) = subplot(211);
h = plot(tabs,medFreq,'.','markersize',markersize); hold on;
ylims = h.Parent.YLim;
xlims = h.Parent.XLim;
medFreqFilt = medfiltSH(medFreq,N);
htmp = plot(tFilt,medFreqFilt,'linewidth',2);
plot_vertical_lines(ha(1),vlines,plotMinAmp,plotMaxAmp); %plot vertical lines if requested
ha(1).XTickLabelRotation = 45;
ylabel('[Hz]');
ylim(ylims);
xlim(xlims);
title('Mean Frequency');
legend(htmp,[num2str(N),' pt. median filter'],'location','best');

% plot bandwidth vs. time
ha(2) = subplot(212);
h = plot(tabs,bandwidth,'.','markersize',markersize,'linewidth',1); hold on;
ylims = h.Parent.YLim;
xlims = h.Parent.XLim;
plot_vertical_lines(ha(2),vlines,plotMinAmp,plotMaxAmp); %plot vertical lines if requested
bwFilt = medfiltSH(bandwidth,N);
htmp = plot(tFilt,bwFilt,'linewidth',2);
ha(2).XTickLabelRotation = 45;
ylabel('[Hz]');
ylim(ylims);
xlim(xlims);
title('Spectral Bandwidth (75\% Level minus 25\% Level)');
legend(htmp,[num2str(N),' pt. median filter'],'location','best');
linkaxes(ha,'x');

%% figure 5, plot maximum amplitude by itself
dcut = cutWindows(maxAmpRMS,N,N-1,false);
rSTD = mean(dcut)./std(dcut); rSTD = 1./rSTD;
t_ = cutWindows(datenum(tabs),N,N-1,false);
median_t_ = median(t_);
median_t_ = dn2dt(median_t_);
mint_ = min([min(median_t_) min(tabs) min(tFilt)]);
maxt_ = max([max(median_t_) max(tabs) max(tFilt)]);
figh(5) = figure('units','normalized','outerposition',[0 0 1 1]);
ha(1) = subplot(211);
h = plot(tabs,maxAmpRMS,'.','markersize',markersize,'linewidth',1); hold on;
ylims = h.Parent.YLim;
plot_vertical_lines(ha(1),vlines,10^floor(log10(min(maxAmpRMS))),10^ceil(log10(max(maxAmpRMS)))); %plot vertical lines if requested
ampFilt = medfiltSH(maxAmpRMS,N);
htmp = plot(tFilt,ampFilt,'linewidth',2);
ha(1).XTickLabelRotation = 45;
ha(1).YScale = 'log';
ylabel('[Counts]');
title('RMS Amplitude');
xlim([mint_ maxt_]);
ylim(ylims);
legend(htmp,[num2str(N),' pt. median filter'],'location','best');

ha(2) = subplot(212);
tSTD = median_t_'; %dn2dt(median_t_);
h = plot(tSTD,rSTD,'linewidth',1); hold on;
ylims = h.Parent.YLim;
plot_vertical_lines(ha(2),vlines,plotMinAmp,plotMaxAmp); %plot vertical lines if requested
ha(2).XTickLabelRotation = 45;
title('RMS CoV');
xlim([mint_ maxt_]);
ylim(ylims);
linkaxes(ha,'x');

% tidy up
clear tplot trate tFilt ampFilt ietFilt P cumN diff_tabs bwFilt level50Filt;

%% save first couple of summary plots
print(figh(1),pFormat,strcat(preFixStr3,'/figure_001_snr'));
print(figh(2),pFormat,strcat(preFixStr3,'/figure_002_iet_cumn'));
print(figh(3),pFormat,strcat(preFixStr3,'/figure_003_rate_pratio'));
print(figh(4),pFormat,strcat(preFixStr3,'/figure_004_freq_bw_time'));
print(figh(5),pFormat,strcat(preFixStr3,'/figure_005_maxAmp'));
indiv_events = normalizeWaveforms(indiv_events);

ccFlag = true;
if ~ccFlag
    maxccp = [];
    maxccn = [];
    plags = [];
    nlags = [];
    flipped_indices = [];
    linear_stack = [];
    family = [];
else
    %% get and plot correlation (similarity) matrix
    disp('get correlation matrix')
    tic;
    [maxccp,plags,maxccn,nlags] = doccFreqCircShift(indiv_events,true);
    %[maxccp,plags,maxccn,nlags] = docc(indiv_events,true);
    toc;
    disp('done getting raw cross correlations');
    %%
    nLarger = maxccn > maxccp;
    %     maxccp(nLarger) = maxccn(nLarger);
    %     plags(nLarger) = nlags(nLarger);
    normcc_down = maxccp;
    
    tic;
    normcc_down = squareform(normcc_down);
    toc;
    disp('done getting correlation matrix');
    
    figh(6) = figure('units','normalized','outerposition',[0 0 1 1]);
    h = imagesc(normcc_down);
    title('Similarity Matrix');
    ylabel('Event Number');
    xlabel('Event Number');
    nSearch = length(normcc_down);
    set(h, 'alphadata', normcc_down >= thresh);
    axis square;
    colorbar;
    caxis([thresh 1]);
    
    %% find repeating events and flipped events
    disp('counting number of repeats')
    nMatches = zeros(nSearch,1);
    n_flipped = 0;
    flipped_indices = cell(nSearch,1);
    nLSquare = squareform(nLarger);
    ccDiff = squareform(maxccn-maxccp);
    diffThresh = 0.15;
    for i = 1:nSearch
        normcc_ = normcc_down(:,i);
        nLSquare_ = nLSquare(:,i);
        nI = normcc_ >= thresh;
        LOCS_ = find(nI);
        nMatches(i) = length(LOCS_);
        
        if ~isempty(LOCS_) % if found events other than yourself, then ...?
            disp(['Event ',num2str(i),'; repeats: ',num2str(nMatches(i)),', ',datestr(tabs(i))])
            LOCS2_ = ccDiff(:,i)>= diffThresh;
            posThresh1Thresh2 = nLSquare_&nI&LOCS2_;
            sumnLSquare_ = sum(posThresh1Thresh2);
            if sumnLSquare_
                n_flipped = n_flipped + sumnLSquare_;
                flipped_indices{i} = find(posThresh1Thresh2);
                disp(['    Number of possible flips: ',num2str(sumnLSquare_)])
                disp(['    Total number of possible flips so far: ',num2str(n_flipped)])
                disp('    Index(es) of possible flips: ');
                disp(flipped_indices{i})
            end
        end
    end
    flipped_indices{nSearch+1} = unique(sort(cat(1,flipped_indices{1:nSearch}))); %store the data
    disp('Display index(es) of all possible flips...');
    disp(flipped_indices{nSearch+1})
    clear normcc_ pol_ nI LOCS_ master posThresh1Thresh2 LOCS2_ ccDiff
    
    %% plot frequency of repeats
    figh(7) = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(211)
    plot(nMatches,'c.','markersize',markersize);
    xlabel('Event Number');
    ylabel('Number of Repetitions');
    
    subplot(212)
    plot(tabs,nMatches,'r.','markersize',markersize);
    ylabel('Number of Repetitions');
    ax = gca;
    ax.XTickLabelRotation = 45;
    
    %% convert similarity matrix to distance vector
    normcc_down = 1-normcc_down; %turn similarity matrix into distance matrix
    dI = normcc_down > 0.99999;
    normcc_down(dI) = 0;
    normcc_down = squareform(normcc_down); %convert distance matrix to distance vector
    clear dI
    
    %% find families (clustering analysis)
    method = 'average';
    [family,l_uniq_indices] = returnClusts(normcc_down,thresh,method);
    ClustTot = length(l_uniq_indices);
    mostMembers = l_uniq_indices(1);
    nFam = find(l_uniq_indices < 2,1) - 1;
    Nsingletons = ClustTot - nFam;
    figh(8) = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(1:nFam,l_uniq_indices(1:nFam),'.','markersize',markersize)
    xlabel('Family Number')
    title('Family Members per Family (min. of 2 members)');
    text(0.8*nFam,0.95*mostMembers,['Number of Singletons: ',num2str(Nsingletons)])
    
    %% save rest of summary plots
    print(figh(6),pFormat,strcat(preFixStr3,'/figure_006_similarity_matrix'));
    print(figh(7),pFormat,strcat(preFixStr3,'/figure_007_frequency_of_repeats'));
    print(figh(8),pFormat,strcat(preFixStr3,'/figure_008_elements_in_family'));
    
    %% print individual families
    l_families = 0;
    maxNFam = 999;
    nFam = min([nFam maxNFam]);
    linear_stack = indiv_events(:,1:nFam);
    
    % convert cc vectors to squareform (2D)
    maxccp = squareform(maxccp);
    maxccn = squareform(maxccn);
    plags = squareform(plags);
    nlags = squareform(nlags);
    
    for i = 1:nFam
        family1 = family{i};
        lfam = length(family1);
        
        indiv_tmp = indiv_events(:,family1);
        
        disp(['Family ',num2str(i),': ',num2str(lfam)]);
        maxccptmp = maxccp(:,family1);
        maxccptmp = maxccptmp(family1,:);
        maxccptmp = squareform(maxccptmp)'; %column vector
        maxccntmp = maxccn(:,family1);
        maxccntmp = maxccntmp(family1,:);
        maxccntmp = squareform(maxccntmp)'; %column vector
        plagstmp = plags(:,family1);
        plagstmp = plagstmp(family1,:);
        plagstmp = squareform(plagstmp)';   %column vector
        nlagstmp = nlags(:,family1);
        nlagstmp = nlagstmp(family1,:);
        nlagstmp = squareform(nlagstmp)';   %column vector
        ccData = [maxccptmp plagstmp maxccntmp nlagstmp];
        disp(['ccData: ',num2str(size(ccData))])
        shifted_data = apply_vdcc(indiv_tmp,ccData);
        linear_stack(:,i) = normalizeWaveforms(mean(shifted_data,2));
        
        if lfam >= 3
            l_families = l_families + 1;
            if i < 10
                fname = strcat(preFixStr0,'/uniq_fams/family_00',num2str(i));
                fname2 = strcat(preFixStr0,'/uniq_fams/time/family_00',num2str(i));
            elseif i < 100
                fname = strcat(preFixStr0,'/uniq_fams/family_0',num2str(i));
                fname2 = strcat(preFixStr0,'/uniq_fams/time/family_0',num2str(i));
            else
                fname = strcat(preFixStr0,'/uniq_fams/family_',num2str(i));
                fname2 = strcat(preFixStr0,'/uniq_fams/time/family_',num2str(i));
            end
            
            fig = figure(1000+i);
            fig.Units = 'normalized';
            fig.OuterPosition = [0 0 1 1];
            
            linear_stack(:,i) = plot_family_and_amps(shifted_data,tabs(family1),maxAmpRMS(family1),1:lfam,5,Fs,-inf,-inf,false,[refTime refTime+length(days)],[plotMinAmp plotMaxAmp]);
            print(pFormat,fname2);
            figure(1000+i);
            print(pFormat,fname2);
            
            figure(10000+i);
            close(10000+i);
            fig = figure(10000+i);
            fig.Units = 'normalized';
            fig.OuterPosition = [0 0 1 1];
            
            ha(1) = subplot(211);
            plot(tabs(family1),maxAmpRMS(family1),'o','color',[0.5 0.5 0.5],'markersize',markersize,'linewidth',2); hold on;
            [P_,S_,MU_] = polyfit(ttmp(family1),log10(maxAmpRMS(family1)),1);
            plot(tabs(family1),10.^(polyval(P_,ttmp(family1),S_,MU_)),'k-','linewidth',2);
            title(['Family ',num2str(i),', N = ',num2str(length(family1))])
            
            plot_vertical_lines(ha(1),vlines,1,plotMaxAmp); %plot vertical lines if requested
            xlim([refTime refTime+length(days)]);
            ylim([plotMinAmp plotMaxAmp]);
            ylabel('Maximum Amplitude (Counts)');
            ha(1).YScale = 'log';
            ha(1).XTickLabelRotation = 45;
            text(refTime+0.1,plotMaxAmp*0.65,['Amp. Ratio (Largest/Smallest): ',...
                num2str(max(maxAmpRMS(family1)) / min(maxAmpRMS(family1)))]);
            
            ha(2) = subplot(212);
            tdum = (0:winlen-1)/Fs;
            plot(tdum,normalizeWaveforms(shifted_data,1,0),'-','color',[0.5 0.5 0.5],'linewidth',1);
            hold(ha(2),'on');
            plot(tdum,normalizeWaveforms(linear_stack(:,i),1,0),'k','linewidth',2);
            xlabel('time [sec.]')
            ylabel('Normalized Amplitude');
            axis(ha(2),'tight');
            
            figure(10000+i);
            pause(1);
            print(pFormat,fname);
            
            close(1000+i);
            close(10000+i);
        end
    end
    
    %% meta analysis
    l_families = min([500 l_families]);
    [maxccp_,plags_,maxccn_,nlags_] = doccFreqDom(linear_stack(:,1:l_families),true);
    %[maxccp_,plags_,maxccn_,nlags_] = docc(linear_stack(:,1:l_families),true);
    ccData_ = [maxccp_,plags_,maxccn_,nlags_];
    shifted_families = apply_vdcc(linear_stack(:,1:l_families),ccData_,false);
    plot_family(shifted_families,1:l_families,8,Fs);
    
    maxccp = squareform(maxccp)';   %column vector
    maxccn = squareform(maxccn)';   %column vector
    plags = squareform(plags)';     %column vector
    nlags = squareform(nlags)';     %column vector
end

function plot_vertical_lines(ax,vlines,plotMinAmp,plotMaxAmp)
hold on;
lvlines = length(vlines);
if lvlines
    for i = 1:lvlines
        vltmp = vlines{i};
        plot(ax,[vltmp vltmp],[plotMinAmp plotMaxAmp],'k--','linewidth',1);
    end
else
    return
end