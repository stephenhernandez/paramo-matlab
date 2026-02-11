function [tnew,newMags,pks,matches,E,pTime,ccnorm,T,stack,templates,S,t,locs_] = findDuplicates(varargin)
%functionDefaults = {"igepn2022gssh.txt",datetime(2022,04,06),128,100,2,32,10,4,16,true,true,true,1,false,30};
functionDefaults = {"igepn2022fzqp.txt",datetime(2022,03,27),64,100,2,32,0.3,2,8,true,false,true,1,false,100};
%functionDefaults = {"igepn2022froh.txt",datetime(2022,03,22:23),20,20,2,24,8,1/5,4/5,true,true,false,1,false,100};
%functionDefaults = {"igepn2016zrom.txt",datetime(2016,12,31),100,20,2,18,8,0.75,12,true,false,true,1,false,100};
%functionDefaults = {"igepn2016cdty.txt",datetime(2016,04,29),25,10,1,24,10,1,8,false,true,false,3,false,100};
%functionDefaults = {"dias2018mkpo.bulletin",datetime(2018,04,15):datetime(2018,09,01),64,20,2,20,12,2,8,false,true,true,3,true,100};
%functionDefaults = {"dias2018mlef.bulletin",datetime(2018,06,26),100,20,2,12,8,2,8,false,true,false,3,true,100};

optsToUse = functionDefaults;
nVarargin = length(varargin);
optsToUse(1:nVarargin) = varargin;
[eventID,loopDays,finalFs,maxNStations,noise,signal,thresh,lfc,hfc,diffFlag,...
    plotFlag,writeFlag,minNStations,diasFlag,maxDist] = deal(optsToUse{:});
warning off signal:findpeaks:largeMinPeakHeight

%%
eventID = char(eventID);
yyyyStr = eventID(6:9);
eventID = string(eventID);
if diasFlag
    E = readSCBulletin(eventID,diasFlag);
    rawDataDir = '~/data/iguana/BROADBAND/';
else
    disp(['~/phaseInformationSC5/',yyyyStr,'/',eventID]);
    E = readSCBulletin(eventID);
    rawDataDir = '~/rawdata';
end

%%
Pphases = E.Pphases;
chans = pull(Pphases,'chan');
for i = 1:length(chans)
    chan_ = char(chans(i));
    chan_(end) = 'Z';
    chans(i) = string(chan_);
    Pphases(i).chan = chans(i);
end

PphasesN = Pphases;
PphasesE = Pphases;
for i = 1:length(chans)
    chan_ = char(chans(i));
    chan_(end) = 'N';
    PphasesN(i).chan = string(chan_);

    %
    chan_ = char(chans(i));
    chan_(end) = 'E';
    PphasesE(i).chan = string(chan_);
end
Pphases = [Pphases; PphasesN; PphasesE];

%for fernandina (FER1) only, delete for anyother type of event!!!
% n = 0;
% lPhases = length(Pphases);
% for i = 1:lPhases
%     n = n + 1;
%     Pphases(lPhases+n) = Pphases(i);
%     tmpChan = char(Pphases(i).chan);
%     Pphases(lPhases+n).chan = string(strcat(tmpChan(1:2),'N'));
%     
%     n = n + 1;
%     Pphases(lPhases+n) = Pphases(i);
%     tmpChan = char(Pphases(i).chan);
%     Pphases(lPhases+n).chan = string(strcat(tmpChan(1:2),'E'));
% end

% for ptgl only, delete for anyother type of event!!!
% kstnms = pull(Pphases,'stnm');
% [lia,locb] = ismember(kstnms,"PTGL");
% 
% if sum(lia)
%     %locb = locb(lia);
%     Pphases = Pphases(lia);
% end
% Pphases.chan = "HHZ";
% Pphases(2) = Pphases(1);
% Pphases(3) = Pphases(1);
% Pphases(2).chan = "HHN";
% Pphases(3).chan = "HHE";
% 
% Pphases(4) = Pphases(1);
% Pphases(5) = Pphases(1);
% Pphases(6) = Pphases(1);
% Pphases(5).chan = "HHN";
% Pphases(6).chan = "HHE";
% 
% Pphases(4).stnm = "HB10";
% Pphases(5).stnm = "HB10";
% Pphases(6).stnm = "HB10";
% 
% Pphases(4).ntwk = "XF";
% Pphases(5).ntwk = "XF";
% Pphases(6).ntwk = "XF";

[yyyy,mm,dd,HH,MM,SS] = datevec(E.t);
origDay = datetime(yyyy,mm,dd);
latitude = E.lat;
longitude = E.lon;
depth = E.depth;
resThresh = 5;
[S,newPhaseStruct] = getEventDayData(Pphases,yyyy,mm,dd,maxNStations,resThresh,maxDist,false,rawDataDir);
kstnms = pull(newPhaseStruct,'stnm');
knetwks = pull(newPhaseStruct,'ntwk');
kcmpnms = pull(newPhaseStruct,'chan');

pTime = pull(newPhaseStruct,'t');

%%
lpTime = length(pTime);
tnew = [];
newMags = [];
pks = [];
matches = [];
if diasFlag
    for i = 1:lpTime
        S(i).knetwk = "9D";
    end
end

%%
if lpTime < minNStations
    disp('Not enough data, sorry');
    return;
end


ttp = min(pTime) - datetime(yyyy,mm,dd,HH,MM,SS); % inferred travel time based on location
if diffFlag
    S = differentiateWaveforms(S);
end
S = resampleWaveforms(S,finalFs);
S = syncWaveforms(S);
S = filterWaveforms(S,lfc,hfc);
S = shiftWaveforms(S,pTime);

minPtime = min(pTime);
Scut = detrendWaveforms(cutWaveforms(S,minPtime-seconds(noise),0,signal));
templates = double(pull(Scut));
maxAmp = max(abs(templates));
maxAmpRMS = rms(templates); 

templates = normalizeWaveforms(templates);
templatesOrig = flipud(templates);
winlen = size(templatesOrig,1);
lS = length(S);
tnew = [];
pks = [];
matches = [];
ampRatio1 = [];
ampRatio2 = [];
nMatch = 0;
minThresh = 0.02; %2;
maxThresh = 0.99;

for j = 1:length(loopDays)
    disp(' ');
    disp(datestr(loopDays(j)));
    nanFlag = true;

    if loopDays(j) == origDay
        T = S;
        goodI = ~isnat(pull(T,'ref'));
        T = T(goodI);
        T = syncWaveforms(T);
    else
        [yyyy_,mm_,dd_] = datevec(loopDays(j));
        tic;
        T = duplicateStructDifferentDay(S,datetime(yyyy_,mm_,dd_));
        goodI = ~isnat(pull(T,'ref'));
        T = T(goodI);
        toc;
        if diffFlag
            T = differentiateWaveforms(T);
        end
        T = resampleWaveforms (T,finalFs);
        T = syncWaveforms(T);
        T = filterWaveforms(T,lfc,hfc);
    end

    badI = isnat(pull(T,'ref'));
    T(badI) = [];
    lT = length(T);

    %%
    if ~lT
        continue;
    end

    %%
    data = double(pull(T));

    %
    templates = templatesOrig(:,goodI);
    nanI = isnan(data);
    disp(['number of nans per station: ',num2str(sum(nanI))]);
    data(nanI) = 0;
    data2 = data.^2;
    t = getTimeVec(T(1));
    fprintf(['nGood: ',num2str(sum(goodI)),'\n']);

    tic;
    box = ones(winlen,sum(goodI));
    norms = fftfilt(box,data2);
    norms = sqrt(abs(norms));
    snorms = norms > 1;
    ccnorm = fftfilt(templates,data);
    ccnorm = snorms.*ccnorm./norms;
    cI = isfinite(ccnorm);
    if nanFlag
        ccnorm(~cI) = NaN;
        cI = 0 == ccnorm;
        ccnorm(cI) = NaN;
        ccnorm(nanI) = NaN; %put nans where they belong...
    else
        ccnorm(~cI) = 0;
    end
    clear cI %T
    toc;

    %stack = nanmean(ccnorm,2); mad_ = nanstd(stack);
    stack = nanmedian(ccnorm,2); medStack = nanmedian(stack); mad_ = nanmedian(abs(stack-medStack)); %mad_ = mad(stack,1);

    if thresh < 1
        newThresh_ = thresh;
    else
        newThresh_ = thresh*mad_;
    end

    fprintf('original threshold is: %f\n',newThresh_);
    newThresh_ = max([minThresh newThresh_]);
    newThresh_ = min([maxThresh newThresh_]);
    fprintf('new threshold is: %f\n',newThresh_);

    [pks_,locs_] = findpeaks(stack,'MinPeakHeight',newThresh_,'MinPeakDistance',winlen);
    locs_ = locs_-winlen+1;
    locsI = locs_ > 0;
    pks_ = pks_(locsI);
    locs_ = locs_(locsI);
    if ~isfinite(pks_)
        fprintf('no events found :( \n\n');
        continue;
    end
    tnew = [tnew; t(locs_)];
    pks = [pks; pks_];
    ll = length(pks_);
    fprintf('number of events found: %d\n',ll);
    dumData = NaN(size(templatesOrig));
    for ii = 1:ll
        newData = data(locs_(ii):locs_(ii)+winlen-1,:);
        dumData(:,goodI) = newData;
        maxAmp_ = max(abs(dumData));
        maxAmpRMS_ = rms(dumData);
        matches = cat(3,matches,dumData);
        ampRatio1 = [ampRatio1; maxAmp_./maxAmp];
        ampRatio2 = [ampRatio2; maxAmpRMS_./maxAmpRMS];
    end
    nMatch = nMatch + ll;

    if nanFlag
        data(nanI) = NaN;
    end
end

%%
logMagRatio1 = log10(nanmedian(ampRatio1,2));
logMagRatio2 = log10(nanmedian(ampRatio2,2));
if E.nMLv
    mType = 'MLv';
    newMags = nanmean(pull(E.MLv,'value'))+logMagRatio2;
elseif E.nMjma
    mType = 'Mjma';
    newMags = nanmean(pull(E.Mjma,'value'))+logMagRatio2;
elseif E.nML
    mType = 'ML';
    newMags = nanmean(pull(E.ML,'value'))+logMagRatio2;
elseif E.nMsBB
    mType = 'MsBB';
    newMags = nanmean(pull(E.MsBB,'value'))+logMagRatio2;
else
    mType = 'Mwp';
    newMags = nanmean(pull(E.Mwp,'value'))+logMagRatio2;
end
newMags = round(newMags*100)/100;
tnew = tnew - ttp - seconds(noise);

%%
if plotFlag
    stack2 = stack(winlen:end);
    ls2 = length(stack2);
    ampFact = 8;

    if nMatch
        for i = 1:lS
            f1 = matches(:,i,:);
            f1 = squeeze(f1);
            f1 = [f1 flipud(templatesOrig(:,i))];
            plot_family(f1,1:nMatch+1,ampFact,finalFs,-inf,-inf,false,false);
            zoom on;
            title(strcat(knetwks(i),'.',kstnms(i),'.',kcmpnms(i)));
        end
    end
    plot_family(flipud(templatesOrig),(1:lS)',ampFact,finalFs);
    zoom on;

    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(tnew,newMags,'.');
    hold on;
    plot(tnew,newMags,'o');
    zoom on;

    figure('units','normalized','outerposition',[0 0 1 1]);
    ha(1) = subplot(211);
    plot(t(1:ls2),stack2);
    hold on;
    plot(t(locs_),pks_,'.','markersize',20);
    plot([t(1) t(end)],[newThresh_ newThresh_],'k--');
    zoom on;

    ha(2) = subplot(212);
    staNumber = 1;
    plot(t,data(:,staNumber)); hold on; plot(t(locs_),pks_,'.','markersize',20);
    zoom on;

    i = staNumber;
    title(strcat(knetwks(i),'.',kstnms(i),'.',kcmpnms(i)));
    linkaxes(ha,'x');
end

%%
if writeFlag
    disp('attempting to write...');
    if nMatch
        eventID_ = char(eventID);
        idBull = regexp(eventID_,'.bulletin','once');
        idTxt = regexp(eventID_,'.txt','once');

        if ~isempty(idTxt)
            eventID_ = eventID_(1:end-4);
        elseif ~isempty(idBull)
            eventID_ = eventID_(1:end-4);
        end
        lon_ = repmat(longitude,nMatch,1);
        lat_ = repmat(latitude,nMatch,1);
        depth_ = repmat(depth,nMatch,1);
        fname = ['~/duplicateEvents/',eventID_,'_ccOutput.txt'];
        formatSpec = '%d %d %d %d %d %5.2f %5.2f %6.5f %f %f %f %s %s';
        [yyyy,mm,dd,HH,MM,SS] = datevec(tnew);
        str = compose(formatSpec,yyyy,mm,dd,HH,MM,SS,newMags,pks,lon_,lat_,depth_,repmat(eventID_,nMatch,1),repmat(mType,nMatch,1));
        str = string(str);
        fileID = fopen(fname,'w');
        fprintf(fileID,'%s\n',str);
        fclose(fileID);

        saveName = ['~/duplicateEvents/',eventID_,'_dataSave.mat'];
        disp(saveName);
        save(saveName,'tnew','newMags','pks','matches','E','kstnms','mType','knetwks','kcmpnms');
    end
end
