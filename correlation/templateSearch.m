function [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,...
    madMain,nUsedMain,ccnorm,tLong,polMain,normers] = ...
    templateSearch(S,T,writeFlag,plotFlag,verboseFlag,thresh,...
    customPrefix,fileVersionNumber,CANUSEGPU)
if nargin < 3
    writeFlag = true;
end

if nargin < 4
    plotFlag = false;
end

if nargin < 5
    verboseFlag = true;
end

if nargin < 6
    thresh = 8;
end

if nargin < 7
    customPrefixFlag = false;
    customPrefix = [];
else
    customPrefixFlag = true;
end

if nargin < 8
    fileVersionNumber = 1;
end

if nargin < 9
    CANUSEGPU = false;
end
warning off signal:findpeaks:largeMinPeakHeight

%%
maxN = 1e6;
ccMain = NaN(maxN,1);
madMain = ccMain;
ampMain = ccMain;
magMain = ccMain;
nUsedMain = ccMain;
dMag = ccMain;
polMain = ccMain;
tMain = NaT(maxN,1);
evidMain = repmat("-",maxN,1);
templateNumber = ccMain;

%%
lS = length(S);
Snpts = pull(S,"npts");
if any(Snpts ~= Snpts(1))
    fprintf("error: day data trace have different lengths\n");
    return;
end
Snpts = Snpts(1);

[nChanT,nTemplateT] = size(T);
if lS ~= nChanT
    fprintf("error: number of channels do not match\n");
    return;
end

%%
winlen = pull(T,"npts");
if any(winlen ~= winlen(1))
    fprintf("error: some templates have different lengths\n");
    return;
end

winlen = winlen(1);
if winlen > Snpts
    fprintf("day data is too short\n");
    return;
end

%% processing of long traces!
tLong = getTimeVec(S);
traceStartTime = tLong(1);
mpd = 0.2*winlen;
box = ones(winlen,nChanT);
origData = pull(S);
lData = size(origData,1);
if CANUSEGPU
    origData = gpuArray(origData);
end

normers = fftfilt(box,origData.^2); %the sum of squares %convn(origData.^2,box);
normers = sqrt(abs(normers)); %this only gets done once
maxLength = lData - winlen + 1;
badNormersI = normers < 1 | ~isfinite(normers);
normers(badNormersI) = NaN;
fprintf("number of bad normers: %d\n",sum(badNormersI,"all"));

%% processing of individual multiplets
istart = 1;
ntot = 0; % sum of all detected events over all possible multiplets
day_ = dateshift(traceStartTime,"start","day");
fprintf("processing %d trace(s) that start(s) on: <strong>%s</strong>\n",...
    nChanT,traceStartTime);

%%
for i = 1:nTemplateT %nTemplateT is number of columns of T struct
    T_ = T(:,i); %
    T_ = detrendWaveforms(T_);
    dT = pull(T_);
    ampWeights = abs(dT).^1;
    ampWeights = ampWeights./sum(ampWeights);
    referenceAmp = sqrt(sum(ampWeights.*abs(dT.^2)));

    dT = flipud(dT);
    dT = dT./rssq(dT);          % normalize my templates
    if CANUSEGPU
        dT = gpuArray(dT);
    end

    cc = fftfilt(dT,origData);  % cross-correlation
    ccnormOrig = cc./normers;
    badNormersI2 = badNormersI | ~isfinite(ccnormOrig);
    missingSNCLs = sum(badNormersI2,2);
    uniqMissingSNCLs = unique(missingSNCLs)';
    weights = ones(nChanT+1,1);

    ccnormOrig = mean(ccnormOrig,2,"omitnan");
    ccnorm = ccnormOrig;
    for ii = uniqMissingSNCLs
        wi = missingSNCLs == ii;
        thisMAD = mad(ccnorm(wi),1);
        we = 1./thisMAD;
        if isfinite(we)
            weights(ii+1) = we;
            ccnorm(wi) = we.*ccnorm(wi);
        end
    end

    if thresh < 1 % CC threshold has been given
        newThresh = thresh;
        ccnorm = abs(ccnormOrig);   %clobbers the sign
    else
        % a MAD threshold has been given
        mad_ = mad(ccnorm,1);
        newThresh = thresh*mad_;
        ccnorm = abs(ccnorm);       %clobbers the sign
    end
    [cc_,locs_] = findpeaks(ccnorm,"MINPEAKDISTANCE",mpd,"MINPEAKHEIGHT",newThresh);
    pols_ = sign(ccnormOrig(locs_));
    locs_ = locs_ - winlen + 1;
    locsI = locs_ <= maxLength & locs_ > 0 & isfinite(cc_);
    cc_ = cc_(locsI);
    locs_ = locs_(locsI);
    pols_ = pols_(locsI);

    n_ = sum(locsI);
    if ~n_
        fprintf("no events found on day: %s\n",day_);
        continue;
    end

    locsOrig = locs_ + winlen - 1;
    t_ = tLong(locs_);
    ntot = ntot + n_;
    thisMag = T_(1).eqmag;
    if ~isfinite(thisMag)
        thisMag = 0;
    end
    thisID = T_(1).evid;

    myMissingSNCLs = missingSNCLs(locsOrig);
    if thresh < 1
        ccMain(istart:istart+n_-1) = cc_;
    else
        ccMain(istart:istart+n_-1) = cc_./weights(myMissingSNCLs+1);
    end

    tMain(istart:istart+n_-1) = t_;
    templateNumber(istart:istart+n_-1) = i;
    nUsedMain(istart:istart+n_-1) = lS - myMissingSNCLs;
    evidMain(istart:istart+n_-1) = thisID;
    if thresh < 1
        madMain(istart:istart+n_-1) = cc_.*weights(myMissingSNCLs+1);
    else
        madMain(istart:istart+n_-1) = cc_;
    end
    polMain(istart:istart+n_-1) = pols_;

    for j = 1:n_
        locs__ = locs_(j);
        if locs__ <= maxLength
            indivEvents = origData(locs__:locs__+winlen-1,:);
            goodStreamsI = ~badNormersI(locs__+winlen-1,:);
        else
            indivEvents = origData(locs__:end,:);
            goodStreamsI = ~badNormersI(end,:);
        end

        %%
        indivEvents = detrend(indivEvents);
        amp_ = sqrt(sum(ampWeights.*abs(indivEvents.^2)));
        amp_ = amp_(goodStreamsI);
        refAmps = referenceAmp(goodStreamsI);

        indiv_station_dMags = log10(amp_./refAmps);
        amp_ = median(amp_,"omitnan");
        ampMain(istart+j-1) = amp_;

        dMag_ = median(indiv_station_dMags,"omitnan");
        dMag(istart+j-1) = dMag_;
        magMain(istart+j-1) = thisMag+dMag_;
    end

    if verboseFlag
        fprintf("template number %d:, %d detection(s) on: %s\n",i,n_,traceStartTime);
    end
    istart = istart+n_; %new future start index
end

ccMain = ccMain(1:ntot);
ampMain = ampMain(1:ntot);
tMain = tMain(1:ntot);
templateNumber = templateNumber(1:ntot);
magMain = magMain(1:ntot);
evidMain = evidMain(1:ntot);
dMag = dMag(1:ntot);
madMain = madMain(1:ntot);
nUsedMain = nUsedMain(1:ntot);
polMain = polMain(1:ntot);

%%
[tMain,sI] = sort(tMain);
dMag = dMag(sI);
ampMain = ampMain(sI);
ccMain = ccMain(sI);
evidMain = evidMain(sI);
templateNumber = templateNumber(sI);
magMain = magMain(sI);
madMain = madMain(sI);
nUsedMain = nUsedMain(sI);
polMain = polMain(sI);

%%
if plotFlag
    tile_spacing = "compact";
    line_width = 1;
    if newThresh < 1
        linkedPlot(tLong,[ccnormOrig origData(:,1)],"k-",tile_spacing,line_width);
    else
        linkedPlot(tLong,[ccnorm origData(:,1)],"k-",tile_spacing,line_width);
    end
end

%%
if ~writeFlag
    if verboseFlag
        fprintf("not writing...\n");
    end
    return;
end

%%
if ~ntot
    if verboseFlag
        fprintf("no repeats found on day: <strong>%s</strong>, cannot write...\n",traceStartTime);
    end
    return;
end

fprintf("attempting to write...\n");
currentDay = dateshift(S(1).ref,"start","day");
currentDay.Format = "yyyy.MM.dd";
kstnm = S(1).kstnm;

if customPrefixFlag
    prefix = fullfile("~","masa","template_search",customPrefix);
else
    prefix = fullfile("~","masa","template_search",kstnm);
end

if ~exist(prefix,"dir")
    mkdir(prefix);
end
fName = fullfile(prefix,strcat("dayTemplateSearch_",string(currentDay),".txt"));
fName2 = fullfile(prefix,strcat("dayTemplateSearch_",string(currentDay),".mat"));

[yyyy,mm,dd,HH,MM,SS] = datevec(tMain);
if fileVersionNumber == 1
    formatSpec = '%04d %02d %02d %02d %02d %05.2f %+7.5f %010.1f %07.5f %04d %07.5f %04.2f %02d %s';
    str = compose(formatSpec,yyyy,mm,dd,HH,MM,SS,dMag,ampMain,ccMain,...
        templateNumber,madMain,magMain,nUsedMain,evidMain);
elseif fileVersionNumber == 2
    formatSpec = '%04d %02d %02d %02d %02d %05.2f %+7.5f %010.1f %07.5f %04d %07.5f %04.2f %02d %s %d';
    str = compose(formatSpec,yyyy,mm,dd,HH,MM,SS,dMag,ampMain,ccMain,...
        templateNumber,madMain,magMain,nUsedMain,evidMain,polMain);
end
str = string(str);

%%
try
    fileID = fopen(fName,"w");
    fprintf(fileID,"%s\n",str);
    fclose(fileID);
    save(fName2,'tMain','dMag','ampMain','ccMain','madMain',...
        'templateNumber','magMain','nUsedMain','evidMain','polMain');
catch
    fprintf("couldnt write file for day: %s, may need to reprocess...!\n",...
        traceStartTime);
    return;
end