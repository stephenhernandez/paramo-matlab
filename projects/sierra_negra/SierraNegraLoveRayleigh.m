clear; close all;

%%
refEllipse = referenceEllipsoid('wgs84');
[stla,stlo] = metaDataFromStationList("VCH1");
load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag');
globalt = t;
globaleqlat = eqlat;
globaleqlon = eqlon;
globaleqdepth = eqdepth;
globaleqmag = eqmag;
clear t eq*

%%
velModel = 'iasp91';
phaseList = 'p,P,Pdiff';
dists = distance(stla,stlo,globaleqlat,globaleqlon,refEllipse)*1e-3;

%%
%Amplitude Constants
VEL = 3.5;                                  % surface wave velocity in km/sec
Elr0 = 1e-8;                                % minimum strain threshold
X0 = log10(Elr0*20*(VEL*1e9)/(2*pi)) + 2;   % strain threshold converted to log_10(amplitude)

%Spatial Constants
km2deg = 1/111.19;                        % coarse conversion factor
Dmin = 200;                               % distance threshold in km

Idist = (dists >= Dmin);
Ithresh = (globaleqmag > X0+1.66*log10(dists*km2deg)); %the +2 is hidden in X0
Itime = globalt >= datetime(2013,07,01) & globalt <= datetime(2019,07,01);
goodI = (Idist & Ithresh & Itime);

amps = 10.^(globaleqmag - 1.66*log10(dists.*km2deg) - 2);
strain = amps/(20*(VEL*1e9)/(2*pi));
SWAT = globalt + seconds(dists./6.5);     % large amplitude arrival time
TphaseTT = seconds(dists./2);         % duration object
fastLove = seconds(dists./4.3);

%%
goodMag = globaleqmag(goodI);
goodDepth = globaleqdepth(goodI);
goodDist = dists(goodI);
goodStrain = strain(goodI);
goodt = globalt(goodI);
goodSWAT = SWAT(goodI);
goodTpTT = TphaseTT(goodI);
goodLove = fastLove(goodI);

[goodStrain,sI] = sort(goodStrain,'descend');
goodt = goodt(sI);
goodSWAT = goodSWAT(sI);
goodTpTT = goodTpTT(sI);
goodDepth = goodDepth(sI);
goodDist = goodDist(sI);
goodMag = goodMag(sI);
goodLove = goodLove(sI);

%%
travTimes = getTravelTimes(velModel,phaseList,goodDepth,goodDist,true);
sI = find(~isfinite(travTimes));
if sI
    goodt(sI) = [];
    goodSWAT(sI) = [];
    goodTpTT(sI) = [];
    goodDepth(sI) = [];
    goodDist(sI) = [];
    goodMag(sI) = [];
    goodStrain(sI) = [];
    goodLove(sI) = [];
    travTimes(sI) = [];
end

%%
maxMinutesAfter = 4000;
maxTriggers = min([length(goodt) 1e2]);

%%
lfc = 1/120;
hfc = 1/15;
newFs = 1;
winlen = 100;

%%
BAZ = [];
pR = [];
pL = [];
t = [];
ergR = [];
ergL = [];
dur = [];
scaledER = [];
scaledEL = [];

%%
for i = 1:100%:maxTriggers
    currentMag = goodMag(i);
    %TA = minutes(maxMinutesAfter);
    %tStart = goodt(i) + seconds(travTimes(i)) - seconds(60);
    tStart = goodt(i) + goodLove(i) - seconds(60);
    tEnd = goodt(i) + seconds(travTimes(i)) + goodTpTT(i);
    TA = tEnd - tStart;
    
    S = extractWaveforms(tStart,TA,"VCH1","HHZ","EC");  %was totalDuration in prev. iteration
    S(2) = extractWaveforms(tStart,TA,"VCH1","HHN","EC");  %was totalDuration in prev. iteration
    S(3) = extractWaveforms(tStart,TA,"VCH1","HHE","EC");  %was totalDuration in prev. iteration
    
    if any(isnat(pull(S,'ref')))
        S = extractWaveforms(tStart,TA,"VCH1","BHZ","EC"); %was totalDuration in prev. iteration
        S(2) = extractWaveforms(tStart,TA,"VCH1","BHN","EC"); %was totalDuration in prev. iteration
        S(3) = extractWaveforms(tStart,TA,"VCH1","BHE","EC"); %was totalDuration in prev. iteration
        if any(isnat(pull(S,'ref')))
            fprintf('skipping event: %s,%f\n',datestr(tStart),currentMag);
            continue;
        end
    end
    
    %%
    %S = interpolateWaveforms(detrendWaveforms(S));
    S = resampleWaveforms(detrendWaveforms(S),newFs); %now at 50 Hz.
    %S = filterWaveforms(taperWaveforms(S),lfc,hfc);
    npts = S.npts;
    if npts < winlen
        warning('data not long enough, skipping');
        continue;
    end
    
    tt = getTimeVec(S(1));
    lt = length(tt);
    
    [estBAZ_,corrcoefs,angs,corrCoeffEstBAZ,loveTT,rayTT,peakRadial_,peakTransverse_] = rayleighBAZsearch(S,lfc,hfc,0);
    BAZ = [BAZ; estBAZ_];
    
    pR = [pR; peak2peak(peakRadial_)];
    pL = [pL; 0.5*peak2peak(peakTransverse_)];
    ergR = [ergR; sum(peakRadial_.^2)];
    ergL = [ergL; sum(peakTransverse_.^2)];
    t = [t; goodt(i)];
    dur = [dur; seconds(TA)];
    
    lenRS = min([lt length(peakRadial_) length(peakTransverse_)]);
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(tt(1:lenRS),peakRadial_,'linewidth',2);
    hold on;
    plot(tt(1:lenRS),abs(peakTransverse_),'linewidth',2);
    title(['pR:',num2str(pR(end)),'; pL: ',num2str(pL(end)),...
        '; dur: ',num2str(dur(end)),...
        '; energyR: ',num2str(ergR(end)),...
        '; energyL: ',num2str(ergL(end)),...
        '; scaled energy R: ',num2str(ergR(end)./dur(end)),...
        '; scaled energy L: ',num2str(ergL(end)./dur(end))]);
    zoom on;
    pause(2);
end

%%
scaledER = ergR./dur;
scaledEL = ergL./dur;



