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
TphaseTT = seconds(dists./1.5);         % duration object

%%
goodMag = globaleqmag(goodI);
goodDepth = globaleqdepth(goodI);
goodDist = dists(goodI);
goodStrain = strain(goodI);
goodt = globalt(goodI);
goodSWAT = SWAT(goodI);
goodTpTT = TphaseTT(goodI);

[goodStrain,sI] = sort(goodStrain,'descend');
goodt = goodt(sI);
goodSWAT = goodSWAT(sI);
goodTpTT = goodTpTT(sI);
goodDepth = goodDepth(sI);
goodDist = goodDist(sI);
goodMag = goodMag(sI);

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
    travTimes(sI) = [];
end

%%
maxMinutesAfter = 400;
maxTriggers = min([length(goodt) 1e2]);
beta = NaN(maxTriggers,1);
ampRatio = beta;

%%
secDur = 10;
lfc = 6;
hfc = 24;
Fs = 100;
winlen = secDur*Fs;

%%
newFs = 50;
minThresh = 0.35;
minDistance = newFs*secDur/2;

load vch1_2018_aligned
d2 = normalizeWaveforms(resample(detrend(d2(1:winlen,:)),1,2)); %now at 50Hz.
winlen = size(d2,1);

%%
masterCat = [];
taAll = [];
tStarts = [];

%%
nn = 0;
for i = 1:maxTriggers
    for iii = 1:2
        tic;
        currentMag = goodMag(i);
        TA = minutes(maxMinutesAfter);
        tStart = goodt(i) + seconds(travTimes(i));
        if iii == 1
            tStart = tStart - minutes(maxMinutesAfter);
        end
        tEnd = tStart + minutes(maxMinutesAfter);
        %tEnd = min([tStart + TA goodt(i) + goodTpTT(i)]);
        TA = tEnd - tStart; %duration object
        %         totalDuration = minutes(maxMinutesAfter);
        %         if TA < maxMinutesAfter
        %             tStart = tStart - (totalDuration - TA);
        %         end
        
        if i > 1
            nearness = tStart - goodt(1:i-1);
            if  any(nearness > 0 & nearness <= days(1/2))
                %skipping because this potential trigger is too close in time
                %to a previous trigger
                %disp(['yikes, skipping: ',datestr(goodSWAT(i))])
                continue;
            end
        end
        
        S = extractWaveforms(tStart,TA,"VCH1","HHZ","EC");  %was totalDuration in prev. iteration
        if isnat(S(1).ref)
            S = extractWaveforms(tStart,TA,"VCH1","BHZ","EC"); %was totalDuration in prev. iteration
            if isnat(S(1).ref)
                fprintf('skipping event: %s,%f\n',datestr(tStart),currentMag);
                continue;
            end
        end
        S = interpolateWaveforms(detrendWaveforms(S));
        S = resampleWaveforms(detrendWaveforms(S),newFs); %now at 50 Hz.
        S = filterWaveforms(taperWaveforms(S),lfc,hfc);
        npts = S.npts;
        if npts < winlen
            warning('data not long enough, skipping');
            continue;
        end
        
        t = getTimeVec(S);
        d = S.d;
        
        if sum(~isfinite(d))
            continue;
        end
        
        %%
        %fprintf('running cross correlation, this can take a while\n');
        cc = fftfilt(flipud(d2),d);
        
        %%
        %fprintf('getting normalization constants\n');
        norms = fftfilt(ones(winlen,1),d.^2);
        norms = sqrt(abs(norms));
        cc = cc./norms;
        
        %%
        %fprintf('finished running cross correlation, now looping through results to find peaks\n');
        Locs = cell(size(d2,2),1);
        Pks = cell(size(d2,2),1);
        I = Pks;
        for j = 1:size(d2,2)
            cc_ = cc(:,j);
            [pks,locs] = findpeaks(cc_,'MINPEAKHEIGHT',minThresh,'MINPEAKDISTANCE',minDistance);
            if sum(locs >= winlen)
                Pks{j} = pks(locs>=winlen);
                Locs{j} = locs(locs>=winlen) - winlen + 1;
                I{j} = repmat(j,sum(locs>=winlen),1);
            end
        end
        Locs = cat(1,Locs{:});
        Pks = cat(1,Pks{:});
        I = cat(1,I{:});
        
        if isempty(Locs)
            continue;
        end
        
        tmatch = t(Locs);
        ampmatch = norms(Locs+winlen-1);
        [tmatch,matchSort] = sort(tmatch);
        ampmatch = ampmatch(matchSort);
        Pks = Pks(matchSort);
        I = I(matchSort);
        [tt,ccc,removeIndices,~] = removeRepeatedMatches(tmatch,Pks,secDur/2);
        %%
        if isempty(tt)
            continue;
        end
        for kk = 1:length(removeIndices)
            ri = removeIndices{kk};
            I(ri) = [];
            ampmatch(ri) = [];
        end
        
        taAll = [taAll; TA];
        tStarts = [tStarts; datenum(tStart)];
        masterCat = [masterCat; datenum(tt) ccc I ampmatch];
        nn = nn + 1;
        fprintf('n: %u; time: %s; mag.: %f; number potential events: %u; log strain: %f; distance: %f; depth: %f; winlen: %f\n',...
            nn,datestr(goodt(i)),currentMag,length(tt),log10(goodStrain(i)),goodDist(i),goodDepth(i),minutes(TA));
        toc;
        fprintf('\n');
    end
end

%%
save('~/research/now/sierra_negra/correlationDetector_v2','masterCat',...
    'taAll','tStarts','goodDist','goodt','goodDepth','goodStrain','goodMag');
