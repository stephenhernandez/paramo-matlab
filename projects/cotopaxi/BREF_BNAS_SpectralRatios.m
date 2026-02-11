tStart = datetime(2012,01,01);
tEnd = datetime(2023,02,09);
dayVec = (tStart:tEnd)';
lDays = length(dayVec);
finalFs = 50;
secDur = 60;

n = 0;
spectralRatio = NaN(2049,lDays);
dayKeep = false(lDays,1);
totalDailyWindows = NaN(lDays,1);
for i = 1:lDays
    day_ = dayVec(i);
    fprintf("processing: %s\n",day_);

    C = loadWaveforms(day_,1,["BREF";"BNAS"],["BHZ";"HHZ"],"EC");
    refs = pull(C,'ref');
    badI = isnat(refs);
    goodI = ~badI;
    if ~sum(goodI)
        fprintf("no data for: %s, continuing...\n",day_);
        continue;
    end

    % get acceleration traces
    Cf = detrendWaveforms(...
        scaleWaveforms(...
        transferWaveforms(...
        detrendWaveforms(differentiateWaveforms(C)),1/100,-inf,4,finalFs,"vel",true,false),1e9));
    Cf = nanGapWaveforms(Cf,0);
    Cf = syncWaveforms(Cf,false,true,false);

    traceDur = Cf(1).e;
    if traceDur < minutes(10)
        fprintf("trace too short for: %s, continuing...\n",day_);
        continue;
    end

    lC = length(Cf);
    if lC < 2
        fprintf("only one trace, continuing...\n");
        continue;
    end

    Crsam = rsamWaveforms(Cf,60,-inf,-inf,4,true,false);
    tRSAM = getTimeVec(Crsam) - minutes(1);
    dRSAM = pull(Crsam);
    dRSAM(dRSAM < 1) = NaN;
    dRSAM = prod(dRSAM,2);

    ampThresh = prctile(dRSAM,25); % 25th percentile; %median(dRSAM,1,"omitnan");
    mI = dRSAM <= ampThresh & isfinite(dRSAM);
    if ~sum(mI)
        fprintf("no valid points to extract spectra from, continuing...\n");
        continue;
    end

    tcut = tRSAM(mI);
    mPxx_ = [];
    d1 = pull(cutWaveforms(Cf(1),tcut,0,minutes(1)));
    d2 = pull(cutWaveforms(Cf(2),tcut,0,minutes(1)));

    pxx1 = pmtm(detrend(d1),[],[],finalFs);
    [pxx2,fxx] = pmtm(detrend(d2),[],[],finalFs);

    n = n + 1;
    spectralRatio(:,n) = sqrt(abs(median(pxx2./pxx1,2,"omitnan")));
    dayKeep(i) = true;
    totalDailyWindows(i) = sum(mI);

    cd ~/research/now/cotopaxi/
    save("SpectralRatios_BNASoverBREF_v2.mat",'spectralRatio','fxx','dayKeep','dayVec','n');
end
spectralRatio = spectralRatio(:,1:n);

cd ~/research/now/cotopaxi/
save("SpectralRatios_BNASoverBREF_v2.mat",'spectralRatio','fxx','dayKeep','dayVec','n');
