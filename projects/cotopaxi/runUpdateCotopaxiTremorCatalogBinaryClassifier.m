function [allAmps_,allFeatures_,featureTime_,nGood_,yesScore_] = ...
    runUpdateCotopaxiTremorCatalogBinaryClassifier(tStart,tEnd,cGentleBoostEnsemble2,goodFeaturesI,experimentalFlag)

kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];
kcmpnms = ["HHZ";"HHZ";"BHZ";"BHZ";"HHZ";"BHZ";"HHZ";"HHZ";"BHZ";"HHZ";"HHZ"];
%kcmpnms = ["HHZ";"BHZ";"BHZ";"BHZ";"BHZ";"BHZ";"SHZ";"SHZ";"BHZ";"HHZ";"SHZ"];
lKstnms = length(kstnms);
minStations = 3;

tw = 500;
lfc = 1/5;
hfc = 10;
newFs = 25;
nF = 101;
newFxx = logspace(log10(lfc),log10(hfc),nF)';
nfft2 = 2^11;
PXX = zeros(nF,lKstnms);
npoles = 4;
Hd = zpkOperator(3,6,newFs,npoles);

dailyMaxSamples  = 1500;
nFeatures = 0.5*lKstnms*(lKstnms-1)*nF;
sampleCount_ = 1;
maxSamples  = 1440*300;
allFeatures_ = NaN(maxSamples,sum(goodFeaturesI));
allAmps_ = NaN(lKstnms,maxSamples);
featureTime_ = NaT(maxSamples,1);
yesScore_ = NaN(maxSamples,1);
nGood_ = yesScore_;

dayInc = 1;
dayVec = (tStart:tEnd)';
lDays = length(dayVec);

%%
for i = 1:lDays
    dailySampleCount = 0;
    allFeaturesTmp = NaN(dailyMaxSamples,nFeatures);
    allAmpsTmp = NaN(lKstnms,dailyMaxSamples);
    featureTimeTmp = NaT(dailyMaxSamples,1);

    day_ = dayVec(i);

    n = 0;
    S = populateWaveforms(lKstnms);
    for j = 1:lKstnms
        kstnm_ = kstnms(j);
        %kcmpnm_ = kcmpnms(j);
        %S_ = loadWaveforms(day_,dayInc,kstnm_,kcmpnm_);
        S_ = loadWaveforms(day_,dayInc,kstnm_,["HHZ";"BHZ";"SHZ"]);
        S_ = S_(1);
        if isnat(S_.ref)
            fprintf("no data for: %s\n",kstnm_);
            continue;
        end

        if strcmp(kstnm_,"BREF") && day_ >= datetime(2025,01,23)
            fprintf("no data for: %s\n",kstnm_);
            continue;
        end

        n = n+1;
        if experimentalFlag
            S_ = detrendWaveforms(S_);
            S_ = taperWaveforms(S_,tw);
            S_ = scaleWaveforms(transferWaveforms(S_,lfc/2,-inf,npoles,100,"vel",1,false),1e9);
            S(j) = S_;
        else
            S(j) = S_;
        end
    end

    if n < minStations
        fprintf("not enough stations for day: %s\n",day_);
        continue;
    end

    if ~experimentalFlag
        Sf = detrendWaveforms(S);
        Sf = resampleWaveforms(Sf,newFs);
        Sf = taperWaveforms(Sf,tw);
        Sf = filterWaveforms(Sf,lfc,hfc);
    else
        Sf = syncWaveforms(S,false,true,true);
        Sf = resampleWaveforms(Sf,newFs);
        Sf = nanGapWaveforms(Sf,0);
        Sf = padWaveforms(Sf);
    end

    n = length(Sf);
    if n ~= lKstnms
        fprintf("synching error: %s\n",day_);
        continue;
    end
    fprintf("read at least %d stations for day: %s\n",n,day_);

    %%
    newStart = dateshift(min(pull(Sf,'ref')),'start','minute');
    newEnd = dateshift(max(pull(Sf,'ref')+pull(Sf,'e')),'start','minute');
    tMinutes = (newStart:minutes(1):newEnd-minutes(1))';

    Sf = cutWaveforms(Sf,newStart,0,newEnd-newStart,false,true);
    Sf = nanGapWaveforms(Sf,0);

    dcut = [];
    badI = find(isnat(pull(Sf,'ref')));
    npts = max(pull(Sf,'npts'));
    for j = 1:lKstnms
        if ismember(j,badI)
            d_ = zeros(npts,1);
        else
            d_ = Sf(j).d;
        end
        dcut_ = cutWindows(d_,newFs*60,0,false); %do not detrend
        nanI = ~isfinite(dcut_);
        dcut_(nanI) = 0;
        dcut = cat(3,dcut,dcut_);
    end

    nGoodTmp = sum(~(squeeze(sum(~isfinite(dcut) | dcut == 0))== newFs*60),2);
    badMinutes = nGoodTmp < minStations;
    tMinutesGood = tMinutes(~badMinutes);
    nGoodTmp(badMinutes) = [];

    lMinutes = length(tMinutesGood);
    if ~lMinutes
        fprintf("no valid minutes for day: %s\n",day_);
        continue;
    end

    dcut = dcut(:,~badMinutes,:);
    for j = 1:lMinutes
        slice = squeeze(dcut(:,j,:));
        tdSliceOrig = slice;
        sliceFiltered = filter(Hd,tdSliceOrig);
        Amp3Hz6Hz_ = rms(sliceFiltered)';

        [slice,fxxN] = pmtm(detrend(slice),4,nfft2,newFs); %do not normalize
        for k = 1:lKstnms
            slice_ = slice(:,k);
            pxxInterp = interp1(fxxN,slice_,newFxx); %<-- the magic
            PXX(:,k) = pxxInterp;
        end

        % get spectral ratios here
        theseFeatures = NaN(nFeatures,1);
        tmpFeatureCount = 1;
        for k = 1:lKstnms-1
            pxx_ = PXX(:,k);
            for kk = k+1:lKstnms
                pxx2_ = PXX(:,kk);
                sr = pxx_./pxx2_; %spectral ratios
                lsr = length(sr);
                theseFeatures(tmpFeatureCount:tmpFeatureCount+lsr-1) = sr;
                zeroI = theseFeatures < 1e-10 | ...
                    theseFeatures > 1e10 | ~isfinite(theseFeatures);
                theseFeatures(zeroI) = NaN;
                tmpFeatureCount = tmpFeatureCount + lsr;
            end
        end

        dailySampleCount = dailySampleCount + 1;
        allFeaturesTmp(dailySampleCount,:) = theseFeatures;
        allAmpsTmp(:,dailySampleCount) = Amp3Hz6Hz_;
        featureTimeTmp(dailySampleCount) = tMinutesGood(j);
    end

    %%
    tic;
    allFeaturesTmp = allFeaturesTmp(1:dailySampleCount,goodFeaturesI);
    allAmpsTmp = allAmpsTmp(:,1:dailySampleCount);

    [~,score] = predict(cGentleBoostEnsemble2,allFeaturesTmp);
    featureTimeTmp = featureTimeTmp(1:dailySampleCount);
    toc;

    yesScoreTmp = score(:,2);
    allFeatures_(sampleCount_:sampleCount_+dailySampleCount-1,:) = allFeaturesTmp;
    allAmps_(:,sampleCount_:sampleCount_+dailySampleCount-1) = allAmpsTmp;
    featureTime_(sampleCount_:sampleCount_+dailySampleCount-1) = featureTimeTmp;
    yesScore_(sampleCount_:sampleCount_+dailySampleCount-1) = yesScoreTmp;
    nGood_(sampleCount_:sampleCount_+dailySampleCount-1) = nGoodTmp;

    sampleCount_ = sampleCount_ + dailySampleCount;
    fprintf("day: %s, today sample count: %d, number yes: %d, cumulative samples: %d\n",...
        day_,dailySampleCount,sum(yesScoreTmp>=0.85),sampleCount_);
    save(fullfile("~","masa","old","research","now","cotopaxi",...
        sprintf("CotopaxiTremorDetectorTmp_%s",datestr(day_,"yyyymmdd"))),...
        "allFeaturesTmp","allAmpsTmp","featureTimeTmp","yesScoreTmp","nGoodTmp");
end

%%
allAmps_ = allAmps_(:,1:sampleCount_-1);
allFeatures_ = allFeatures_(1:sampleCount_-1,:);
featureTime_ = featureTime_(1:sampleCount_-1);
nGood_ = nGood_(1:sampleCount_-1);
yesScore_ = yesScore_(1:sampleCount_-1);
end