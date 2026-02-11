clear; close all;

cd ~/research/now/cotopaxi/
load("CotopaxiTremorCatalogReduced_30MAY2023")
kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];
kcmpnms = ["HHZ";"BHZ";"BHZ";"BHZ";"BHZ";"HHZ";"SHZ";"SHZ";"BHZ";"HHZ";"SHZ"];
lKstnms = length(kstnms);
minStations = 3;

cd ~/research/now/cotopaxi/tremor_binary_classification
%load('tremor_binary_classification/CotopaxiGentleBoostEnsembleCompact_v4'); experimentalFlag = true;
% load('tremor_binary_classification/CotopaxiGentleBoostEnsembleCompact_10Splits25Learners_FS1_min10min20'); experimentalFlag = true;
load('~/research/now/cotopaxi/tremor_binary_classification/test04'); experimentalFlag = true;

%
close all;
tw = 500;
lfc = 1/5;
hfc = 10;
newFs = 25;
nF = 101;
newFxx = logspace(log10(lfc),log10(hfc),nF)';
nfft2 = 2^11;
PXX = zeros(nF,lKstnms);
npoles = 4;
Hd = zpkOperator(3,6,newFs,npoles); %for plotting purposes

dailyMaxSamples  = 1500;
nFeatures = 0.5*lKstnms*(lKstnms-1)*nF;
sampleCount_ = 1;
maxSamples  = 1440*300;
allFeatures_ = NaN(maxSamples,sum(goodFeaturesI));
allAmps_ = NaN(lKstnms,maxSamples);
featureTime_ = NaT(maxSamples,1);
yesScore_ = NaN(maxSamples,1);
nGood_ = yesScore_;

tStart = datetime(2020,12,22);
tEnd = datetime(2020,12,22);
dayInc = 1;
plotFlag = true;
saveFlag = false; %fNamePrefix = 'CotopaxiTremorCatalog_BinaryClassifier_v04';

dayVec = (tStart:tEnd)';
lDays = length(dayVec);

%
for i = 1:lDays
    fprintf("\n");
    dailySampleCount = 0;
    allFeaturesTmp = NaN(dailyMaxSamples,nFeatures);
    allAmpsTmp = NaN(lKstnms,dailyMaxSamples);
    featureTimeTmp = NaT(dailyMaxSamples,1);

    day_ = dayVec(i);
    dI = t2 >= day_ & t2 < day_ + 1;
    tNow = t2(dI);

    n = 0;
    S = populateWaveforms(lKstnms);
    for j = 1:lKstnms
        kstnm_ = kstnms(j);
        kcmpnm_ = kcmpnms(j);
        S_ = loadWaveforms(day_,dayInc,kstnm_,kcmpnm_);
        if isnat(S_.ref)
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

    nOrig = n;
    if ~experimentalFlag
        Sf = detrendWaveforms(S);
        Sf = resampleWaveforms(Sf,newFs);
        Sf = taperWaveforms(Sf,tw);
        Sf = filterWaveforms(Sf,lfc,hfc);
    else
        Sf = resampleWaveforms(S,newFs);
        Sf = detrendWaveforms(Sf);
    end

    n = length(Sf);
    if n ~= lKstnms
        fprintf("synching error: %s\n",day_);
        continue;
    end

    %
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

    goodIndex = ismember(tMinutesGood,t2);
    dcut = dcut(:,~badMinutes,:);

    thisDayTrueCount = 0;
    thisDayFalseCount = thisDayTrueCount;
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

    %
    tic;
    allFeaturesOrig = allFeaturesTmp;
    allFeaturesTmp = allFeaturesTmp(1:dailySampleCount,goodFeaturesI);
    allAmpsTmp = allAmpsTmp(:,1:dailySampleCount);

    [Yfit,score] = predict(cGentleBoostEnsemble2,allFeaturesTmp);
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
        datestr(day_),dailySampleCount,sum(yesScoreTmp>=0.99),sampleCount_);
end

%%
allAmps_ = allAmps_(:,1:sampleCount_-1);
allFeatures_ = allFeatures_(1:sampleCount_-1,:);
featureTime_ = featureTime_(1:sampleCount_-1);
nGood_ = nGood_(1:sampleCount_-1);
yesScore_ = yesScore_(1:sampleCount_-1);

if saveFlag
    clearvars -except allFeatures allAmps featureTime yesScore nGood sampleCount
    save('hello');
end

%%
if plotFlag
    close all;
    figure();
    ax(1) = subplot(211);
    plot(featureTimeTmp,Yfit,'.'); zoom on; grid on;
    ax(2) = subplot(212);
    plot(featureTimeTmp,score,'.'); zoom on; grid on;
    linkaxes(ax,'x');

    figure('units','normalized','outerposition',[0 0 0.6 1]);
    plot(featureTimeTmp(Yfit),1:sum(Yfit),'.'); zoom on; hold on;
    plot(featureTimeTmp(score(:,2)>=0.9),1:sum(score(:,2)>=0.9),'.'); zoom on;
    plot(featureTimeTmp(score(:,2)>=0.99),1:sum(score(:,2)>=0.99),'.'); zoom on;
    legend('$\ge 0.5$','$\ge 0.9$','$\ge 0.99$','Location','Best');

    figure();
    plot(featureTimeTmp,score(:,2),'.'); zoom on; grid on;

    %     figure();
    %     semilogy((1:sum(goodFeaturesI))',allFeaturesTmp(930:931,:)','.'); zoom on; grid on;
    %     xlim([1 sum(goodFeaturesI)]);

    AM = load('~/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v28JUL2023.mat','allEffects');
    allEffects = AM.allEffects;

    figure();
    spax(1) = subplot(311);
    semilogy(featureTimeTmp,allAmpsTmp,'.'); zoom on; grid on;
    spax(2) = subplot(312);
    allAmpsTmp2 = allAmpsTmp;
    for i = 1:length(allEffects)
        allAmpsTmp2(i,:) = allEffects(i).*allAmpsTmp(i,:);
    end
    semilogy(featureTimeTmp,allAmpsTmp2,'.'); zoom on; grid on;
    spax(3) = subplot(313);
    semilogy(featureTimeTmp,mad(log10(allAmpsTmp2),1,1),'.'); zoom on; grid on;
    linkaxes(spax,'x');

    figure();
    plot(tMinutesGood,nGoodTmp,'.'); zoom on; grid on;
end
