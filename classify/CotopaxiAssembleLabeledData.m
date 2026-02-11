clear; close all; clc;

cd ~/research/now/cotopaxi/
load("CotopaxiTremorCatalogReduced_30MAY2023")
kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];
kcmpnms = ["HHZ";"BHZ";"BHZ";"BHZ";"BHZ";"HHZ";"SHZ";"SHZ";"BHZ";"HHZ";"SHZ"];
lKstnms = length(kstnms);

%% too many stations
% kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB";"BRRN";"VCES";"SRAM";"TOMA"];
% kcmpnms = ["HHZ";"BHZ";"BHZ";"BHZ";"BHZ";"HHZ";"SHZ";"SHZ";"BHZ";"HHZ";"SHZ";"BHZ";"HHZ";"BHZ";"BHZ"];
%%

% we need CO1V, so filter beyond December 16
t2 = t2 - minutes(1);
tI = t2 >= datetime(2022,12,16);
M2 = M2(tI);
t2 = t2(tI);

tw = 200;
lfc = 1/5;
hfc = 10;
newFs = 25;
nF = 101;
newFxx = logspace(log10(lfc),log10(hfc),nF)';
nfft2 = 2^11;
PXX = zeros(nF,lKstnms);
npoles = 4;

featureCountTrue = 0;
featureCountFalse = 0;
maxFeatures  = 3e4;
nFeatures = 0.5*lKstnms*(lKstnms-1)*nF;

featuresTrue = NaN(maxFeatures,nFeatures);
featuresFalse = featuresTrue;

timeTrue = NaT(maxFeatures,1);
timeFalse = timeTrue;

uniqDays = unique(dateshift(t2,'start','day'));
lUniqDays = length(uniqDays);
for i = 1:lUniqDays
    day_ = uniqDays(i);
    dI = t2 >= day_ & t2 < day_ + 1;
    tNow = t2(dI);

    n = 0;
    S = populateWaveforms(lKstnms);
    for j = 1:lKstnms
        kstnm_ = kstnms(j);
        kcmpnm_ = kcmpnms(j);
        S_ = loadWaveforms(day_,1,kstnm_,kcmpnm_);
        if isnat(S_.ref)
            continue;
        end
        n = n+1;
        S_ = detrendWaveforms(S_);
        S_ = taperWaveforms(S_,tw);
        S_ = scaleWaveforms(transferWaveforms(S_,lfc/2,-inf,npoles,100,"vel",1,false),1e9);
        S(j) = S_;
    end

    if n ~= lKstnms
        fprintf("not enough stations for day: %s\n",day_);
        continue;
    end

    Sf = syncWaveforms(S,0,1,0);
    n = length(Sf);
    if n ~= lKstnms
        fprintf("synching error: %s\n",day_);
        continue;
    end

%     Sf = detrendWaveforms(...
%         intWaveforms(...
%         taperWaveforms(...
%         filterWaveforms(...
%         taperWaveforms(...
%         syncWaveforms(...
%         detrendWaveforms(...
%         differentiateWaveforms(Sf))),tw),lfc,hfc),tw)));
    Sf = interpolateWaveforms(Sf);
    Sf = resampleWaveforms(Sf,newFs);

    allGaps = [];
    for j = 1:lKstnms
        gapInfo = Sf(j).gapInfo;
        allGaps = [allGaps; gapInfo];
    end

    if size(allGaps,1) > 0
        allGaps = sortrows(allGaps);
        for j = 1:lKstnms
            Sf(j).gapInfo = allGaps;
            Sf(j).gapFlag = true;
        end
        Sf = resampleWaveforms(Sf,newFs*2);
        Sf = resampleWaveforms(Sf,newFs);
    end
    Sf = detrendWaveforms(Sf);

    %%
    newStart = dateshift(min(pull(Sf,'ref')),'start','minute');
    newEnd = dateshift(max(pull(Sf,'ref')+pull(Sf,'e')),'end','minute');
    tMinutes = (newStart:minutes(1):newEnd-minutes(1))';

    Sf = cutWaveforms(Sf,newStart,0,newEnd-newStart,false,true);
    Sf = nanGapWaveforms(Sf,0);

    dcut = [];
    for j = 1:lKstnms
        d_ = Sf(j).d;
        dcut_ = cutWindows(d_,newFs*60,0,false); %do not detrend
        nanI = ~isfinite(dcut_);
        dcut_(nanI) = 0;
        dcut = cat(3,dcut,dcut_);
    end

    badMinutes = max(squeeze(sum(~isfinite(dcut) | dcut == 0)),[],2) == newFs*60;
    tMinutesGood = tMinutes(~badMinutes);

    lMinutes = length(tMinutesGood);
    if ~lMinutes
        fprintf("no valid minutes for day: %s\n",day_);
        continue;
    end

    goodIndex = ismember(tMinutesGood,t2);
    trueLabelsI = find(goodIndex);
    dcut = dcut(:,~badMinutes,:);

    thisDayTrueCount = 0;
    thisDayFalseCount = thisDayTrueCount;
    for j = 1:lMinutes
        slice = squeeze(dcut(:,j,:));
        [slice,fxxN] = pmtm(detrend(slice),4,nfft2,newFs); %do not normalize
        for k = 1:lKstnms
            slice_ = slice(:,k);
            pxxInterp = interp1(fxxN,slice_,newFxx);
            PXX(:,k) = pxxInterp;
        end

        % get spectral ratios here
        theseFeatures = NaN(nFeatures,1);
        tmpFeatureCount = 1;
        for k = 1:lKstnms-1
            pxx_ = PXX(:,k);
            for kk = k+1:lKstnms
                pxx2_ = PXX(:,kk);
                sr = pxx_./pxx2_;
                lsr = length(sr);
                theseFeatures(tmpFeatureCount:tmpFeatureCount+lsr-1) = sr;
                tmpFeatureCount = tmpFeatureCount + lsr;
            end
        end

        %
        if ismember(j,trueLabelsI)
            featureCountTrue = featureCountTrue + 1;
            %
            if featureCountTrue > maxFeatures
                featureCountTrue = featureCountTrue - 1;
                fprintf("cannot add anymore true labels\n");
                continue;
            end
            featuresTrue(featureCountTrue,:) = theseFeatures;
            timeTrue(featureCountTrue) = tMinutesGood(j);
            thisDayTrueCount = thisDayTrueCount + 1;
        else
            featureCountFalse = featureCountFalse + 1;
            if featureCountFalse > maxFeatures
                featureCountFalse = featureCountFalse - 1;
                fprintf("cannot add anymore false labels\n");
                continue;
            end
            featuresFalse(featureCountFalse,:) = theseFeatures;
            timeFalse(featureCountFalse) = tMinutesGood(j);
            thisDayFalseCount = thisDayFalseCount + 1;
        end
    end
    fprintf("%s: %d true labels added, <strong>%d</strong> cumulative\n",day_,thisDayTrueCount,featureCountTrue);
    fprintf("%s: %d false labels added, %d cumulative\n",day_,thisDayFalseCount,featureCountFalse);
end

%%
saveFlag = true;
if saveFlag
    cd ~/research/now/cotopaxi/;
    clear S*
    clear p*
    clear PXX
    clear t2
    clear d_ dc*
    clear dI
    clear M2
    clear tI
    clear nanI
    save("CotopaxiFeatures_v1");
end
