clear; close all; clc;
%cd ~/research/now/cotopaxi/;
%load('tremor_binary_classification/CotopaxiTremorCatalog_BinaryClassifier_M5L50');

tremorHome = fullfile("~","masa","old","research","now","cotopaxi");
cd(tremorHome);
load('CotopaxiTremorCatalog_BinaryClassifier_M5L50');

%%
goodI = (yesScore >= 0.96 & nGood == 10) | (yesScore >= 0.98 & nGood >= 9) | ...
    (yesScore >= 0.99 & nGood >= 7) | (yesScore >= 0.999 & nGood >= 5) | (yesScore >= 0.9999);

% close all;
figure(); plot(featureTime(goodI),1:sum(goodI),'.'); zoom on; grid on;
figure(); plot(featureTime(goodI),nGood(goodI),'.'); zoom on; grid on;
rate = t2r(featureTime(goodI),days(1));
figure(); plot(featureTime(goodI),rate/60,'.'); zoom on; grid on
disp([unique(nGood(goodI)) groupcounts(nGood(goodI))])
allAmps(allAmps < 10) = NaN;

figure(); semilogy(featureTime(goodI),allAmps(:,goodI),'.'); zoom on; grid on;
kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];%"BRRN";"VCES";"SRAM";"TOMA"];

%%
t2 = featureTime(goodI) + minutes(1);
gaps = minutes(diff(t2));
minGap = 3;
minDur = 3;
gapI = gaps >= minGap;
episodeStarts = [t2(1)-minutes(1); t2(find(gapI)+1)-minutes(1)];
episodeEnds = [episodeStarts(2:end) - minutes(gaps(gapI))+minutes(1); t2(end)];
quiescence = episodeStarts(2:end) - episodeEnds(1:end-1);
episodeDurations = episodeEnds - episodeStarts;
durI = episodeDurations > minutes(minDur);
eStarts = episodeStarts(durI);
eEnds = episodeEnds(durI);
eDurs = episodeDurations(durI);
reposeTimes = eStarts(2:end) - eEnds(1:end-1);
eStarts = eStarts(1:end-1);
eEnds = eEnds(1:end-1);
eDurs = eDurs(1:end-1);
TremorStartTime = eStarts;
TremorEndTime = eEnds;
TremorDuration = eDurs;
ReposeTime = reposeTimes;
TotalCycleDuration = diff(episodeStarts(durI));
timetable(TremorStartTime,TremorEndTime,TremorDuration,ReposeTime,TotalCycleDuration)
figure(); semilogy(TremorStartTime,minutes(TremorDuration),'.'); grid on; zoom on;

%%
tBC = featureTime(goodI)+minutes(1);
AMPMODEL = load("CotopaxiTremorCatalogReduced_30MAY2023");
tAM = AMPMODEL.t2;
[lia,locb] = ismember(tBC,tAM);
amp1 = AMPMODEL.M2;
amp2 = allAmps(:,goodI);
amp2 = amp2(:,lia);
amp1 = amp1(locb(lia));
tBC = tBC(lia);
tAM = tAM(locb(lia));

%%
close all;
for i = 1:length(kstnms)
    rawAmp = amp2(i,:)';
    tmpI = isfinite(rawAmp) & rawAmp >= 1;
    tmpI = find(tmpI);

    figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(rawAmp(tmpI),amp1(tmpI),'.'); zoom on; grid on;
    title(kstnms(i)); ylabel('Corrected Amp Model Amplitude'); xlabel('Raw');
    xlim([1e0 1e6]); ylim([1e2 1e5]);

    b_ = flipud(robustfit(log10(rawAmp(tmpI)),log10(amp1(tmpI))));
    aSample = [min(log10(rawAmp(tmpI))); log10(rawAmp(tmpI(randsample(length(tmpI),100)))); max(log10(rawAmp(tmpI)))];
    yq = polyval(b_,aSample);
    hold on; ll = loglog(10.^aSample,10.^yq,'.');
    legStr = sprintf('$b_{1}$: %f, $b_{2}$: %f\n',b_(1),b_(2));
    legend(ll,legStr);

    clear binv;
    binv(1) = 1./b_(1);
    binv(2) = -b_(2)./b_(1);

    PredRaw = polyval(binv,log10(amp1));
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(log10(rawAmp),'.'); hold on; plot(PredRaw,'.'); zoom on;
    residual = log10(rawAmp)-PredRaw;
    madScore = (residual-median(residual,'omitnan'))./mad(residual,1);
    outlierFraction = sum(abs(madScore)>= 3)/length(madScore);
    title(sprintf('%s, %f\n',kstnms(i),outlierFraction));
    outlierI = abs(madScore) >= 3;
    tdum = (1:length(rawAmp))'; plot(tdum(outlierI),log10(rawAmp(outlierI)),'o');
end

%%
mdl = fitlm(log10(amp2([1 2 3 4 6 7 8 9 10 11],:)'),log10(amp2(5,:)'),'RobustOpts','on');
BMORPredictions = (mdl.Coefficients.Estimate)' * [ones(1,size(amp2,2)); log10(amp2([1 2 3 4 6 7 8 9 10 11],:))];
BMORPredictions = BMORPredictions';
bmorI = sum(isfinite([amp2(5,:)' 10.^BMORPredictions]),2) == 2;
figure(); semilogy([amp2(5,bmorI)' 10.^BMORPredictions(bmorI)],'.'); zoom on; grid on;
residual = log10(amp2(5,bmorI)') - BMORPredictions(bmorI);
figure(); plot(residual,'.'); zoom on;
madScore = (residual-median(residual,'omitnan'))./mad(residual,1);
outlierI = abs(madScore) >= 3;
figure(24); hold on; tdum = (1:length(residual))'; plot(tdum(outlierI),residual(outlierI),'o');
bmorOutlierFraction = sum(outlierI)/length(outlierI);

%%
mdl = fitlm(log10(amp2(2:11,:)'),log10(amp2(1,:)'),'RobustOpts','on');
BMORPredictions = (mdl.Coefficients.Estimate)' * [ones(1,size(amp2,2)); log10(amp2(2:11,:))];
BMORPredictions = BMORPredictions';
bmorI = sum(isfinite([amp2(1,:)' 10.^BMORPredictions]),2) == 2;
figure(); semilogy([amp2(1,bmorI)' 10.^BMORPredictions(bmorI)],'.'); zoom on; grid on;
residual = log10(amp2(1,bmorI)') - BMORPredictions(bmorI);
figure(); plot(residual,'.'); zoom on;
madScore = (residual-median(residual,'omitnan'))./mad(residual,1);
outlierI = abs(madScore) >= 3;
figure(26); hold on; tdum = (1:length(residual))'; plot(tdum(outlierI),residual(outlierI),'o');
co1vOutlierFraction = sum(outlierI)/length(outlierI);

%%
close all;
clear ax;
ignoreI = true(11,1);
allOutliers = [];
for i = 1:11
    thisI = i;
    ignoreI(thisI) = false;
    mdl = fitlm(log10(amp2(ignoreI,:)'),log10(amp2(~ignoreI,:)'),'RobustOpts','on');
    BMORPredictions = (mdl.Coefficients.Estimate)' * [ones(1,size(amp2,2)); log10(amp2(ignoreI,:))];
    BMORPredictions = BMORPredictions';
    bmorI = sum(isfinite([amp2(~ignoreI,:)' 10.^BMORPredictions]),2) == 2;
    tdum = tAM(bmorI);

    figure('units','normalized','outerposition',[0 0 1 1]);
    ax(1) = subplot(211);
    semilogy(tdum,[amp2(~ignoreI,bmorI)' 10.^BMORPredictions(bmorI)],'.'); zoom on; grid on;
    residual = log10(amp2(~ignoreI,bmorI)') - BMORPredictions(bmorI);
    title(kstnms(~ignoreI));

    ax(2) = subplot(212);
    plot(tdum,residual,'.'); zoom on; hold on;

    madScore = (residual-median(residual,'omitnan'))./mad(residual,1);
    outlierI = abs(madScore) >= 3;
    sum(outlierI)/length(outlierI)
    allOutliers = [allOutliers; find(outlierI)];

    plot(tdum(outlierI),residual(outlierI),'o');
    co1vOutlierFraction = sum(outlierI)/length(outlierI);
    ignoreI(thisI) = true;
    linkaxes(ax,'x');
end

%%
close all;
goodInversionI = true(sum(bmorI),1);
goodInversionI(unique(allOutliers)) = false;
amp3 = amp2(:,bmorI);
t3 = tAM(bmorI);
figure();
semilogy(t3(goodInversionI),amp3(:,goodInversionI)','.'); zoom on; grid on;

t4 = t3(goodInversionI);
amp4 = amp3(:,goodInversionI);

sourceHeightsASL = 5900; %5050; %(5900:-50:1000)'; %898;
%center_lat = -0.6836; center_lon = -78.4361;
center_lon = -78.4361;
center_lat = -0.6836;
[stla,stlo,stel] = metaDataFromStationList(kstnms);
refEllipse = referenceEllipsoid('wgs84');
lHeights = length(sourceHeightsASL);
meanMAD = NaN(lHeights,1);
mbest = NaN(8,lHeights);
rdum = logspace(log10(0.1),log10(200),251)';
madFilter = 0.03;
maxNumberOfAmpsForInversion = 50000;
d = amp4';
t = t4;
plotFlag = true;
mbest = NaN(length(kstnms)+1,lHeights);
lKstnms = length(kstnms);

for jj = 1:lHeights
    fprintf("<strong>%d</strong>\n",sourceHeightsASL(jj));
    d_ = distance(stla,stlo,center_lat,center_lon,refEllipse)*1e-3;
    d_ = sqrt(d_.^2 + ((sourceHeightsASL(jj) - stel)/1000).^2);

    [d_,sI] = sort(d_);
    stla = stla(sI);
    stlo = stlo(sI);
    stel = stel(sI);
    kstnms = kstnms(sI);

    sampleIndex = sort(randsample(length(t),min([maxNumberOfAmpsForInversion length(t)])));
    tsample = t(sampleIndex);
    ampsScrambled = d(sampleIndex,:);

    nStations = size(ampsScrambled,2);
    G_ = full(Gvdcc(nStations));
    G_ = G_(1:end-1,:);

    fixFlag = true;
    if fixFlag
        %n = 1/2; % median mad: 0.0176244, Q = 30, Vs = 2.3;
        n = 1; % median mad: 0.0176993, Q=60, Vs = 2.3;
        G_ = [getDD(d_) G_];
    else
        G_ = [getDD(log10(d_)) getDD(d_) G_];
    end

    tic;
    ampsScrambled2 = ampsScrambled;
    sizeSamples = size(ampsScrambled,1);
    lCombos = 0.5*lKstnms*(lKstnms-1);
    nn = 1;
    nnMax = sizeSamples*lCombos;
    dd = zeros(nnMax,1);
    G = NaN(nnMax,lKstnms+1);
    for i = 1:sizeSamples
        disp(i);
        amps_ = ampsScrambled(i,:); %[CASC_WA(findI(i)) BONI_WA(findI(i)) ANTS_WA(findI(i)) ANTG_WA(findI(i))];
        amps_ = amps_./median(amps_);
        %G = [G; G_];
        G(nn:nn+lCombos-1,:) = G_;
        if fixFlag
            %dd = [dd; -getDD(log10(amps_')) - n*getDD(log10(d_))];
            dd(nn:nn+lCombos-1) = -getDD(log10(amps_')) - n*getDD(log10(d_));
        else
            dd = [dd; -getDD(log10(amps_'))];
        end
        %disp(i)
        ampsScrambled2(i,:) = amps_;
        nn = nn+lCombos;
    end
    toc;

    if fixFlag
        G = [G; 0 ones(1,nStations)];
    else
        G = [G; 0 0 ones(1,nStations)];
    end
    dd = [dd; 0];

    mbest_ = pinv(G)*dd;
    format long g
    disp(mbest_);

    if fixFlag
        att1 = n*log10(rdum) + mbest_(1)*rdum;              %hernandez
    else
        att1 = mbest_(1)*log10(rdum) + mbest_(2)*rdum;      %hernandez
    end
    att2 = 1.11*log10(rdum) + 0.00189*rdum + 0.591;     %uhrhammer
    gamma1 = -interp1(rdum,att1,1);
    gamma17 = -interp1(rdum,att1,17) + 2;
    gamma100 = -interp1(rdum,att1,100) + 3;
    attGamma1 = att1 + gamma1; %gamma17;
    attGamma100 = att1 + gamma100;

    if plotFlag
        figure('units','normalized','outerposition',[0 0 1 1]);
        semilogx(rdum,att1,'.'); zoom on; grid on; hold on;
        semilogx(rdum,att2,'.');
        ylim([0 6]);

        hold on;
        semilogx(rdum,attGamma1,'.','color',[0.5 0.5 0.5]); grid on;
        hold on;
        semilogx(rdum,attGamma100,'.','color',[0.3 0.3 0.3]); grid on;
        legend('new','uhrhammer','fix17','fix100','location','northwest');
    end

    ampsScrambled3 = ampsScrambled;
    McorrOrig = ampsScrambled;
    Mcorr = McorrOrig;
    for i = 1:nStations
        if fixFlag
            Mcorr(:,i) = log10(Mcorr(:,i)) + (n*log10(d_(i)) + mbest_(1)*d_(i) + gamma1 + mbest_(i+1));
        else
            Mcorr(:,i) = log10(Mcorr(:,i)) + (mbest_(1)*log10(d_(i)) + mbest_(2)*d_(i) + gamma1 + mbest_(i+2));
        end
        ampsScrambled3(:,i) = ampsScrambled3(:,i)*(10^mbest_(i+1));
    end
    fprintf('mean mad: %g\n',mean(mad(Mcorr,1,2)));

    M1 = median(Mcorr,2,"omitnan");
    magErr = mad(Mcorr,1,2);
    if plotFlag
        figure();
        mI = isfinite(M1); %M1 >= 0.2 & magErr <= 0.5;
        ax(1) = subplot(211);
        semilogy(tsample(mI),10.^M1(mI),'.'); zoom on; grid on;
        ax(2) = subplot(212);
        plot(tsample(mI),magErr(mI),'.');
        zoom on; grid on;
        linkaxes(ax,'x');

        figure('units','normalized','outerposition',[0 0 1/2 1]);
        ax(1) = subplot(211);
        loglog(d_,ampsScrambled2,'.'); zoom on;
        grid on; ax = gca; ax.LineWidth = 1.5;
        ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
        xlim([1 20]);

        ax(2) = subplot(212);
        ampsScrambled3 = ampsScrambled3./median(ampsScrambled3,2);
        loglog(d_,ampsScrambled3,'.'); zoom on;
        grid on; ax = gca; ax.LineWidth = 1.5;
        ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
        xlim([1 20]);
        linkaxes(ax,'xy');

        figure();
        semilogy(t,d); zoom on; grid on;
    end

    %
    sIndex2 = magErr <= madFilter;
    ampsScrambled = d(sampleIndex(sIndex2),:);

    nStations = size(ampsScrambled,2);
    G_ = full(Gvdcc(nStations));
    G_ = G_(1:end-1,:);

    fixFlag = true;
    if fixFlag
        %n = 1/2; % median mad: 0.0176244, Q = 30, Vs = 2.3;
        n = 1; % median mad: 0.0176993, Q=50, Vs = 2.3;
        G_ = [getDD(d_) G_];
    else
        G_ = [getDD(log10(d_)) getDD(d_) G_];
    end

    ampsScrambled2 = ampsScrambled;
    sizeSamples = size(ampsScrambled2,1);
    lCombos = 0.5*lKstnms*(lKstnms-1);
    nn = 1;
    nnMax = sizeSamples*lCombos;
    dd = zeros(nnMax,1);
    G = NaN(nnMax,lKstnms+1);
    for i = 1:sizeSamples
        disp(i);
        amps_ = ampsScrambled(i,:); %[CASC_WA(findI(i)) BONI_WA(findI(i)) ANTS_WA(findI(i)) ANTG_WA(findI(i))];
        amps_ = amps_./median(amps_);
        G(nn:nn+lCombos-1,:) = G_;
        if fixFlag
            %dd = [dd; -getDD(log10(amps_')) - n*getDD(log10(d_))];
            dd(nn:nn+lCombos-1) = -getDD(log10(amps_')) - n*getDD(log10(d_));
        else
            dd = [dd; -getDD(log10(amps_'))];
        end
        %disp(i)
        ampsScrambled2(i,:) = amps_;
        nn = nn+lCombos;
    end
    toc;

    if fixFlag
        G = [G; 0 ones(1,nStations)];
    else
        G = [G; 0 0 ones(1,nStations)];
    end
    dd = [dd; 0];

    mbest_ = pinv(G)*dd;
    mbest(:,jj) = mbest_;
    format long g
    disp(mbest_);

    if fixFlag
        att1 = n*log10(rdum) + mbest_(1)*rdum;      %hernandez
    else
        att1 = mbest_(1)*log10(rdum) + mbest_(2)*rdum;      %hernandez
    end
    att2 = 1.11*log10(rdum) + 0.00189*rdum + 0.591;     %uhrhammer
    gamma1 = -interp1(rdum,att1,1);
    gamma17 = -interp1(rdum,att1,17) + 2;
    gamma100 = -interp1(rdum,att1,100) + 3;
    attGamma1 = att1 + gamma1;
    attGamma100 = att1 + gamma100;

    if plotFlag
        figure('units','normalized','outerposition',[0 0 1 1]);
        semilogx(rdum,att1,'.'); zoom on; grid on; hold on;
        semilogx(rdum,att2,'.');
        ylim([0 6]);
        hold on;
        semilogx(rdum,attGamma1,'.','color',[0.5 0.5 0.5]); grid on;
        hold on;
        semilogx(rdum,attGamma100,'.','color',[0.3 0.3 0.3]); grid on;
        legend('new','uhrhammer','fix17','fix100','location','northwest');
    end

    ampsScrambled3 = ampsScrambled;
    McorrOrig = ampsScrambled;
    Mcorr = McorrOrig;
    for i = 1:nStations
        if fixFlag
            Mcorr(:,i) = log10(Mcorr(:,i)) + (n*log10(d_(i)) + mbest_(1)*d_(i) + gamma1 + mbest_(i+1));
        else
            Mcorr(:,i) = log10(Mcorr(:,i)) + (mbest_(1)*log10(d_(i)) + mbest_(2)*d_(i) + gamma1 + mbest_(i+2));
        end
        ampsScrambled3(:,i) = ampsScrambled3(:,i)*(10^mbest_(i+1));
    end
    meanMAD(jj) = mean(mad(Mcorr,1,2));
    fprintf('mean mad: %g\n',meanMAD(jj));

    %
    tsample2 = tsample(sIndex2);
    M1 = median(Mcorr,2,"omitnan");
    magErr = mad(Mcorr,1,2);
    if plotFlag
        figure();
        mI = isfinite(M1); %M1 >= 0.2 & magErr <= 0.5;
        ax(1) = subplot(211);
        semilogy(tsample2(mI),10.^M1(mI),'.'); zoom on; grid on;
        ax(2) = subplot(212);
        plot(tsample2(mI),magErr(mI),'.');
        zoom on; grid on;
        linkaxes(ax,'x');

        %
        figure('units','normalized','outerposition',[0 0 1/2 1]);
        ax(1) = subplot(211);
        loglog(d_,ampsScrambled2,'.'); zoom on;
        grid on; ax = gca; ax.LineWidth = 1.5;
        ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
        xlim([1 20]);

        ax(2) = subplot(212);
        ampsScrambled3 = ampsScrambled3./median(ampsScrambled3,2);
        loglog(d_,ampsScrambled3,'.'); zoom on;
        grid on; ax = gca; ax.LineWidth = 1.5;
        ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
        xlim([1 20]);
        linkaxes(ax,'xy');
    end

    %
    allEffects = 10.^(n*log10(d_) + mbest_(1).*d_ + gamma1 + mbest_(2:end));
    table(kstnms,stla,stlo,stel,d_,10.^(mbest_(2:end)),1./(10.^(mbest_(2:end))),allEffects,...
        'VariableNames',...
        {'SNCL';'Lat.';'Lon.';'Elev.';'Distance';'Correction';'Amplification';'All Effects'})
    toc;
end

amp5 = amp4;
for i = 1:lKstnms
    amp5(i,:) = allEffects(i)*amp5(i,:);
end

figure();
semilogy(t4,median(amp5',2,"omitnan"),'.'); zoom on; grid on;

figure();
semilogy(t4,mad(log10(amp5'),1,2),'.'); zoom on; grid on;

C = loadWaveforms(datetime(2023,02,14),1,kstnms,["BHZ";"HHZ";"SHZ"],"EC");
Cf = scaleWaveforms(transferWaveforms((cutWaveforms(C,dateshift(C(1).ref,'start','day')+hours(0),0,hours(24))),3,6,4,20,"vel",1,false),1e9);
for i = 1:length(Cf)
    kstnm_ = Cf(i).kstnm;
    [lia,locb] = ismember(kstnm_,kstnms);
    if lia
        Cf(i) = scaleWaveforms(Cf(i),allEffects(locb));
    end
end
plotWaveforms(Cf);
