clear;
close all;

lfc = 3;
hfc = 6;
f = sqrt(lfc*hfc);
Q = 60;
Vs = 2.3; %km/s
B = (pi*f)/(Q*Vs);
alpha = 1/2;
tw = 0.0002;
newFs = 50; %8*hfc;
npoles = 4;
waFlag = false;
deconvolveFlag = true;

plotFlag = ~true;
madFilter = 0.02;
sourceHeightsASL = 5900; %5050; %(5900:-50:1000)'; %898;
%center_lat = -0.6836; center_lon = -78.4361;
center_lon = -78.4361;
center_lat = -0.6836;
kstnmsOrig = ["CO1V";"BREF";"BVC2";"BTAM";...
    "BMOR";"BNAS";"SLOR"];

C = loadWaveforms(datetime(2023,02,14),1,kstnmsOrig,...
    ["BHZ";"HHZ"],"EC");

%
if plotFlag
    close all;
end

secDur = 60;
% Ccut = detrendWaveforms(cutWaveforms(C,dateshift(C(1).ref,'start','day')+...
%     hours(7)+minutes(30)+seconds(00),0,...
%     hours(5)));
Ccut = detrendWaveforms(cutWaveforms(C,dateshift(C(1).ref,'start','day')+...
    hours(9)+minutes(30)+seconds(00),0,...
    hours(2)));

Cf = detrendWaveforms(...
    scaleWaveforms(...
    transferWaveforms(...
    taperWaveforms(...
    detrendWaveforms(medfiltWaveforms(differentiateWaveforms(Ccut),3)),tw),...
    lfc,hfc,npoles,newFs,"disp",deconvolveFlag,waFlag),1e9));
Cf = nanGapWaveforms(Cf,0);

newFs2 = 1;
technique1 = true;
if technique1
    Cfenv2 = convWaveforms(powWaveforms(Cf,-1),newFs*secDur+1);
    Cfenv2 = resampleWaveforms(Cfenv2,newFs2);
    Cfenv2 = medfiltWaveforms(Cfenv2,newFs2*secDur+1,true);
else
    Cfenv2 = envelopeWaveforms(Cf);
    Cfenv2 = resampleWaveforms(Cfenv2,newFs2);
    Cfenv2 = medfiltWaveforms(Cfenv2,9,true);
    Cfenv2 = filterWaveforms(Cfenv2,-inf,1/newFs2/secDur,1,[],true);
end
Cfenv2 = syncWaveforms(Cfenv2,true,true,true);

kstnms = pull(Cfenv2,'kstnm');
[stla,stlo,stel] = metaDataFromStationList(kstnms);
refEllipse = referenceEllipsoid('wgs84');

lHeights = length(sourceHeightsASL);
meanMAD = NaN(lHeights,1);
mbest = NaN(8,lHeights);
rdum = logspace(log10(0.1),log10(200),251)';
for jj = 1:lHeights
    fprintf("<strong>%d</strong>\n",sourceHeightsASL(jj));
    d_ = distance(stla,stlo,center_lat,center_lon,refEllipse)*1e-3;
    d_ = sqrt(d_.^2 + ((sourceHeightsASL(jj) - stel)/1000).^2);

    [d_,sI] = sort(d_);
    stla = stla(sI);
    stlo = stlo(sI);
    stel = stel(sI);
    kstnms = kstnms(sI);
    Cfenv2 = Cfenv2(sI);

    d = pull(Cfenv2);
    t = getTimeVec(Cfenv2);

    sampleIndex = sort(randsample(length(t),min([5000 length(t)])));
    winlen = secDur*newFs2;
    sI2 = sampleIndex > winlen;
    sampleIndex = sampleIndex(sI2);
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

    G = [];
    dd = [];
    ampsScrambled2 = ampsScrambled;
    for i = 1:size(ampsScrambled,1)
        amps_ = ampsScrambled(i,:); %[CASC_WA(findI(i)) BONI_WA(findI(i)) ANTS_WA(findI(i)) ANTG_WA(findI(i))];
        amps_ = amps_./median(amps_);
        G = [G; G_];
        if fixFlag
            dd = [dd; -getDD(log10(amps_')) - n*getDD(log10(d_))];
        else
            dd = [dd; -getDD(log10(amps_'))];
        end
        %disp(i)
        ampsScrambled2(i,:) = amps_;
    end

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

    G = [];
    dd = [];
    ampsScrambled2 = ampsScrambled;
    for i = 1:size(ampsScrambled,1)
        amps_ = ampsScrambled(i,:); %[CASC_WA(findI(i)) BONI_WA(findI(i)) ANTS_WA(findI(i)) ANTG_WA(findI(i))];
        amps_ = amps_./median(amps_);
        G = [G; G_];
        if fixFlag
            dd = [dd; -getDD(log10(amps_')) - n*getDD(log10(d_))];
        else
            dd = [dd; -getDD(log10(amps_'))];
        end
        %disp(i)
        ampsScrambled2(i,:) = amps_;
    end
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
end

%%
processDifferentDay = false;
if processDifferentDay
    tic;
    C = loadWaveforms(datetime(2023,02,14),1,kstnmsOrig,...
        ["BHZ";"HHZ";"SHZ"],"EC");

    Cf = detrendWaveforms(...
        scaleWaveforms(...
        transferWaveforms(...
        taperWaveforms(...
        detrendWaveforms(medfiltWaveforms(differentiateWaveforms(C),3)),tw),...
        lfc,hfc,npoles,newFs,"disp",deconvolveFlag,waFlag),1e9));
    Cf = nanGapWaveforms(Cf,0);
    if technique1
        Cfenv2 = envelopeWaveforms(Cf);
        Cfenv2 = resampleWaveforms(Cfenv2,newFs2);
        Cfenv2 = medfiltWaveforms(Cfenv2,newFs2*secDur,false);
    else
        Cfenv2 = envelopeWaveforms(Cf);
        Cfenv2 = resampleWaveforms(Cfenv2,newFs2);
        Cfenv2 = medfiltWaveforms(Cfenv2,9,true);
        Cfenv2 = filterWaveforms(Cfenv2,-inf,1/newFs2/secDur,1,[],true);
    end
    Cfenv2 = syncWaveforms(Cfenv2,true,true,true);
    Cf = nanGapWaveforms(Cf,0);
    
    winlen = secDur*newFs2;
    d2 = pull(Cfenv2);
    t2 = getTimeVec(Cfenv2);
    ref = t2(1);
    tref = dateshift(ref,'end','minute');
    iStart = t2i(tref,ref,1/newFs2);

    %%
    d2 = d2(iStart:winlen:end,:);
    npts = size(d2,1);
    t2 = tref + seconds(secDur*(0:npts-1)');

    dI = d2 <= 1;
    d2(dI) = 1;
    McorrOrig2 = d2;
    Mcorr2 = McorrOrig2;
    for i = 1:nStations
        if fixFlag
            Mcorr2(:,i) = log10(Mcorr2(:,i)) + (n*log10(d_(i)) + mbest_(1)*d_(i) + gamma1 + mbest_(i+1));
        else
            Mcorr2(:,i) = log10(Mcorr2(:,i)) + (mbest_(1)*log10(d_(i)) + mbest_(2)*d_(i) + gamma1 + mbest_(i+2));
        end
    end
    fprintf('mean mad: %g\n',mean(mad(Mcorr,1,2)));

    %
    M2 = median(Mcorr2,2,"omitnan");
    magErr2 = mad(Mcorr2,1,2);
    mI2 = magErr2 <= 0.05 & isfinite(magErr2) & isfinite(M2);
    figure('units','normalized','outerposition',[0 0 1/2 1]);
    M2 = median(Mcorr2,2,"omitnan");
    magErr2 = mad(Mcorr2,1,2);
    clear ax;
    ax(1) = subplot(311);
    semilogy(t2(mI2),10.^M2(mI2),'.'); zoom on; grid on; hold on;
    semilogy(t2(~mI2),10.^M2(~mI2),'.');
    ax(2) = subplot(312);
    plot(t2(mI2),magErr2(mI2),'.');
    zoom on; grid on; hold on;
    plot(t2(~mI2),magErr2(~mI2),'.');
    ax(3) = subplot(313);
    plot(getTimeVec(Cf),Cf(1).d);
    zoom on; grid on; hold on;
    linkaxes(ax,'x');

    figure('units','normalized','outerposition',[0 0 1/2 1]);
    plot(t2(mI2),1:sum(mI2),'.'); zoom on; grid on;

    figure('units','normalized','outerposition',[0 0 1/2 1]);
    plot(t2(mI2),t2r(t2(mI2),hours(1)),'.'); zoom on; grid on;
    toc;

    knetwks = pull(Cfenv2,'knetwk');
    kstnms = pull(Cfenv2,'kstnm');
    kholes = pull(Cfenv2,'khole');
    kcmpnms = pull(Cfenv2,'kcmpnm');
    allMySNCLs = strcat(knetwks,kstnms,kholes,kcmpnms);

    table(allMySNCLs,stla,stlo,stel,d_,10.^(mbest_(2:end)),1./(10.^(mbest_(2:end))),...
        'VariableNames',...
        {'SNCL';'Lat.';'Lon.';'Elev.';'Distance';'Correction';'Amplification'})
end
