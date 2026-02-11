% function [M,magErr,t,d,Mcorr,Lats,Lons,Depths] = ...
%     loopCotopaxiTremorLocation3D(tStart,tEnd)
% if nargin < 1
%     tStart = datetime(2023,01,13);
% end
%
% if nargin < 2
%     tEnd = dateshift(datetime('now'),'start','day');
% end
%
% if nargin < 2
%     plotFlag = false;
% end

clear; close all;
%compare 13 and 30 of january, i have nice long time-lapses from both days.
% tStart = datetime(2022,02,17);
% tEnd = datetime(2022,02,17);

%tStart = datetime(2023,01,26);
%tEnd = datetime(2023,01,26);
% tStart = datetime(2024,03,21);
% tEnd = datetime(2024,03,21);
%tStart = datetime(2023,05,12);
%tEnd = datetime(2023,05,12);
tStart = datetime(2024,04,03);
tEnd = datetime(2024,04,03);

% tStart = datetime(2015,08,14);
% tEnd = datetime(2015,08,14);
dayInc = 1;
plotFlag = ~false;
saveFlag = ~true;
center_lon = -78.4361;
center_lat = -0.6836;

%
%dataHome = '~/research/now/cotopaxi/';
dataHome = '~/igdata';
Model = load(fullfile(dataHome,'CotopaxiTremorAttenuationModel_v15022023.mat'));
%Model = load('~/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v30012023.mat');
%Model = load('~/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v11022023.mat');

dOrig = Model.d_;
mbest = Model.mbest;
allMySNCLs= Model.allMySNCLs;
deconvolveFlag = Model.deconvolveFlag;
waFlag = Model.waFlag;
npoles = Model.npoles;
lfc = Model.lfc;
hfc = Model.hfc;
kstnmsOrig = Model.kstnmsOrig;
secDur = Model.secDur;
n = Model.n;
newFs = Model.newFs;
fixFlag = true; %Model.fixFlag;
gamma1 = Model.gamma1;
newFs2 = Model.newFs2;
charSNCL = char(allMySNCLs);
allMySNCLs = string(charSNCL(:,1:6));

stationCorrectionsOrig = mbest(2:end);
%load("~/research/now/cotopaxi/Cotopaxi3DCoarseTremorModel");
load(fullfile(dataHome,"Cotopaxi3DCoarseTremorModel_2Steps"));
lX = length(X);
distsCoarseOrig = distsCoarse;
stlaOrig = stla;
stloOrig = stlo;
stelOrig = stel;

%%
M = [];
magErr = M;
t = M;
d = M;
correctedAmps = M;
Lats = M;
Lons = M;
Depths = M;

dayVec = (tStart:tEnd)';
lDays = length(dayVec);
refEllipse = referenceEllipsoid('wgs84');
for i = 1:lDays
    day_ = dayVec(i);
    C = loadWaveforms(day_,dayInc,kstnmsOrig,...
        ["BHZ";"HHZ"],"EC");

    Cf = detrendWaveforms(...
        scaleWaveforms(...
        transferWaveforms(...
        detrendWaveforms(...
        C),lfc,hfc,npoles,newFs,"vel",deconvolveFlag,waFlag),1e9));

    Cf = nanGapWaveforms(Cf,0);
    Cf = syncWaveforms(Cf,0,0,true);
    Cf = nanGapWaveforms(Cf,0);

    Tsncls = allMySNCLs;
    Csncls = strcat(pull(Cf,'knetwk'),pull(Cf,'kstnm'));
    [lia,locb] = ismember(Csncls,Tsncls);

    if sum(lia) < 3
        continue;
    end

    distsCoarse = distsCoarseOrig(:,locb(lia));
    stla = stlaOrig(locb(lia));
    stlo = stloOrig(locb(lia));
    stel = stelOrig(locb(lia));

    Cf = Cf(lia);
    nStations = length(Cf);
    stationCorrections = stationCorrectionsOrig(locb(lia));

    Cfenv2 = envelopeWaveforms(Cf);
    Cfenv2 = resampleWaveforms(Cfenv2,newFs2);
    Cfenv2 = medfiltWaveforms(Cfenv2,newFs2*secDur,false);
    Cfenv2 = syncWaveforms(Cfenv2,0,0,true);

    Cfenv2 = cutWaveforms(Cfenv2,dateshift(Cfenv2(1).ref,'start','day')+hours(0)+minutes(0),...
        00,hours(4));
    Ccut2 = cutWaveforms(Cf,dateshift(Cf(1).ref,'start','day')+hours(0)+minutes(0),...
        00,hours(4));
    kstnmsTmp = pull(Ccut2,'kstnm');
    [stlaTmp,stloTmp] = metaDataFromStationList(kstnmsTmp);

    winlen = secDur*newFs2;
    d2 = pull(Cfenv2);
    t2 = getTimeVec(Cfenv2);
    newD2 = [];
    for p = 1:size(d2,2)
        d2_ = d2(:,p);
        [d2_,~,endIndex] = cutWindows(d2_,winlen,0.5,false);
        d2_ = d2_(end,:)';
        newD2 = cat(2,newD2,d2_);
    end
    d2 = newD2;
    clear newD2 d2_;
    t2 = t2(endIndex);
    %ref = t2(1);
    %tref = dateshift(ref,'end','minute');
    %iStart = t2i(tref,ref,1/newFs2);
    %d2 = d2(iStart:winlen:end,:);
    %t2 = tref + seconds(secDur*(0:npts-1)');

    npts = size(d2,1);

    dI = d2 <= 1 | ~isfinite(d2);
    d2(dI) = 1;
    ObservedAmplitudes = d2;
    CorrectedAmplitudes = d2;

    lats = NaN(size(d2,1),1);
    lons = lats;
    minErr = lats;
    depths = lats;
    MeanCorrectedAmplitude = lats;
    tic;

    for j = 1:size(d2,1) %loop over time
        observedAmplitudes_ = ObservedAmplitudes(j,:);
        logCorrAmplitudes = NaN(lX,nStations);
        for k = 1:nStations
            d_ = distsCoarse(:,k);
            observedAmplitudes__ = observedAmplitudes_(k);
            logCorrAmplitudes(:,k) = log10(observedAmplitudes__) + ...
                (n*log10(d_) + ...
                mbest(1)*d_ + ...
                gamma1 + ...
                stationCorrections(k));
        end

        %
        badAmpsI = ~isfinite(logCorrAmplitudes);
        nGood = sum(lia) - sum(badAmpsI,2);
        meanLogCorrectedAmp = mean(logCorrAmplitudes,2,"omitnan");
        [~,maxMI] = max(abs(logCorrAmplitudes-meanLogCorrectedAmp),[],2,"omitnan");
        for k = 1:size(logCorrAmplitudes,1)
            if nGood(k) > 5
                logCorrAmplitudes(k,maxMI(k)) = NaN;
            end
        end

        badAmpsI = ~isfinite(logCorrAmplitudes);
        nGood = sum(lia)-sum(badAmpsI,2);
        meanLogCorrectedAmp = mean(logCorrAmplitudes,2,"omitnan");
        [~,maxMI] = max(abs(logCorrAmplitudes-meanLogCorrectedAmp),[],2,"omitnan");
        for k = 1:size(logCorrAmplitudes,1)
            if nGood(k) > 5
                logCorrAmplitudes(k,maxMI(k)) = NaN;
            end
        end

        meanLogCorrectedAmp = mean(logCorrAmplitudes,2,"omitnan");
        magErrTmp = mad(logCorrAmplitudes,0,2);
        magErrTmp1 = magErrTmp;

        [minErr_,minI] = min(magErrTmp);
        minErr1 = minErr_;
        lats(j) = Y(minI);
        lons(j) = X(minI);
        depths(j) = Z(minI);
        minErr(j) = minErr_;
        MeanCorrectedAmplitude(j) = 10.^meanLogCorrectedAmp(minI);
        CorrectedAmplitudes(j,:) = 10.^logCorrAmplitudes(minI,:);

        %% add finer grid search here
        X2 = lons(j)+Xf;
        Y2 = lats(j)+Yf;
        Z2 = depths(j)+Zf;

        lX2 = length(X2);
        uniqX2 = unique(X2);
        uniqY2 = unique(Y2);
        distsCoarse2 = NaN(lX2,nStations);
        for k = 1:nStations
            stla_ = stla(k);
            stlo_ = stlo(k);
            stel_ = stel(k);
            for l = 1:length(uniqX2)
                x_ = uniqX2(l);
                for m = 1:length(uniqY2)
                    y_ = uniqY2(m);
                    horI = X2 == x_ & Y2 == y_;
                    d_ = distance(y_,x_,stla_,stlo_,refEllipse)*1e-3;
                    distsCoarse2(horI,k) = sqrt(d_.^2 + ((Z2(horI) - stel_)/1000).^2);
                end
            end
        end

        observedAmplitudes_ = ObservedAmplitudes(j,:);
        logCorrAmplitudes = NaN(lX2,nStations);
        for k = 1:nStations
            d_ = distsCoarse2(:,k);
            observedAmplitudes__ = observedAmplitudes_(k);
            logCorrAmplitudes(:,k) = log10(observedAmplitudes__) + ...
                (n*log10(d_) + ...
                mbest(1)*d_ + ...
                gamma1 + ...
                stationCorrections(k));
        end

        %
        badAmpsI = ~isfinite(logCorrAmplitudes);
        nGood = sum(lia) - sum(badAmpsI,2);
        meanLogCorrectedAmp = mean(logCorrAmplitudes,2,"omitnan");
        [~,maxMI] = max(abs(logCorrAmplitudes-meanLogCorrectedAmp),[],2,"omitnan");
        for k = 1:size(logCorrAmplitudes,1)
            if nGood(k) > 5
                logCorrAmplitudes(k,maxMI(k)) = NaN;
            end
        end

        badAmpsI = ~isfinite(logCorrAmplitudes);
        nGood = sum(lia)-sum(badAmpsI,2);
        meanLogCorrectedAmp = mean(logCorrAmplitudes,2,"omitnan");
        [~,maxMI] = max(abs(logCorrAmplitudes-meanLogCorrectedAmp),[],2,"omitnan");
        for k = 1:size(logCorrAmplitudes,1)
            if nGood(k) > 5
                logCorrAmplitudes(k,maxMI(k)) = NaN;
            end
        end

        meanLogCorrectedAmp = mean(logCorrAmplitudes,2,"omitnan");
        magErrTmp = mad(logCorrAmplitudes,0,2);

        [minErr_,minI] = min(magErrTmp);
        lats(j) = Y2(minI);
        lons(j) = X2(minI);
        depths(j) = Z2(minI);
        minErr(j) = minErr_;

        MeanCorrectedAmplitude(j) = 10.^meanLogCorrectedAmp(minI);
        CorrectedAmplitudes(j,:) = 10.^logCorrAmplitudes(minI,:);
        disp(j);

        %%
        frameFlag = ~true;
        if frameFlag
            [minErr_,minI] = min(magErrTmp1);
            minErr1 = minErr_;
            lats(j) = Y(minI);
            lons(j) = X(minI);
            depths(j) = Z(minI);
            fix_xI = X == lons(j);
            fix_yI = Y == lats(j);
            fix_zI = Z == depths(j);

            close all;
            clear ax;
            fig = figure('units','normalized','outerposition',[0 0 0.55 1]);
            fig.Visible = 'off';
            ax(1) = subplot(421);
            scatter(X(fix_yI),Z(fix_yI),[],magErrTmp1(fix_yI),'o');
            zoom on; grid on; cinvis = colorbar; cinvis.Visible = 'off';
            hold(ax(1),'on'); plot(ax(1),lons(j),depths(j),'p','markerfacecolor','w');
            %xlabel('e-w (lon.)'); ylabel('elevation [m.]');

            ax(2) = subplot(422);
            scatter(Y(fix_xI),Z(fix_xI),[],magErrTmp1(fix_xI),'o'); zoom on; grid on; colorbar;
            hold(ax(2),'on'); plot(ax(2),lats(j),depths(j),'p','markerfacecolor','w');
            ax(2).XDir = 'reverse';
            xlabel('n-s (lat.)');  %ylabel('elevation [m.]');

            ax(3) = subplot(4,2,[3 4 5 6]);
            scatter(X(fix_zI),Y(fix_zI),[],magErrTmp1(fix_zI),'o'); zoom on; grid on; colorbar;
            hold(ax(3),'on'); plot(ax(3),lons(j),lats(j),'p','markerfacecolor','w');
            ylabel('n-s (lat.)'); xlabel('e-w (lon.)');

            minLat = min(Y);
            maxLat = max(Y);
            minLon = min(X);
            maxLon = max(X);
            axis(ax(3),[minLon maxLon minLat maxLat])
            axis(ax(3),'equal');

            [minErr_,minI] = min(magErrTmp);
            lats(j) = Y2(minI);
            lons(j) = X2(minI);
            depths(j) = Z2(minI);

            fix_xI2 = X2 == lons(j);
            fix_yI2 = Y2 == lats(j);
            fix_zI2 = Z2 == depths(j);

            scatter(ax(1),X2(fix_yI2),Z2(fix_yI2),[],magErrTmp(fix_yI2),'filled');
            zoom on; grid on; plot(ax(1),lons(j),depths(j),'p','markerfacecolor','r');
            set(ax(1),'ColorScale','log'); clim(ax(1),[0.005 0.5]);

            scatter(ax(2),Y2(fix_xI2),Z2(fix_xI2),[],magErrTmp(fix_xI2),'filled');
            zoom on; grid on; plot(ax(2),lats(j),depths(j),'p','markerfacecolor','r');
            set(ax(2),'ColorScale','log'); clim(ax(2),[0.005 0.5]);

            ylim(ax(1),[-5000 6000]);
            ylim(ax(2),[-5000 6000]);

            %scatter(ax(3),X2(fix_zI2),Y2(fix_zI2),[],magErrTmp(fix_zI2),'filled'); zoom on; grid on; colorbar;
            %plot(ax(3),lons(j),lats(j),'p','markerfacecolor','r');
            plot(ax(3),stloTmp,stlaTmp,'kv','markerfacecolor','w');
            plot(ax(3),center_lon,center_lat,'kd','markerfacecolor','m');
            scatter(ax(3),X2(fix_zI2),Y2(fix_zI2),[],magErrTmp(fix_zI2),'filled'); zoom on; grid on; colorbar;
            plot(ax(3),lons(j),lats(j),'p','markerfacecolor','r');
            ylabel('n-s (lat.)'); xlabel('e-w (lon.)');
            set(ax(3),'ColorScale','log');
            clim(ax(3),[0.005 0.5]);
            title(sprintf("Depth: %d",depths(j)));

            ax(4) = subplot(4,2,[7 8]);
            tdum = getTimeVec(Ccut2);
            ddum = Ccut2(1).d;
            plot(tdum,ddum,'k');
            zoom on; grid on; hold on;
            c4 = colorbar;
            c4.Visible = 'off';
            dumI = tdum >= t2(j) - minutes(1) & tdum < t2(j);
            plot(tdum(dumI),ddum(dumI),'r');

            titleStr = sprintf("%s",string(t2(j)));
            sgtitle(titleStr);
            frameStr = sprintf("/home/shernandez/research/now/cotopaxi/frames/frame_%04d.jpg",j);
            print('-djpeg',frameStr);
            pause(1);
            %continue;
        end
    end
    toc;
    CorrectedAmplitudes(dI) = NaN;

    %%
    fprintf("day: %s, number of stations: %d\n",datestr(day_),nStations);
    M = [M; MeanCorrectedAmplitude];
    magErr = [magErr; minErr];
    t = [t; t2];

    Lats = [Lats; lats];
    Lons = [Lons; lons];
    Depths = [Depths; depths];

    d2_ = NaN(size(d2,1),7);
    d2_(:,locb(lia)) = d2;
    d = [d; d2_];

    Mcorr2_ = NaN(size(CorrectedAmplitudes,1),7);
    Mcorr2_(:,locb(lia)) = CorrectedAmplitudes;
    correctedAmps = [correctedAmps; Mcorr2_];

    %%
    if plotFlag
        distFromSummit = distance(Lats,Lons,center_lat,center_lon,refEllipse)*1e-3;
        figure(); semilogy(t,M,'.'); zoom on;
        figure(); semilogy(t,magErr,'.'); zoom on;
        [center_lat,center_lon] = get_region_dimensions("cotopaxi");

        figure();
        plot(t,deg2km(Lons-center_lon),'.');
        zoom on; grid on; hold on;
        plot(t,medfiltSH(deg2km(Lons-center_lon),61,true),'-','linewidth',3);
        ylabel("Lon.");

        figure();
        plot(t,deg2km(Lats-center_lat),'.');
        zoom on; grid on; hold on;
        plot(t,medfiltSH(deg2km(Lats-center_lat),61,true),'-','linewidth',3);
        ylabel("Lat.");

        figure();
        plot(t,Depths,'.');
        zoom on; grid on; hold on;
        plot(t,medfiltSH(Depths,61,true),'-','linewidth',2);
        ylabel("Depth");

        figure();
        SS = scatter(Lons+randn(length(Lons),1)*0.0005,Lats+randn(length(Lats),1)*0.0005,...
            2*exp(log10(M)),Depths,'filled'); zoom on; grid on; colorbar; axis equal;
        SS.MarkerFaceAlpha = 0.5;
        SS.MarkerEdgeColor = 'k';
        SS.MarkerEdgeAlpha = 0.5;
        xlabel("Lon."); ylabel("Lat.");
        hold on; plot(center_lon,center_lat,'rv');

        figure();
        SS = scatter(Lons+randn(length(Lons),1)*0.0005,Lats+randn(length(Lats),1)*0.0005,...
            2*exp(log10(M)),M,'filled'); zoom on; grid on; colorbar; axis equal;
        SS.MarkerFaceAlpha = 0.5;
        SS.MarkerEdgeColor = 'k';
        SS.MarkerEdgeAlpha = 0.5;
        set(gca,'ColorScale','log');
        xlabel("Lon."); ylabel("Lat.");
        hold on; plot(center_lon,center_lat,'rv');

        figure();
        SS = scatter(Lons+randn(length(Lons),1)*0.0005,Lats+randn(length(Lats),1)*0.0005,...
            2*exp(log10(M)),datenum(t),'filled');
        zoom on; grid on; colorbar; axis equal;
        SS.MarkerFaceAlpha = 0.5;
        SS.MarkerEdgeColor = 'k';
        SS.MarkerEdgeAlpha = 0.5;
        xlabel("Lon."); ylabel("Lat.");
        hold on; plot(center_lon,center_lat,'rv');

        figure();
        SS3 = scatter3(deg2km(Lons+randn(length(Lons),1)*0.0005-center_lon),...
            deg2km(Lats+randn(length(Lats),1)*0.0005-center_lat),...
            (Depths+randn(length(Depths),1))/1000,2*exp(log10(M)),...
            datenum(t),'filled'); zoom on; grid on; colorbar; axis equal;
        SS3.MarkerFaceAlpha = 0.5;
        SS3.MarkerEdgeColor = 'k';
        SS3.MarkerEdgeAlpha = 0.5;
        xlabel("Lon."); ylabel("Lat.");
        hold on; plot3(0,0,5.9,'rv');

        magErr3 = std(log10(correctedAmps),0,2,"omitnan");
        badAmpsI = ~isfinite(correctedAmps);
        nBad = sum(badAmpsI,2);
        goodI = 7-nBad >= 5;
        %mI2 = magErr <= 0.08 & isfinite(magErr) & isfinite(M);
        mI2 = magErr <= 0.01 & isfinite(magErr) & isfinite(M) & ...
            mean(correctedAmps,2,"omitnan") >= 1e2 & distFromSummit <= 6.5 & ...
            t2 >= datetime(2023,05,12,09,00,00) & t2 < datetime(2023,05,12,15,00,00);

        figure('units','normalized','outerposition',[0 0 1/2 1]);
        plot(nBad,'.'); zoom on;

        figure('units','normalized','outerposition',[0 0 1/2 1]);
        semilogy(t(goodI),magErr3(goodI).*mad(log10(correctedAmps(goodI,:)),0,2),'.');
        zoom on; grid on; hold on;
        semilogy(t(goodI),mad(correctedAmps(goodI,:),0,2),'.');
        semilogy(t(goodI),magErr3(goodI),'.');

        figure('units','normalized','outerposition',[0 0 1/2 1]);
        clear ax;
        ax(1) = subplot(311);
        semilogy(t(mI2),correctedAmps(mI2,:),'.'); zoom on; grid on; hold on;
        semilogy(t(mI2),mean(correctedAmps(mI2,:),2,"omitnan"),'ko','markersize',14);
        zoom on; grid on;
        ax(2) = subplot(312);
        plot(t(mI2),magErr(mI2),'.'); %blue
        zoom on; grid on; hold on;
        plot(t(~mI2),magErr(~mI2),'.'); %red
        ax(3) = subplot(313);
        plot(getTimeVec(Cf),Cf(1).d);
        zoom on; grid on; hold on;
        linkaxes(ax,'x');
        axis tight;

        clear ax;
        figure('units','normalized','outerposition',[0 1/4 1 1/2]);
        ax(1) = subplot(121);
        semilogy(t2,ObservedAmplitudes,'.');
        zoom on; grid on;
        ax(2) = subplot(122);
        semilogy(t2,CorrectedAmplitudes,'.');
        zoom on; grid on;
        linkaxes(ax,'x');
        disp(sum(mI2));

        figure();
        plotLons = Lons(mI2);
        plotLats = Lats(mI2);
        plotDepths = Depths(mI2);
        SS4 = scatter3(deg2km(plotLons+randn(length(plotLons),1)*0.0005-center_lon),...
            deg2km(plotLats+randn(length(plotLats),1)*0.0005-center_lat),...
            (plotDepths+randn(length(plotDepths),1))/1000,2*exp(log10(M(mI2))),plotDepths,'filled');
        zoom on; grid on; colorbar; axis equal;
        SS4.MarkerFaceAlpha = 0.5; SS4.MarkerEdgeColor = 'k'; SS4.MarkerEdgeAlpha = 0.5;
        xlabel("Lon."); ylabel("Lat.");
        hold on; plot3(0,0,5.9,'rv');
    end
end

if saveFlag
    cd ~/research/now/cotopaxi/;
    correctedAmpsOrig = correctedAmps;
    save("CotopaxiTremorLocationResults_v1.4.1.mat",'M','magErr','t','d',...
        'correctedAmps','Lats','Lons','Depths','correctedAmpsOrig');
end
