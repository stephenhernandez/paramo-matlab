function [M,magErr,t,d,Mcorr,mbest,Cfenv2,Cf,C] = ...
    loopCotopaxiTremorLocation(tStart,tEnd,plotFlag,saveFlag)
if nargin < 1
    tStart = datetime(2023,01,13);
end

if nargin < 2
    tEnd = dateshift(datetime('now'),'start','day');
end

if nargin < 3
    plotFlag = false;
end

if nargin < 4
    saveFlag = true;
end

%%
Model = load('~/masa/old/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v15022023.mat');
%Model = load('~/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v11022023.mat');
%Model = load('~/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v30012023.mat');

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
fixFlag = Model.fixFlag;
gamma1 = Model.gamma1;
newFs2 = Model.newFs2;
charSNCL = char(allMySNCLs);
allMySNCLs = string(charSNCL(:,1:6));

if fixFlag
    stationCorrectionsOrig = mbest(2:end);
else
    stationCorrectionsOrig = mbest(3:end);
end

%%
M = [];
magErr = M;
t = M;
d = M;
Mcorr = M;
Cenv = [];

dayVec = (tStart:tEnd)';
lDays = length(dayVec);

for j = 1:lDays
    day_ = dayVec(j);
    C = loadWaveforms(day_,1,kstnmsOrig,...
        ["BHZ";"HHZ"],"EC");

    Cf = detrendWaveforms(C);
    Cf = syncWaveforms(Cf,false,true);
    Cf = nanGapWaveforms(Cf,0);
    Cf = resampleWaveforms(Cf,newFs);
    Cf = scaleWaveforms(...
        transferWaveforms(Cf,lfc,hfc,npoles,newFs,"vel",deconvolveFlag,waFlag),1e9);
    Cf = nanGapWaveforms(Cf,0);
    Cf = padWaveforms(Cf);

    Tsncls = allMySNCLs;
    Csncls = strcat(pull(Cf,'knetwk'),pull(Cf,'kstnm'));
    [lia,locb] = ismember(Csncls,Tsncls);

    if sum(lia) < 3
        continue;
    end

    Cf = Cf(lia);
    nStations = length(Cf);
    d_ = dOrig(locb(lia));
    stationCorrections = stationCorrectionsOrig(locb(lia));

    Cfenv2 = envelopeWaveforms(Cf);
    Cfenv2 = resampleWaveforms(Cfenv2,newFs2);
    Cfenv2 = medfiltWaveforms(Cfenv2,newFs2*secDur,false);
    Cfenv2 = syncWaveforms(Cfenv2);

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
            Mcorr2(:,i) = log10(Mcorr2(:,i)) + ...
                (n*log10(d_(i)) + mbest(1)*d_(i) + gamma1 + stationCorrections(i));
        else
            Mcorr2(:,i) = log10(Mcorr2(:,i)) + ...
                (mbest(1)*log10(d_(i)) + mbest(2)*d_(i) + gamma1 + stationCorrections(i));
        end
    end

    %
    Mcorr2(dI) = NaN;
    M2 = median(Mcorr2,2,"omitnan");
    magErr2 = mad(Mcorr2,1,2);

    if plotFlag
        mI2 = magErr2 <= 0.08 & isfinite(magErr2) & isfinite(M2);
        figure('units','normalized','outerposition',[0.05 0.05 0.8 1]);
        TL = tiledlayout(3,1,"Padding","compact","TileSpacing","none");
        clear ax;

        ax(1) = nexttile();
        M2 = median(Mcorr2,2,"omitnan");
        semilogy(t2,10.^Mcorr2,'.'); zoom on; grid on; hold on;
        semilogy(t2,median(10.^Mcorr2,2,"omitnan"),'k.','markersize',14);
        zoom on; grid on;

        ax(2) = nexttile();  % subplot(312);
        plot(t2(mI2),magErr2(mI2),'.');
        zoom on; grid on; hold on;
        plot(t2(~mI2),magErr2(~mI2),'.');

        ax(3) = nexttile(); % subplot(313);
        plot(getTimeVec(Cf),Cf(1).d);
        zoom on; grid on; hold on;
        linkaxes(ax,"x");
        ax(1).XTickLabel = [];
        ax(2).XTickLabel = [];
        axis tight;

        disp(sum(mI2));
    end

    %%
    fprintf("day: %s, number of stations: %d\n",string(day_),nStations);
    M = [M; M2];
    magErr = [magErr; magErr2];
    t = [t; t2];

    d2_ = NaN(size(d2,1),7);
    d2_(:,locb(lia)) = d2;
    d = [d; d2_];

    Mcorr2_ = NaN(size(Mcorr2,1),7);
    Mcorr2_(:,locb(lia)) = Mcorr2;
    Mcorr = [Mcorr; Mcorr2_];
end

if saveFlag
    cd ~/masa/old/research/now/cotopaxi/;
    save("CotopaxiTremorLocationResults_v5",'M','magErr','t','d','Mcorr');
end