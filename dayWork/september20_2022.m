clear; close all; clc;
%loadDay = datetime(2024,01,04); eqlat = -5.848; eqlon = -80.917; tCutStart = loadDay + hours(13)+minutes(30); tCutDuration = minutes(60);
%loadDay = datetime(2023,02,06); eqlat = 37.225; eqlon = 37.02; tCutStart = loadDay + hours(01); tCutDuration = minutes(120);
%loadDay = datetime(2024,01,01); eqlat = 37.498; eqlon = 137.242; tCutStart = loadDay + hours(07); tCutDuration = minutes(240);
%loadDay = datetime(2024,12,05); eqlat = 40.370; eqlon = -125.025; tCutStart = loadDay + hours(18)+minutes(50); tCutDuration = minutes(120);
%loadDay = datetime(2025,03,28); eqlat = 22.013; eqlon = 95.922; tCutStart = loadDay + hours(06)+minutes(20); tCutDuration = minutes(240);
loadDay = datetime(2025,07,29); eqlat = 52.530; eqlon = 160.165; tCutStart = loadDay + hours(23)+minutes(0); tCutDuration = hours(12);

%loadDay = datetime(2013,02,06); eqlat = -10.799; eqlon = 165.114; tCutStart = loadDay + hours(01); tCutDuration = minutes(120);
%loadDay = datetime(2021,07,29); eqlat = 55.364; eqlon = -157.888; tCutStart = loadDay + hours(06); tCutDuration = minutes(120);
%loadDay = datetime(2021,08,12); eqlat = -58.375; eqlon = -25.264; tCutStart = loadDay + hours(18); tCutDuration = minutes(120); 
%loadDay = datetime(2021,03,04); eqlat = -29.723; eqlon = -177.279; tCutStart = loadDay + hours(19); tCutDuration = minutes(150); 
%loadDay = datetime(2023,05,10); eqlat = -15.628; eqlon = -174.493; tCutStart = loadDay + hours(16); tCutDuration = minutes(150); 
%loadDay = datetime(2015,04,15); eqlat = 32.791; eqlon = 130.754; tCutStart = loadDay + hours(16); tCutDuration = minutes(120);
%loadDay = datetime(2014,04,01); eqlat = -19.610; eqlon = -70.769; tCutStart = loadDay + hours(23); tCutDuration = minutes(120);
%loadDay = datetime(2022,09,19); eqlat = 18.367; eqlon = -103.252;
%loadDay = datetime(2017,09,08);

S = seedData2(loadDay,tCutStart,tCutStart+tCutDuration,false,true,["BHZ";"HHZ"],[],[],"~/rawdata_cotopaxi");

%%
%S = cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(01)+minutes(15),0,hours(3/4));
%S = cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(18)+minutes(00),0,hours(3/4));
%S = cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(04)+minutes(50),0,hours(3/4));

badI = isnat(pull(S,'ref')) | strcmp(pull(S,'kcmpnm'),"SHZ");
S(badI) = [];

%knets = pull(S,'knetwk');
%goodI = strcmp(pull(S,'knetwk'),"EC") | strcmp(pull(S,'knetwk'),"CM") | strcmp(pull(S,'knetwk'),"OP") | strcmp(pull(S,'knetwk'),"C1");
%S(~goodI) = [];

badI = strcmp(pull(S,'kcmpnm'),"BDF");
S(badI) = [];

%badI = strcmp(pull(S,'kcmpnm'),"HNZ") | strcmp(pull(S,'kcmpnm'),"ENZ");
%S(badI) = [];

%%
refEllipse = referenceEllipsoid('wgs84');
newFs = 0.1;
rmhpSecDur = 100;
units = 'vel';
% Supdated = resampleWaveforms(...
%     detrendWaveforms(...
%     syncWaveforms(...
%     differentiateWaveforms(S))),100); %prefilt
Supdated = detrendWaveforms(...
    syncWaveforms(...
    rmhpWaveforms(...
    nanGapWaveforms(...
    taperWaveforms(...
    detrendWaveforms(differentiateWaveforms(S)),200*rmhpSecDur*2),0),...
    rmhpSecDur,2*200*rmhpSecDur),false,true,true)); %prefilt
Supdated = resampleWaveforms(Supdated,newFs);

%
Supdated = nanGapWaveforms(Supdated,0);
Sorig = Supdated;

e = pull(Supdated,'e');
refs = pull(Supdated,'ref');
npts = pull(Supdated,'npts');
deltas = pull(Supdated,'delta');
knetwks = pull(Supdated,'knetwk');
kstnms = pull(Supdated,'kstnm');
kholes = pull(Supdated,'khole');
kcmpnms = pull(Supdated,'kcmpnm');
tEnd = refs + e;

updatedSNCLs = strcat(knetwks,kstnms,kholes,kcmpnms);
lUpdatedSNCLs = length(updatedSNCLs);

%
Fs = 1./deltas;
FsI = Fs >= 1;
Fs(FsI) = round(Fs(FsI));
mySNCLs = [kstnms kcmpnms knetwks kholes];
%responseStructure = singleSNCLFreqResponse(mySNCLs,loadDay,loadDay+1,npts,Fs,units);
uniqFs = unique(Fs);

S2 = Supdated;
lfc = 1/500;
hfc = 1/200; %-inf;
tw = 100;
npoles = 4;
scaleFactor = 1e9;

lS = length(Supdated);
gapFlags = pull(S,'gapFlag'); 
Ssnr = zeros(lS,1); 
for i = 1:length(Supdated)
    npts_ = npts(i);
    S_ = Supdated(i);
    kcmpnm_ = kcmpnms(i);
    R = singleSNCLFreqResponse(mySNCLs(i,:),loadDay,loadDay+1,npts(i),Fs(i),units); %responseStructure(i);
    ref_ = refs(i);
    tEnd_ = tEnd(i);
    Fs_ = 1./S_.delta;
    Tresp = R.Tstart;
    if isnat(Tresp)
        continue;
    end

    Hresp = R.H;
    Hresp = Hresp(:,1);
    dOrig = S_.d;           % time-domain
    dOrig = detrend(taper(detrend(dOrig),tw));
    Dorig = fft(dOrig);     % frequency-domain
    Hbu = freqOperator(npts_,lfc,hfc,Fs_,npoles);
    D = Dorig.*Hresp.*Hbu;
    df = scaleFactor*ifft(D,'symmetric');
    S_ = dealHeader(S_,df);

    disp(i);
    S_.stel = R.Stel;
    S_.stlo = R.Stlo;
    S_.stla = R.Stla;
    S2(i) = S_;
    [locs,snr,staOlta,sosSTA] = stalta(S_,60,60,10,false,0,0,-inf,true);
    
    if ~isempty(snr)
        maxSNR = max(snr);
        Ssnr(i,1) = maxSNR;
    end
end

%%
close all;
stla = pull(S2,'stla');
stlo = pull(S2,'stlo');

S2orig = S2;
badI = ~isfinite(stla) | ~isfinite(stlo) | gapFlags > 0 | Ssnr < 10;
S2(badI) = [];
stla(badI) = [];
stlo(badI) = [];

%%
d_ = distance(stla,stlo,eqlat,eqlon,refEllipse);
d_ = d_*1e-3;

%
S2 = syncWaveforms(S2,~true,true,true);
d = pull(S2);
d(~isfinite(d)) = 0;
maxRMS = prctile(log10(rms(d)'),95);
minRMS = prctile(log10(rms(d)'),05);
badI = rms(d)' > 10.^maxRMS | rms(d)' < 10.^minRMS; % | d_ < 3000;
S2(badI) = [];

% run above code again
S2 = syncWaveforms(S2,~true,true,true);
d = pull(S2);
d(~isfinite(d)) = 0;
maxRMS = prctile(log10(rms(d)'),95);
minRMS = prctile(log10(rms(d)'),05);
badI = rms(d)' > 10.^maxRMS | rms(d)' < 10.^minRMS; % | d_ < 3000;
S2(badI) = [];

%
stla = pull(S2,'stla');
stlo = pull(S2,'stlo');

refEllipse = referenceEllipsoid('wgs84');
[d_,az] = distance(eqlat,eqlon,stla,stlo,refEllipse);
[~,baz] = distance(stla,stlo,eqlat,eqlon,refEllipse);
d_ = d_*1e-3;

[d_,sortI] = sort(d_);
S2 = S2(sortI);
stla = pull(S2,'stla');
stlo = pull(S2,'stlo');

%
S2 = nanGapWaveforms(S2,0);
d = pull(S2);
d(~isfinite(d)) = 0;
scaleFactor2 = 1./median(rms(d)');

%
close all;
figure('units','normalized','outerposition',[0 0 0.6 1]);
hold on;
plot(stlo,stla,'.'); zoom on; grid on; axis equal;
geodesicN = 301;
for i = 1:length(S2)
    latout = eqlat;
    lonout = eqlon;
    geolat = latout;
    geolon = lonout;
    stla_ = stla(i);
    stlo_ = stlo(i);
    az_ = az(i);
    d__ = d_(i);
    for j = 1:geodesicN
        [latout,lonout] = reckon(latout,lonout,1e3*d__/geodesicN,az_,refEllipse);
        geolat = [geolat; latout];
        geolon = [geolon; lonout];
        [~,az_] = distance(latout,lonout,stla_,stlo_,refEllipse);
    end
    latout = eqlat;
    lonout = eqlon;
    plot(geolon,geolat,'k-');
end

load ~/igdata/ec_boundaries.mat
load ~/igdata/soam_noec.mat;
hold on
plot(lonEC,latEC,'-','color',[0.5 0.5 0.5],'LineWidth',2);
plot(lon_noec,lat_noec,'-','color',[0.5 0.5 0.5],'LineWidth',2);
text(stlo,stla,pull(S2,'kstnm'),'FontSize',15); grid on;
axis([-82.5 -74.5 -6 2]);

%
figure('units','normalized','outerposition',[0 0 0.6 1]);
hold on;
scaleFactor2 = 3/median(rms(d)');
tdum = (0:size(d,1)-1)'/newFs; %getTimeVec(S2);
for i = 1:length(S2)
    dd = d(:,i);
    pp = plot(tdum,4*scaleFactor2*dd + d_(i),'linewidth',2);
    pp.Color(4) = 0.8;
end
zoom on;
text(repmat(max(tdum),length(d_),1),d_,pull(S2,'kstnm'),'FontSize',20); grid on;
ylabel('Distance from Epicenter [km.]');
grid on;

%
figure('units','normalized','outerposition',[0 0 0.6 1]);
plot(d_,az,'.'); zoom on;
grid on;
%xlim([0 10000]);
ylim([0 360]);

%
figure('units','normalized','outerposition',[0 0 0.6 1]);
hold on;
scaleFactor2 = 2000./median(rms(d)');
for i = 1:length(S2)
    dd = zpkFilter(log10(abs(hilbert(d(:,i)))),-inf,1/100,1,2,1);
    pp = plot(tdum,16*scaleFactor2*dd + d_(i),'linewidth',2);
    pp.Color(4) = 0.8;
end
zoom on;
text(repmat(max(tdum),length(d_),1),d_,pull(S2,'kstnm'),'FontSize',20); grid on;
ylabel('Distance from Epicenter [km.]');
grid on;

%% code valid for september 2022 mexican earthquake
S3 = cutWaveforms(S2,dateshift(S2(1).ref,'start','day')+hours(18)+minutes(53)+seconds(00),0,seconds(120));
d4 = pull(S3);

[shiftedData1,~,~,~,rawShifts1] = apply_vdcc(d4,[],0,true,true);
d2 = d4;
tmp_stack = plot_family(d2,1:size(d2,2),30,newFs);
d2 = shiftedData1;
shiftedStack1 = plot_family(d2,1:size(d2,2),30,newFs);

[xcorr1,lags1] = doCrossCorrFreqDom(repmat(normalizeWaveforms(shiftedStack1(25*newFs:50*newFs,:)),1,size(d4,2)),...
    normalizeWaveforms(shiftedData1(25*newFs:50*newFs,:)));
[~,maxXCI1] = max(abs(xcorr1));
rawShifts2 = round(lags1(maxXCI1));
signXC = sign(xcorr1);
flipFlag = false(size(d4,2),1);
for i = 1:size(d4,2)
    signxc_ = signXC(:,i);
    if signxc_(maxXCI1(i)) == -1
        flipFlag(i) = true;
        d4(:,i) = -d4(:,i);
    end
end

figure(5);
hold on;
plot(rawShifts1+rawShifts2,'.'); grid on; zoom on;

shiftedData2 = apply_shifts(d4,rawShifts1+rawShifts2);
d2 = shiftedData2;
shiftedStack2 = plot_family(d2,1:size(d2,2),30,newFs);
[xcorr2,lags2] = doCrossCorrFreqDom(repmat(normalizeWaveforms(shiftedStack2(30*newFs:50*newFs,:)),1,size(d4,2)),...
    normalizeWaveforms(shiftedData2(30*newFs:50*newFs,:)));
[~,maxXCI2] = max(abs(xcorr2));
rawShifts3 = round(lags2(maxXCI2));
signXC = sign(xcorr2);
flipFlag2 = false(size(d4,2),1);
for i = 1:size(d4,2)
    signxc_ = signXC(:,i);
    if signxc_(maxXCI2(i)) == -1
        flipFlag2(i) = true;
        d4(:,i) = -d4(:,i);
    end
end

figure(5);
hold on;
plot(rawShifts1+rawShifts2+rawShifts3,'.'); grid on; zoom on;

shiftedData3 = apply_shifts(d4,rawShifts1+rawShifts2+rawShifts3);
d2 = shiftedData3;
shiftedStack3 = plot_family(d2,1:size(d2,2),30,newFs);
figure();
plot([rawShifts1 rawShifts2 rawShifts3],'.');
zoom on;

[xcorr3,lags3] = doCrossCorrFreqDom(repmat(normalizeWaveforms(shiftedStack3(30*newFs:50*newFs,:)),1,size(d4,2)),...
    normalizeWaveforms(shiftedData2(30*newFs:50*newFs,:)));
[maxXC3,maxXCI3] = max(abs(xcorr3));

%
d2 = apply_shifts(d,rawShifts1+rawShifts2+rawShifts3);
figure('units','normalized','outerposition',[0 0 0.6 1]);
Nroot = 1;
d2 = normalizeWaveforms(d2,true,true);
Vroot = mean((abs(d2).^(1/Nroot)).*(sign(d2)),2,"omitnan");
plot(tdum,(abs(Vroot).^Nroot).*(sign(Vroot)),'linewidth',2); zoom on;

figure(11);
hold on;
Nroot = 2;
d2 = apply_shifts(d,rawShifts1+rawShifts2+rawShifts3);
d2 = normalizeWaveforms(d2,true,true);
Vroot = mean((abs(d2).^(1/Nroot)).*(sign(d2)),2,"omitnan");
plot(tdum,(abs(Vroot).^Nroot).*(sign(Vroot)),'linewidth',2);
zoom on;

Nroot = 3;
d2 = apply_shifts(d,rawShifts1+rawShifts2+rawShifts3);
d2 = normalizeWaveforms(d2,true,true);
Vroot = mean((abs(d2).^(1/Nroot)).*(sign(d2)),2,"omitnan");
plot(tdum,(abs(Vroot).^Nroot).*(sign(Vroot)),'linewidth',2);
zoom on;

Nroot = 4;
d2 = apply_shifts(d,rawShifts1+rawShifts2+rawShifts3);
d2 = normalizeWaveforms(d2,true,true);
Vroot = mean((abs(d2).^(1/Nroot)).*(sign(d2)),2,"omitnan");
plot(tdum,(abs(Vroot).^Nroot).*(sign(Vroot)),'linewidth',2);
zoom on;

plot(tdum,pws(normalizeWaveforms(apply_shifts(d,rawShifts1+rawShifts2+rawShifts3),true,true)),'linewidth',2);
zoom on;
