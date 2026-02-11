clear; close all; clc;
srcLat = -5.807;
srcLon = -75.264;
srcT = datetime(2019,05,26,07,41,15);

stnmList = ["BONI","PIAT","CASC","BV15","ISPT","BTER"];

%%
lfc = 1/160;
hfc = 1/20;
npoles = 4;
units = 'disp';

lS = length(stnmList);
[wa_poles,wa_zeros,wa_sensitivity] = readSensorPolesZeros(10);
S = populateSacStruct(lS);
parfor i = 1:lS
    disp(i);
    kstnm_ = stnmList(i);
    S(i) = loadSacData(dn2dt(floor(datenum(srcT))),1,kstnm_,"HHZ");
    [zeros,poles,constant] = getPolesZeros(S(i));
    S(i) = transfer(S(i),zeros,poles,constant,lfc,hfc,npoles,units,1); % 1 = deconvolve
    %S(i) = transfer(S(i),wa_zeros,wa_poles,wa_sensitivity,-inf,-inf,npoles,'disp',0); % 1 = deconvolve
end
S = scaleSacData(S,1e3);

%%
Scut = cutSacData(S,srcT,0,10*60);
[stla,stlo,stel] = metaDataFromStationList(stnmList);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(srcLat,srcLon,stla,stlo,refEllipse)*1e-3;
[dSort,dI] = sort(d_);
Scut = Scut(dI);
stnmList = stnmList(dI);

close all;
for i = 1:lS
    ax = plotSacData(Scut(i));
end
ax = plotSacData(Scut);

%%
Scut = synchSacData(Scut);
d = pull(Scut);
maxAmp = max(d);
%d = normalizeTraces(d,true,true);
ampFact = 40;
dOrig = d;
for i = 1:lS
    d(:,i) = dSort(i) + ampFact*d(:,i)/max(maxAmp);
end
t = getTimeVec(Scut);
figure(); plot(t,d,'linewidth',4); zoom on;

