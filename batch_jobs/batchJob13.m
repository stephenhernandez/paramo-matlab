clear
cd ~/templateSearch/PINO/
files = dir('dayTemplateSearch*mat');
lfiles = length(files);
tMain = []; ccMain = []; ampMain = []; dMag = []; magMain = []; templateNumber = []; madMain = []; evidMain = [];
for i = 1:lfiles
    T = load(files(i).name);
    t_ = T.tMain;
    tMain = [tMain; t_];
    cc_ = T.ccMain;
    ccMain = [ccMain; cc_];
    amp_ = T.ampMain;
    ampMain = [ampMain; amp_];
    dmag_ = T.dMag;
    dMag = [dMag; dmag_];
    mag_ = T.magMain;
    magMain = [magMain; mag_];
    tnum_ = T.templateNumber;
    templateNumber = [templateNumber; tnum_];
    mad_ = T.madMain;
    madMain = [madMain; mad_];
    id_ = T.evidMain;
    evidMain = [evidMain; id_];
    clear t_ cc_ amp_ dmag_ mad_ mag_ tnum_ id_
end
[tMain,sI] = sort(tMain);
dMag = dMag(sI);
ampMain = ampMain(sI);
ccMain = ccMain(sI);
evidMain = evidMain(sI);
templateNumber = templateNumber(sI);
magMain = magMain(sI);
madMain = madMain(sI);
clear files i lfiles sI T
close all; [t,cc,amp_,dMag_,evid_,mag_,templateNumber_,mad_] = filterCatalog(tMain,ccMain,30,ampMain,dMag,evidMain,magMain,templateNumber,madMain);
tI = cc >= 0.18 & mad_ >= 8 & amp_ >= 80 & templateNumber_ == 1; % & t >= datetime(2023,01,01); % & t < datetime(2022,10,23);
figure(); nDays = 1; plot(t(tI),t2r(t(tI),days(nDays))/nDays,'.'); zoom on;
figure(); semilogy(t(tI),amp_(tI),'.'); zoom on;
figure(); plot(t(tI),cc(tI),'o'); zoom on; hold on; plot(t(~tI),cc(~tI),'.'); grid on;
figure(); plot(t(tI),1:sum(tI),'.'); zoom on; grid on;

A_ = extractWaveforms(sort(t(tI)),seconds(30),["BTER"],"HHZ","EC",[""],true,true,1,true,[1,16,false,true,false]);
clearvars -except A_
save('~/research/now/pichincha/BTER_mainRepeater_waveforms');
