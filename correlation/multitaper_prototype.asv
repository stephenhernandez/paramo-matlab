clear; close all; clc;

secDur = 600;
newFs = 10; %64;
totN = secDur*newFs;
nfft = 2^(nextpow2(totN)+1); %*2;
nOverlap = 6/8;
df = newFs/nfft;
f1 = (0:nfft-1)'*df;
fshift = (-nfft/2:nfft/2-1)'*(newFs/nfft);
fOrig = fftshift(fshift,1);
nw = 4;
lfc = 0.1;
hfc = 0.2;
npoles = 4;
zeroPhaseFlag = true;
nu = 4;

%rawdatadir = '/Users/stephen/data/galapagos/iguana/BROADBAND';
rawdatadir = '/Volumes/KINGSTON/data/galapagos/iguana/BROADBAND/';
%dayStart = datetime(2018,05,01);
%dayEnd = datetime(2019,03,10);
%dayStart = datetime(2023,07,15);
%dayEnd = datetime(2023,07,15);
%dayStart = datetime(2024,08,23);
%dayEnd = datetime(2024,08,23);
% dayStart = datetime(2013,08,23);
% dayEnd = datetime(2013,08,23);
dayStart = datetime(2018,01,280);
dayEnd = datetime(2018,01,280);
% dayStart = datetime(2018,06,26);
% dayEnd = datetime(2018,06,26);
%dayStart = datetime(2018,08,11);
%dayEnd = datetime(2018,08,11);
% dayStart = datetime(2024,10,21);
% dayEnd = datetime(2024,10,21);
%dayStart = datetime(2023,05,12);
%dayEnd = datetime(2023,05,12);
% dayStart = datetime(2024,07,09);
% dayEnd = datetime(2024,07,09);
%dayStart = datetime(2015,08,14);
%dayEnd = datetime(2015,08,14);
%dayStart = datetime(2022,10,22);
%dayEnd = datetime(2022,10,22);
%dayStart = datetime(2018,04,25);
%dayEnd = datetime(2018,04,25);
% dayStart = datetime(2021,12,02);
% dayEnd = datetime(2021,12,02);
% dayStart = datetime(2022,07,31);
% dayEnd = datetime(2022,07,31);
% dayStart = datetime(2020,01,12);
% dayEnd = datetime(2020,01,12);

dayInc = 3;
dayVec = (dayStart:dayInc:dayEnd)';
lDays = length(dayVec);
dayDecon = NaN(nfft,lDays);
plotFlag = true;

%S = loadWaveforms(datetime(2024,08,30),1,["BREF";"PUYO"],["BHZ";"HHZ"],"EC");
%S = loadWaveforms(datetime(2024,03,05),1,["BREF";"CASC"],["BHZ";"HHZ"],"EC");
for i = 1:lDays
    day_ = dayVec(i);
    %S = loadWaveforms(day_,dayInc,["SN14";"SN11"],"HHN","9D","",true,true,rawdatadir);
    S = loadWaveforms(day_,dayInc,["FER2";"FER1"],["BHZ";"HHZ"],"9D","",true,true,rawdatadir);
    %S = loadWaveforms(day_,1,["BREF";"CASC"],["BHZ";"HHZ"],"EC");
    %S = loadWaveforms(day_,1,["BREF";"BNAS"],"BHZ","EC");
    %S = loadWaveforms(day_,5,"PINO",["HHZ";"SHZ"],"EC",""); S(2) = scaleWaveforms(S(2),-1);
    %S = loadWaveforms(day_,5,"PINO",["SHZ";"HHZ"],"EC",""); S(1) = scaleWaveforms(S(1),-1);
    %S = loadWaveforms(day_,1,["BREF";"CO1V"],["BHZ";"HHZ"],"EC",""); S(1) = scaleWaveforms(S(1),-1);
    % S = loadWaveforms(day_,2,"VCH1",["BDF";"HHZ"],"EC",["";"01"]);
    % S = (detrendWaveforms(cutWaveforms((S),dateshift(S(1).ref,'start','day')+hours(12)+minutes(0)+seconds(0),0,hours(24))));
    % S = loadWaveforms(day_,2,"PIKA",["HDF";"HHZ"],"EC",["";"01"]);
    % S = detrendWaveforms(cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(0)+minutes(0)+seconds(0),0,hours(26)));
    % S = loadWaveforms(day_,1,"RUNA",["BDF";"HHZ"],"EC",["";"01"]);
    % S = detrendWaveforms(cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(0)+minutes(0)+seconds(0),0,hours(24)));
    % S = loadWaveforms(day_,2,"BNAS",["BDF";"BHZ"],"EC",["";"01"]);
    % S = detrendWaveforms(cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(12)+minutes(0)+seconds(0),0,hours(3)));

    %S = loadWaveforms(day_,2,"SAG1",["BDF";"HHN"],"EC",["";"05"]);
    %S = detrendWaveforms(cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(0)+minutes(0)+seconds(0),0,hours(24*2)));

    %S = loadWaveforms(day_,1,"BNAS",["BDF";"HHZ"],"EC",["";"01"]);
    %S = (detrendWaveforms(cutWaveforms((S),dateshift(S(1).ref,'start','day')+hours(0)+minutes(0)+seconds(0),0,hours(6))));
    %S = loadWaveforms(day_,2,"REVS",["BDF";"HHZ"],"EC",["";"01"]);
    %S = (detrendWaveforms(cutWaveforms((S),dateshift(S(1).ref,'start','day')+hours(20)+minutes(0)+seconds(0),0,hours(24))));
    % S = loadWaveforms(day_,1,["SAGA"],["HHZ";"BDF"],"EC",["";"01"]);
    % S = (detrendWaveforms(cutWaveforms((S),dateshift(S(1).ref,'start','day')+hours(2)+minutes(0)+seconds(0),0,hours(10))));
    %S = loadWaveforms(day_,dayInc,"PINO","SHZ","EC",""); %,true,true,rawdatadir); S(2,1) = S;
    [pws_stack,decon_fd,Sxx,Syy,Sxy,xcohere] = process_(S,newFs,totN,nfft,nOverlap,fshift,nw,lfc,hfc,npoles,zeroPhaseFlag,nu,plotFlag);
    dayDecon(:,i) = pws_stack;
end
m = size(Sxy,1)/2;
lags = (-m+1:m-1)';
lags = lags/newFs;

function [pws_stack,decon,Sxx,Syy,Sxy,xcohere] = process_(S,newFs,totN,nfft,nOverlap,fshift,nw,lfc,hfc,npoles,zeroPhaseFlag,nu,plotFlag)
S = nanGapWaveforms(resampleWaveforms(detrendWaveforms(S),newFs),0);
S = nanGapWaveforms(syncWaveforms((S),true,true,true),0);

tOrig = getTimeVec(S);
detrendFlag = true;
[dcut1,~,endIndex] = cutWindows(S(1).d,totN,nOverlap,detrendFlag);
dcut2 = cutWindows(S(2).d,totN,nOverlap,detrendFlag);

%
[npts,ntimewin] = size(dcut1);
[dpsSeq,lambda] = dpss(npts,nw,2*nw-1);
w = lambda./sum(lambda);
Sxx = NaN(nfft,size(dcut1,2));
Syy = Sxx;
Sxy = Sxx;
for i = 1:size(dcut1,2)
    dcut_1 = detrend(dcut1(:,i));
    dcut_1 = dcut_1.*dpsSeq;
    dcut_2 = detrend(dcut2(:,i));
    dcut_2 = dcut_2.*dpsSeq;
    D_1 = fft(dcut_1,nfft);
    D_2 = fft(dcut_2,nfft);
    D_1 = w'.*D_1;
    D_2 = w'.*D_2;
    Sxx_1 = sum(D_1.*conj(D_1),2);
    Sxx(:,i) = Sxx_1;
    Syy_1 = sum(D_2.*conj(D_2),2);
    Syy(:,i) = Syy_1;
    Sxy_ = sum(conj(D_1).*D_2,2);
    Sxy(:,i) = Sxy_;
end
tdown = tOrig(endIndex);

%%
tic;
wienerWindow = [7,3];
SxxSmooth = wiener2(Sxx,wienerWindow);
SyySmooth = wiener2(Syy,wienerWindow);
xcohere = Sxy./(sqrt(Sxx).*sqrt(Syy));
%xcohere = Sxy./(sqrt(SxxSmooth).*sqrt(SyySmooth));
wldFlag = false;
if wldFlag
    wlthresh = 0.1;
    maxS = median(abs(Sxx),1); %,[],1);
    wl = wlthresh*maxS;
    Sxx_ = Sxx;
    for j = 1:ntimewin
        Sxx_1 = Sxx_(:,i);
        Sxx_1(Sxx_1<=wl(i)) = wl(i);
        Sxx_(:,i) = Sxx_1;
    end
    %decon = xcohere;
    decon = Sxy./Sxx_;
else
    %decon = Sxy./Sxx;
    decon = xcohere;
    %decon = Sxy./SxxSmooth;
end
msc = xcohere.*conj(xcohere);
mscSmooth = wiener2(msc,wienerWindow);

close all;
clear ax;

decon_dB = 10*log10(decon.*conj(decon));
if plotFlag
    figure(); imagesc(tdown,fshift,fftshift(10*log10(SxxSmooth),1)); zoom on; colorbar; title(sprintf('%s.%s psd',S(1).kstnm,S(1).kcmpnm)); axis xy; ylim([0 newFs/2]); ax(1) = gca;
    figure(); imagesc(tdown,fshift,fftshift(10*log10(SyySmooth),1)); zoom on; colorbar; title(sprintf('%s.%s psd',S(2).kstnm,S(2).kcmpnm)); axis xy; ylim([0 newFs/2]); ax(2) = gca;
    figure(); imagesc(tdown,fshift,fftshift(10*log10(Sxy.*conj(Sxy)),1)); zoom on; colorbar; title('cross psd'); axis xy; ylim([0 newFs/2]); ax(3) = gca;
    figure(); imagesc(tdown,fshift,fftshift(decon_dB,1)); zoom on; colorbar; title('decon psd'); axis xy; ylim([0 newFs/2]); ax(4) = gca;
    phase = rad2deg(angle(decon));
    figure(); imagesc(tdown,fshift,fftshift(phase,1)); zoom on; colorbar; title('phase'); axis xy; ylim([0 newFs/2]); ax(10) = gca;
    figure(); plot(fshift,fftshift(mean(phase,2,'omitnan'),1)); zoom on; title('average phase');
    figure(); plot(fshift,fftshift(msc(:,3),1)); zoom on; title('example msc of one time window');
    figure(); plot(tdown,sum(10*log10(decon.*conj(decon))),'.'); zoom on; title('example time evolution of psd of deconvolution');
    figure(); histogram(sum(10*log10(decon.*conj(decon)))); zoom on; title('histogram of cross-spectrum amps');

    figure(); imagesc(tdown,fshift,fftshift(msc,1)); zoom on; colorbar; clim([0 1]);  axis xy; ylim([0 newFs/2]);  ax(5) = gca;
    figure(); imagesc(tdown,fshift,fftshift(mscSmooth,1)); zoom on; colorbar; clim([0 1]);  axis xy; ylim([0 newFs/2]); ax(6) = gca;
end

decon_dB_sort = decon_dB;
msc_sort = msc; %Smooth;
deconSort = decon;
maxWindows = floor(0.5*ntimewin); %keep a certain percentage, discard a certain percentage
msc_sort_index = NaN(nfft/2,maxWindows);

nboot = ntimewin;
deconScrambled = NaN(nfft/2,nboot);

for i = 1:nfft/2 %<-- loop over all frequencies...
    decon_dB_sort(i,:) = sort(decon_dB(i,:));
    [msc_sort_,msc_sort_index_] = sort(10*log10(abs(Sxy(i,:))));
    %[msc_sort_,msc_sort_index_] = sort(10*log10(Syy(i,:)));
    %[msc_sort_,msc_sort_index_] = sort(msc(i,:),"descend");
    %[msc_sort_,msc_sort_index_] = sort(msc(i,:),"ascend");
    %[msc_sort_,msc_sort_index_] = sort(decon_dB(i,:),"descend");
    %[msc_sort_,msc_sort_index_] = sort(decon_dB(i,:),"ascend");
    msc_sort(i,:) = msc_sort_;
    msc_sort_index_ = msc_sort_index_(1:maxWindows);
    msc_sort_index(i,:) = msc_sort_index_;
    decon_ = decon(i,:);
    decon_ = decon_(msc_sort_index_);
    booti = randi(maxWindows,nboot,1);
    deconScrambled(i,:) = decon_(booti);
end

if plotFlag
    figure(); imagesc(tdown,fshift,fftshift(decon_dB_sort,1)); zoom on; colorbar; title('decon psd'); axis xy; ax(7) = gca;
    figure(); imagesc(tdown,fshift,fftshift(msc_sort,1)); zoom on; colorbar; title('msc'); axis xy; ax(8) = gca;
    linkaxes(ax,'xy');
end

deconSort = deconSort(1:nfft/2,1:maxWindows);
msc_sort = msc_sort(1:nfft/2,1:maxWindows);
toc;

deconScrambled = [deconScrambled; flipud(conj(deconScrambled))];
decon2 = fftshift(ifft(deconScrambled,[],1,'symmetric'),1); % now in time domain
d2 = (decon2);

if plotFlag
    figure();
    m = (size(d2,1)+1)/2;
    lags = (-m+1:m-1)';
    plot(lags,normalizeWaveforms(pws(normalizeWaveforms(zpkFilter(d2,lfc,hfc,newFs,npoles,zeroPhaseFlag)),true,true,nu)),'linewidth',2); zoom on; grid on;
    %stackOfStacks = plot_family(zpkFilter(normalizeWaveforms(detrend(d2)),1/20,1/10,newFs,2,true),1:size(d2,2),50,newFs); title('new');

    decon2 = (fftshift(ifft(decon,[],1,'symmetric'),1)); % now in time domain
    d2 = (decon2);
    figure(14); hold on;
    plot(lags,normalizeWaveforms(pws(normalizeWaveforms(zpkFilter(d2,lfc,hfc,newFs,npoles,zeroPhaseFlag)),true,true,nu)),'linewidth',2); zoom on; grid on;
    legend('scrambled decon','regular WL decon');
end
decon2 = fftshift(ifft(deconScrambled,[],1,'symmetric'),1); % now in time domain
d2 = (decon2);
pws_stack = normalizeWaveforms(pws(normalizeWaveforms(zpkFilter(d2,lfc,hfc,newFs,npoles,zeroPhaseFlag)),true,true,nu));
end

%d2 = decon2; stackOfStacks = plot_family(zpkFilter(normalizeWaveforms(detrend(d2)),1/20,1/10,newFs,2,true),1:size(d2,2),50,newFs); title('new');

%%
% figure(); semilogy(sum(mscSmooth),'.'); zoom on; %colorbar; clim([0 1])
% mscSmooth = normalizeWaveforms(mscSmooth);
% figure(); imagesc(mscSmooth); zoom on; colorbar; %clim([0 1])
% figure(); imagesc(1./mscSmooth); zoom on; colorbar; %clim([0 1])
% figure(); imagesc(log10(1./mscSmooth)); zoom on; colorbar; %clim([0 1])
% mscSmooth2 = wiener2(1./msc,[5,5]);
% mscSmooth2 = normalizeWaveforms(mscSmooth2);
% figure(); imagesc(log10(mscSmooth2)); zoom on; colorbar; %clim([0 1])
% figure(); plot(mean(mscSmooth,2),'.'); zoom on;