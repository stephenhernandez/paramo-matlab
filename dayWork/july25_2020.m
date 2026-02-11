clear; close all; clc;
cd ~/research/now/pichincha_nxcorr/
load ggp_2008_2020_PLD.mat

%%
Ntot = size(tabs(:),1);
trainFrac = 0.5;
trainDatasetNum = floor(trainFrac*Ntot);
npoles = 4;
nMinutes = 5;

%%
I = (1:2:Ntot)'; %unique(randsample(Ntot,trainDatasetNum)); %the unique command is not necessary
tTrain = tabs(I);
S = extractWaveforms(tTrain - seconds(nMinutes*60),seconds(nMinutes*60 + 10),"PINO","SHZ","EC","",false,true);

%%
refs = pull(S,'ref');
rI = isnat(refs);
S(rI) = [];
tabs(rI) = [];

%%
newFs = 40;
S = detrendWaveforms(S);
S = taperWaveforms(S,0.002);
S = filterWaveforms(S,1/64);
S = resampleWaveforms(S,newFs);
Sorig = S;

%%
S = Sorig;
close all;
clearvars -except S tabs nMinutes npoles newFs I Sorig
clipFlag = false;
if clipFlag
    S = clipWavforms(S,2,true);
end
S = pull(S);
Fs = newFs; %sampling rate
maxAmpRMS = rssq(detrend(S(end-10*Fs:end)))';

%
windowSize = 10; %seconds
envMean = [];
envStd = envMean;
%stackMed = stackMean;
stdMean = envMean;
stdStd = envMean;
skewMean = envMean;
skewStd = envMean;
kurtMean = envMean;
kurtStd = envMean;

%
S = S(1:nMinutes*60*Fs + 075,:);
S = detrend(S);
boxSamples = windowSize*Fs;
box = ones(boxSamples,1)/boxSamples;

%
n = 0;
close all;
lfc = [1/8; 1/4]; %; 1/2; 1; 2; 4; 8; 16];
%lfc = [1/16; 1/8; 1/4; 1/2; 1; 2; 4; 8; 16; 32];

%%
l_lfc = length(lfc);
ll = l_lfc - 1;
figure('units','normalized','outerposition',[0 0 1 1]);
tic;
for i = 1:ll
    lfc_ = lfc(i);
    for j = i+1:l_lfc
        disp([i,j]);
        n = n+1;
        
        %%
        hfc_ = lfc(j);
        Sf = zpkFilter(S,lfc_,hfc_,Fs,npoles);
        %Sf = normalizeWaveforms(Sf);
        Sf = fftfilt(box,Sf.^2); %var
        std_ = sqrt(abs(Sf));
        std_ = std_(3*windowSize*Fs:end,:);
        stdMean = [stdMean nanmedian(std_./nanmedian(std_),2)];
        stdStd = [stdStd mad(std_./nanmedian(std_),1,2)];
        toc;
        
        %%
        Sf = zpkFilter(S,lfc_,hfc_,Fs,npoles);
        Sf = fftfilt(box,Sf.^3);
        Sf = Sf(3*windowSize*Fs:end,:)./(std_.^3);
        %Sf = normalizeWaveforms(Sf);
        skewMean = [skewMean nanmedian(Sf,2)];
        skewStd = [skewStd mad(Sf,1,2)];
        toc;
        
        %%
        Sf = zpkFilter(S,lfc_,hfc_,Fs,npoles);
        Sf = fftfilt(box,Sf.^4);
        Sf = Sf(3*windowSize*Fs:end,:)./(std_.^4);
        Sf = Sf./nanmedian(Sf);
        %Sf = normalizeWaveforms(Sf);
        kurtMean = [kurtMean nanmedian(Sf,2)];
        kurtStd = [kurtStd mad(Sf,1,2)];
        toc;
        
        %%
        Sf = zpkFilter(S,lfc_,hfc_,Fs,npoles);
        Sf = abs(hilbert(Sf));
        Sf = zpkFilter(Sf,-inf,1/windowSize,Fs,1);
        Sf = Sf(3*windowSize*Fs:end,:);
        Sf = Sf./nanmedian(Sf);
        %Sf = normalizeWaveforms(Sf);
        envMean = [envMean nanmedian(Sf,2)];
        envStd = [envStd mad(Sf,1,2)];
        %stackMed = [stackMed nanmedian(Sf,2)];
        toc;
        
        %%
        spn = ll*(i-1) + (j-1);
        ax_(n) = subplot(ll,ll,spn);
        tdum = (0:size(envMean,1)-1)'/Fs - size(envMean,1)/Fs;
        p1 = plot(tdum,envMean(:,n),'linewidth',3);
        p1.Color(4) = 0.7;
        %hold on; p2 = plot(tdum,stackMed(:,n),'linewidth',3); p2.Color(4) = 0.7;
        
        %%
        if lfc_ < 1
            titleStr = [num2str(1/lfc_),' sec. - '];
        else
            titleStr = [num2str(lfc_),' hz. - '];
        end
        if lfc(j) < 1
            titleStr = [titleStr,num2str(1/lfc(j)),' sec.'];
        else
            titleStr = [titleStr,num2str(lfc(j)),' Hz.'];
        end
        title(titleStr);
        % ax = gca;
        ax_(n).YScale = 'log';
        xlabel('time to failure [sec.]'); zoom on;
        toc;
    end
end
linkaxes(ax_,'x');
toc;

%%
figure('units','normalized','outerposition',[0 0 0.5 1]);
meanax(1) = subplot(311);
plot(tdum,envMean,'linewidth',3);
meanax(1).YScale = 'log';
grid on; zoom on; title('Env. Mean');
meanax(2) = subplot(312);
plot(tdum,envStd,'linewidth',3);
grid on; zoom on; title('Env. Std');
meanax(2).YScale = 'log';
meanax(3) = subplot(313);
plot(tdum,envMean./envStd,'linewidth',3);
grid on; zoom on; title('Env. SNR');
linkaxes(meanax,'x');

%
figure('units','normalized','outerposition',[0 0 0.5 1]);
stdax(1) = subplot(311);
plot(tdum,stdMean,'linewidth',3);
grid on;
zoom on; title('Std Mean');
stdax(1).YScale = 'log';
stdax(2) = subplot(312);
plot(tdum,stdStd,'linewidth',3);
grid on; zoom on; title('Std Std');
stdax(2).YScale = 'log';
stdax(3) = subplot(313);
plot(tdum,stdMean./stdStd,'linewidth',3);
grid on; zoom on; title('Std SNR');
linkaxes(stdax,'x');

%
figure('units','normalized','outerposition',[0 0 0.5 1]);
skewax(1) = subplot(311);
plot(tdum,skewMean,'linewidth',3); grid on; zoom on; title('Skewness Mean');
skewax(2) = subplot(312);
plot(tdum,skewStd,'linewidth',3);
grid on; zoom on; title('Skewness Std');
skewax(2).YScale = 'log';
skewax(3) = subplot(313);
plot(tdum,skewMean./skewStd,'linewidth',3);
grid on; zoom on; title('Skewness SNR');
linkaxes(skewax,'x');

%
figure('units','normalized','outerposition',[0 0 0.5 1]);
kurtax(1) = subplot(311);
plot(tdum,kurtMean,'linewidth',3); grid on; zoom on; title('Kurtosis Mean');
kurtax(1).YScale = 'log';
kurtax(2) = subplot(312);
plot(tdum,kurtStd,'linewidth',3);
grid on; zoom on; title('Kurtosis Std');
kurtax(2).YScale = 'log';
kurtax(3) = subplot(313);
plot(tdum,kurtMean./kurtStd,'linewidth',3);
grid on; zoom on; title('Kurtosis SNR');
linkaxes(kurtax,'x');
