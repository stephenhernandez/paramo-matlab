clear; close all; clc;

cd ~/research/now/reventador/
[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v6',10); %,10);

%%
minKurt = min(kurt(:,1:2),[],2,'omitnan');
maxKurt = max(kurt(:,1:2),[],2,'omitnan');
kurtRatio = maxKurt./minKurt;

%%
z2p = z2p(:,1); %max(z2p,[],2);
kurt = kurt(:,1); %max(kurt,[],2);

%%
minAmp = 400;
winlen = 50;
zeroPhaseFlag = false;
goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.04) & ~(Neff == 1) & kurt < 25 & z2p < 1e4;

%%
[N,edges] = histcounts(tabs(goodI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
N = N';
edges = edges(1:end-1)';
figure();
plot(edges,N,'.');
zoom on; grid on;

%%
figure();
semilogy(tabs(goodI),z2p(goodI),'.'); zoom on; grid on;
%hold on;
%semilogy(tabs(~goodI),z2p(~goodI),'.'); zoom on; grid on;

%%
figure();
plot(tabs(goodI),1:sum(goodI),'.'); zoom on; grid on;

tt = tabs(goodI);
close all;
clearvars -except tt

figure();
plot(tt,1:length(tt),'.');
zoom on;

k = 2e3;
n = length(tt);

yI = randsample(n,k);
yI = sort(yI);

%%
SZ = extractWaveforms(tt(yI)-seconds(30),seconds(90),"CASC","HHZ","EC","",true,true);
SN = extractWaveforms(tt(yI)-seconds(30),seconds(90),"CASC","HHN","EC","",true,true);
SE = extractWaveforms(tt(yI)-seconds(30),seconds(90),"CASC","HHE","EC","",true,true);

%%
refsZ = pull(SZ,'ref');
refsN = pull(SN,'ref');
refsE = pull(SE,'ref');

goodI = ~(isnat(refsZ) | isnat(refsN) | isnat(refsE));

SZ = interpolateWaveforms(SZ);
SN = interpolateWaveforms(SN);
SE = interpolateWaveforms(SE);

dz = pull(SZ);
dn = pull(SN);
de = pull(SE);

dzi = sum(~isfinite(dz))';
dni = sum(~isfinite(dn))';
dei = sum(~isfinite(de))';

goodI = goodI & ~(dzi | dni | dei);

%%
SZc = SZ(goodI);
SNc = SN(goodI);
SEc = SE(goodI);

refsZc = pull(SZc,'ref');
refsNc = pull(SNc,'ref');
refsEc = pull(SEc,'ref');

%%
SZ = extractWaveforms(tt(yI)-seconds(30),seconds(90),"ANTS","HHZ","EC","",true,true);
SN = extractWaveforms(tt(yI)-seconds(30),seconds(90),"ANTS","HHN","EC","",true,true);
SE = extractWaveforms(tt(yI)-seconds(30),seconds(90),"ANTS","HHE","EC","",true,true);

%
refsZ = pull(SZ,'ref');
refsN = pull(SN,'ref');
refsE = pull(SE,'ref');

goodI = ~(isnat(refsZ) | isnat(refsN) | isnat(refsE));

SZ = interpolateWaveforms(SZ);
SN = interpolateWaveforms(SN);
SE = interpolateWaveforms(SE);

dz = pull(SZ);
dn = pull(SN);
de = pull(SE);

dzi = sum(~isfinite(dz))';
dni = sum(~isfinite(dn))';
dei = sum(~isfinite(de))';

goodI = goodI & ~(dzi | dni | dei);


SZa = SZ(goodI);
SNa = SN(goodI);
SEa = SE(goodI);

refsZa = pull(SZ,'ref');
refsNa = pull(SN,'ref');
refsEa = pull(SE,'ref');

clear SZ SN SE

%%
tic;
tw = 0.02;
newFs = 100;
d = double(pull(taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SNa),newFs))),tw)));
toc;

secDur = 45;
winlen = secDur*newFs;

dnoise = detrend(d(1:winlen,:));
d = detrend(d(winlen+1:2*winlen,:));
toc;

extraPadding = 1;
nfft = 2^(nextpow2(winlen) + extraPadding);
toc;

[pxxN_,~] = pmtm(dnoise,{2,'trace'},nfft,newFs);
toc;

clear dnoise;
[pxxS_,fxx] = pmtm(d,{2,'trace'},nfft,newFs);
toc;

clear d
pxxS_ = pxxS_./pxxN_;
toc;

snr = pxxS_;
clear pxxS_;
toc;

close all;
figure(); semilogx(fxx(2:end),nanmedian(snr(2:end,:),2)); zoom on; 
%after experimenting with these curves, i landed on a filter of 0.3 - 2.4

%%
close all;
tw = 0.02;
lfc = 0.4;
hfc = 1.6;

plotWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(SZc(1:15)),tw),lfc,hfc),[],[],[],[],1); xlim([40 90]);
plotWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(SNc(1:15)),tw),lfc,hfc),[],[],[],[],1); xlim([40 90]);
plotWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(SEc(1:15)),tw),lfc,hfc),[],[],[],[],1); xlim([40 90]);

plotWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(SZa(1:15)),tw),lfc,hfc),[],[],[],[],1); xlim([40 90]);
plotWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(SNa(1:15)),tw),lfc,hfc),[],[],[],[],1); xlim([40 90]);
plotWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(SEa(1:15)),tw),lfc,hfc),[],[],[],[],1); xlim([40 90]);

%%
% i want to concatenate all CASC filtered data 
% then i want to align everything, and save those shifts and apply them to
% new time vector. with that new time vector, i will reload X number of
% waveforms with confidence they will be aligned at CAYR, CASC, and ANTS

close all; clc;
tw = 0.04;
lfc = 0.4;
hfc = 1.6;

tic;
clear fxx snr S
S = taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SZc),newFs))),tw);
S = filterWaveforms(taperWaveforms(detrendWaveforms(S),tw),lfc,hfc);
d = double(pull(S));
dz = detrend(d(winlen+1:2*winlen,:));
toc;

S = taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SNc),newFs))),tw);
S = filterWaveforms(taperWaveforms(detrendWaveforms(S),tw),lfc,hfc);
d = double(pull(S));
dn = detrend(d(winlen+1:2*winlen,:));
toc;

S = taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SEc),newFs))),tw);
S = filterWaveforms(taperWaveforms(detrendWaveforms(S),tw),lfc,hfc);
d = double(pull(S));
de = detrend(d(winlen+1:2*winlen,:));
toc;

clear d S;
%d = normalizeWaveforms([normalizeWaveforms(dz); normalizeWaveforms(dn); normalizeWaveforms(de)]);
d = [dz; dn; de];
amps = rms(d)';
p2rms = peak2rms(d)';
d = normalizeWaveforms([normalizeWaveforms(dz); normalizeWaveforms(dn); normalizeWaveforms(de)]);
toc;


%%
SZa2 = extractWaveforms(t4,seconds(45),"ANTS","HHZ","EC","",true,true);
dza_good = ~isnat(pull(SZa2,'ref')) & sum(isfinite(double(pull(SZa2))))' >= 4500;
close all
figure(); plot(dza_good,'o'); zoom on;
dza = double(pull(filterWaveforms(detrendWaveforms(SZa2(dza_good)),0.4,1.6)));

%%
SNa2 = extractWaveforms(t5,seconds(45),"ANTS","HHN","EC","",true,true);
SEa2 = extractWaveforms(t5,seconds(45),"ANTS","HHE","EC","",true,true);

%%
clear SZa SNa SEa
SZc = extractWaveforms(t5,seconds(45),"CASC","HHZ","EC","",true,true);
SNc = extractWaveforms(t5,seconds(45),"CASC","HHN","EC","",true,true);
SEc = extractWaveforms(t5,seconds(45),"CASC","HHE","EC","",true,true);

%%
SZcay = extractWaveforms(t5,seconds(45),"CAYR","SHZ","EC","",true,true);
SNcay = extractWaveforms(t5,seconds(45),"CAYR","SHN","EC","",true,true);
SEcay = extractWaveforms(t5,seconds(45),"CAYR","SHE","EC","",true,true);















