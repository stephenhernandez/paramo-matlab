clear; close all;

%%
close all; 
cd ~/research/now/sangay/
load sangaySubspaceDetectorSAGA_v2.mat
z2p(tabs <= datetime(2018,11,11)) = 4*z2p(tabs <= datetime(2018,11,11));
gI = ~(z2p <= 2e3 | p2rms >= 10 | kurt < 10 | NCC < 0.1 | z2p > 1e6);

figure(); 
semilogy(tabs(gI),t2r(tabs(gI),hours(24)),'.'); 
zoom on;
grid on;

tt = tabs(gI);
k = 6e3;
n = sum(gI);

yI = randsample(n,k);
yI = sort(yI);

figure(); plot(sort(tt(yI)),1:length(yI),'.'); zoom on; grid on;

%%
kstnm = "SAGA";
totDur = 60;
noiseWin = -30;
SZ = extractWaveforms(tt(yI)-seconds(noiseWin),seconds(totDur),kstnm,"HHZ","EC","",true,true);
SN = extractWaveforms(tt(yI)-seconds(noiseWin),seconds(totDur),kstnm,"HHN","EC","",true,true);
SE = extractWaveforms(tt(yI)-seconds(noiseWin),seconds(totDur),kstnm,"HHE","EC","",true,true);

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

%
SZs = SZ(goodI);
SNs = SN(goodI);
SEs = SE(goodI);

refsZc = pull(SZs,'ref');
refsNc = pull(SNs,'ref');
refsEc = pull(SEs,'ref');

%%
close all; clc;
tw = 0.04;
lfc = 0.5;
hfc = 4.0;
newFs = 20;

tic;
clear fxx snr S
S = taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SZs),newFs))),tw);
S = filterWaveforms(taperWaveforms(detrendWaveforms(S),tw),lfc,hfc);
dz = double(pull(S));
toc;

S = taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SNs),newFs))),tw);
S = filterWaveforms(taperWaveforms(detrendWaveforms(S),tw),lfc,hfc);
dn = double(pull(S));
toc;

S = taperWaveforms(detrendWaveforms(interpolateWaveforms(resampleWaveforms(detrendWaveforms(SEs),newFs))),tw);
S = filterWaveforms(taperWaveforms(detrendWaveforms(S),tw),lfc,hfc);
de = double(pull(S));
toc;

clear S;
d = [dz; dn; de];
amps = rms(d)';

%%
amps(refsZc <= datetime(2018,11,11)) = 4*amps(refsZc <= datetime(2018,11,11));
p2rms = peak2rms(d)';
d = normalizeWaveforms([normalizeWaveforms(dz); normalizeWaveforms(dn); normalizeWaveforms(de)]);
toc;

%%
kstnm = "SAGA";
totDur = 60;
noiseWin = 0;
SZ = extractWaveforms(t4-seconds(noiseWin),seconds(totDur),kstnm,"HHZ","EC","",true,true);
SN = extractWaveforms(t4-seconds(noiseWin),seconds(totDur),kstnm,"HHN","EC","",true,true);
SE = extractWaveforms(t4-seconds(noiseWin),seconds(totDur),kstnm,"HHE","EC","",true,true);

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

%
SZs = SZ(goodI);
SNs = SN(goodI);
SEs = SE(goodI);

refsZc = pull(SZs,'ref');
refsNc = pull(SNs,'ref');
refsEc = pull(SEs,'ref');
