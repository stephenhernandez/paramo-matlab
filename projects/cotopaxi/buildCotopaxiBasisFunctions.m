clear; 
cd ~/research/now/cotopaxi/
load cotopaxiSubspaceDetectorBREF_v7.mat

rmsEstimate = nanmean(z2p./p2rms,2);
zI = nanmean(z2p,2) >= 500 & Neff > 2 & nanmean(z2p,2) < 1e5 & nanmean(p2rms,2) <= 3.5 & nanmean(kurt,2) <= 6.5;
close all;
figure(); semilogy(tabs(zI),nanmean(z2p(zI,:),2),'.'); zoom on;
figure(); plot(tabs(zI),1:length(tabs(zI)),'.'); zoom on;
[N,edges] = histcounts(tabs(zI),dateshift(min(tabs(zI)),'start','day'):days(1):dateshift(max(tabs(zI)),'end','day')+1);
N = N';
edges = edges(1:end-1)';
figure(); stairs(edges,N,'linewidth',2); zoom on; grid on;
sum(zI)

newFs = 20;
tw = 0.1;
lfc = 12/16;
hfc = 3;
npoles = 4;

%%
S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BREF","BHZ","EC","",true,true);
BREFZ = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BREF","BHN","EC","",true,true);
BREFN = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BREF","BHE","EC","",true,true);
BREFE = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BTAM","BHZ","EC","",true,true);
BTAMZ = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BTAM","BHN","EC","",true,true);
BTAMN = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BTAM","BHE","EC","",true,true);
BTAME = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BMOR","BHZ","EC","",true,true);
BMORZ = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BMOR","BHN","EC","",true,true);
BMORN = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BMOR","BHE","EC","",true,true);
BMORE = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BVC2","BHZ","EC","",true,true);
BVC2Z = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BVC2","BHN","EC","",true,true);
BVC2N = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BVC2","BHE","EC","",true,true);
BVC2E = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BNAS","BHZ","EC","",true,true);
BNASZ = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BNAS","BHN","EC","",true,true);
BNASN = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BNAS","BHE","EC","",true,true);
BNASE = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"VC1","SHZ","EC","",true,true);
VC1 = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"NAS2","SHZ","EC","",true,true);
NAS2 = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"TAMB","SHZ","EC","",true,true);
TAMB = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BVC2","HHZ","EC","",true,true);
BVC2HZ = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BVC2","HHN","EC","",true,true);
BVC2HN = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BVC2","HHE","EC","",true,true);
BVC2HE = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BNAS","HHZ","EC","",true,true);
BNASHZ = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BNAS","HHN","EC","",true,true);
BNASHN = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

S = extractWaveforms(tabs(zI)-seconds(5),seconds(35),"BNAS","HHE","EC","",true,true);
BNASHE = resampleWaveforms(taperWaveforms(detrendWaveforms(interpolateWaveforms(S)),tw),newFs);

%%
tic;
BREFZ = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BREFZ)),0.75,3));
BREFN = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BREFN)),0.75,3));
BREFE = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BREFE)),0.75,3));
BTAMZ = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BTAMZ)),0.75,3));
BTAMN = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BTAMN)),0.75,3));
BTAME = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BTAME)),0.75,3));
BMORZ = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BMORZ)),0.75,3));
BMORN = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BMORN)),0.75,3));
BMORE = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BMORE)),0.75,3));
BVC2Z = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BVC2Z)),0.75,3));
BVC2N = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BVC2N)),0.75,3));
BVC2E = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BVC2E)),0.75,3));
BNASZ = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BNASZ)),0.75,3));
BNASN = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BNASN)),0.75,3));
BNASE = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BNASE)),0.75,3));

VC1 = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(VC1)),0.75,3));
NAS2 = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(NAS2)),0.75,3));
TAMB = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(TAMB)),0.75,3));

BVC2HZ = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BVC2HZ)),0.75,3));
BVC2HN = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BVC2HN)),0.75,3));
BVC2HE = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BVC2HE)),0.75,3));
BNASHZ = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BNASHZ)),0.75,3));
BNASHN = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BNASHN)),0.75,3));
BNASHE = interpolateWaveforms(filterWaveforms(interpolateWaveforms(detrendWaveforms(BNASHE)),0.75,3));
toc;

%%
tic;
maxBasisFunctions = 50;
maxN = 700; % 35 seconds @ 20 Hz == 700 points
S = [BREFZ BREFN BREFE...
    BTAMZ BTAMN BTAME...
    BMORZ BMORN BMORE...
    BVC2Z BVC2N BVC2E...
    BNASZ BNASN BNASE...
    VC1 NAS2 TAMB...
    BVC2HZ BVC2HN BVC2HE...
    BNASHZ BNASHN BNASHE];

lS = size(S,2);
Uall = [];
for i = 1:lS
    S_ = S(:,i);
    d = double(pull(S_));
    refs = pull(S_,'ref');
    badI = isnat(refs) | sum(~isfinite(d))' ~= 0 | rms(d)' == 0;
    d(:,badI) = [];
    disp(['number of good columns: ',num2str(size(d,2))]);
    if size(d,1) > maxN
        d = d(1:maxN,:);
    end
    refs(badI) = [];
    S_(badI) = [];
    disp(strcat(S_(1).kstnm,S_(1).kcmpnm));
    disp(refs(1));
    d(~isfinite(d)) = 0;
    d = taper(detrend(d),tw);
    d = normalizeWaveforms(d);
    [U_,~,~] = svd(d);
    U_ = U_(:,1:maxBasisFunctions);
    U_ = reshape(U_,[maxN,1,maxBasisFunctions]);
    Uall = cat(2,Uall,U_);
end
toc;


