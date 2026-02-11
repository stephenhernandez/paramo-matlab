clear; close all; clc;
cd ~/research/now/sangay/snr/
load brun_data.mat

S = table2struct(S);
S = resampleWaveforms(intWaveforms(filterWaveforms(taperWaveforms(detrendWaveforms(differentiateWaveforms(S)),0.04),0.6,1.2,6)),8);
S = detrendWaveforms(S);
p2p = 0.5*(pull(S,'depmax') - pull(S,'depmin'));
t = pull(S,'ref');

figure(); semilogy(t,p2p,'.'); zoom on;
figure(); plot(log10(sort(p2p(t>=datetime(2015,01,01)))),1:length(t((t>=datetime(2015,01,01)))),'.'); zoom on;
pI = log10(p2p) >= 1.5 & log10(p2p) < 1.8 & t >= datetime(2015,01,01);
sum(pI)

%%
[maxccp,plags,maxccn,nlags] = doccFreqCircShift(double(pull(resampleWaveforms(S(pI),8))),true);

%%
[newEvents,newFamilies,maxccp,plags] = pruneAndMergeEvents(double(pull(resampleWaveforms(S(pI),8))),maxccp,plags,[0.7 0.9 0.9],'weighted');
ccParams = [0.7 0.9 0.9];
data = double(pull(resampleWaveforms(S(pI),8)));

%%
t = t(pI);
p2p = p2p(pI);
for i = 1:length(newFamilies)
lF(i) = length(newFamilies{i});
end
figure(); loglog(lF); zoom on;
figure(); plot(cumsum(lF),'o'); zoom on;

%%
clearvars -except data maxccp plags ccParams lF newEvents newFamilies p2p t
newFs = 8;
lfc = 0.6;
hfc = 1.2;
npoles = 6;
cd ..
save('BRUN_big_events_2015_2020_forSVD','-v7.3');

%%
[shifted_data,maxccp_,G,plags_,raw_shifts] = apply_vdcc(newEvents(:,1:200),[],false,true);
[U,SVs,V] = svd(shifted_data,0);
close all; figure(); plot(U(:,1:5)); zoom on;

%%
clear; clear all; load BRUN_big_events_2015_2020_forSVD.mat
[shifted_data,maxccp_,G,plags_,raw_shifts] = apply_vdcc(newEvents(:,1:200),[],false,true);
[U,SVs,V] = svd(shifted_data,0);

%%
clearvars -except U newFs lfc hfc kstnm chan ntwk
kstnm = "BRUN";
chan = "BHZ";
ntwk = "EC";
snr = [];
U = U(:,1:100);
save('brun_svd_basis_functions');

%%
clear
load('sangay_svd_basis_functions');
U_ = U;
load('brun_svd_basis_functions','U');
U_(:,11,1:100) = U(1:1440,:);
U = U_;
chan = [chan; "BHZ"];
ntwk = [ntwk; "EC"];
kstnm = [kstnm; "BRUN"];
clear U_;

%%
chan
ntwk
kstnm
%save('sangay_svd_basis_functions');