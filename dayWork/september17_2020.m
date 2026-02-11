clear; close all; clc; 
cd ~/research/now/sangay/

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"PUYO","HHZ","EC","",true,true);
z2pWA_PUYO_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"TAIS","HHZ","EC","",true,true);
z2pWA_TAIS_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"PORT","HHZ","EC","",true,true);
z2pWA_PORT_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"PKYU","HHZ","EC","",true,true);
z2pWA_PKYU_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"TAMH","HHZ","EC","",true,true);
z2pWA_TAMH_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"BMAS","BHZ","EC","",true,true);
z2pWA_BMAS_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"BULB","BHZ","EC","",true,true);
z2pWA_BULB_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"BPAT","BHZ","EC","",true,true);
z2pWA_BPAT_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
save('sangayWoodAndersonAmplitudes');

%%
load('sangayWoodAndersonAmplitudes');
load('SAGA_WoodAnderson.mat','refs');
S = extractWaveforms(refs,seconds(180),"BRUN","BHZ","EC","",true,true);
z2pWA_BRUN_raw = max(abs(double(pull(filterWaveforms(detrendWaveforms(S),0.6,1.2)))))';
clear S;
save('sangayWoodAndersonAmplitudes');
