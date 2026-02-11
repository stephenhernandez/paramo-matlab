clear; close all; clc;

cd ~/research/now/sangay/;

% lfc = filterObject(1);
% hfc = filterObject(2);
% waFlag = filterObject(3);
% S_ = transferWaveforms(S_,lfc,hfc,6,newFs,'disp',true,waFlag);
                    
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"PUYO","HHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_PUYO_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"TAIS","HHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_TAIS_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"PORT","HHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_PORT_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"PKYU","HHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_PKYU_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"TAMH","HHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_TAMH_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','refs');
S = extractWaveforms(t6,seconds(180),"BMAS","BHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_BMAS_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"BULB","BHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_BULB_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"BPAT","BHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_BPAT_raw = S;
clear S;
save('sangayDisplacementAmplitudes');

%%
load('sangayDisplacementTimes.mat','t6');
S = extractWaveforms(t6,seconds(180),"BRUN","BHZ","EC","",true,true,1);
S = resampleWaveforms(S,10);
z2pWA_BRUN_raw = S;
clear S;
save('sangayDisplacementAmplitudes');
