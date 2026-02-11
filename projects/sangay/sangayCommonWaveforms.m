
%jsch_data.mat  port_data.mat  tais_data.mat
%brtu_data.mat  chsh_data.mat  piat_data.mat  puyo_data.mat  tamh_data.mat
%brun_data.mat  cohc_data.mat  pkyu_data.mat  saga_data.mat  tCommon.mat
clear; close all; clc;

% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"PUYO","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('puyo_data');

%%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"SAGA","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('saga_data');
% 
% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BBIL","BHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('bbil_data');
% 
% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BMAS","BHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('bmas_data');
% 
% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BREF","BHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('bref_data');
% 
% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BULB","BHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('bulb_data');
% 
% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BPAT","BHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('bpat_data');
% 
% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BRUN","BHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('brun_data');

% %%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"TAIS","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('tais_data');
% 
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"PKYU","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('pkyu_data');
% 
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"BRTU","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('brtu_data');
% 
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"TAMH","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('tamh_data');
% 
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"CHSH","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('chsh_data');

%%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"PORT","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('port_data');
% 
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"PIAT","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('piat_data');
% 
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"JSCH","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('jsch_data');

%%
% clear; close all; clc;
% cd ~/research/now/sangay/snr/
% load tCommon.mat
% S = extractWaveforms(tCommon - seconds(30),seconds(180),"COHC","HHZ","EC","",true,true);
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% refs(rI) = [];
% refs = refs + seconds(30);
% S = resampleWaveforms(S,100);
% S = struct2table(S);
% save('cohc_data');

%%
clear; close all; clc;
cd ~/research/now/sangay/snr/
load tCommon.mat
S = extractWaveforms(tCommon - seconds(30),seconds(180),"BOSC","HHZ","EC","",true,true);
refs = pull(S,'ref');
rI = isnat(refs);
S(rI) = [];
refs(rI) = [];
refs = refs + seconds(30);
S = resampleWaveforms(S,100);
S = struct2table(S);
save('bosc_data');

clear; close all; clc;
cd ~/research/now/sangay/snr/
load tCommon.mat
S = extractWaveforms(tCommon - seconds(30),seconds(180),"PIS1","HHZ","EC","",true,true);
refs = pull(S,'ref');
rI = isnat(refs);
S(rI) = [];
refs(rI) = [];
refs = refs + seconds(30);
S = resampleWaveforms(S,100);
S = struct2table(S);
save('pis1_data');