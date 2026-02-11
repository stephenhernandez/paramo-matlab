clear;
cd '/home/shernandez/research/now/pichincha/pichincha_nxcorr';
load('GGP20082020.mat','tabs','p2p');
tI = tabs >= datetime(2015,01,07) & p2p >= 500;

t2 = tabs(tI); 
S = extractWaveforms(t2-seconds(5),seconds(20),"PINO",["SHZ";"SHN";"SHE"],"EC","",true,true,2,true,[1/4,-inf,false,false]);
sizeS = size(S);
S = S(:);
S = struct2table(S);

clearvars -except S sizeS
save('ggp_templates_25OCT2022');