tw = 0.1;
noiseWin = 5;
secDur = 150;
newFs = 8;
lfc = 0.6;
hfc = 1.2;

S = extractWaveforms(tt(tI)-seconds(noiseWin),seconds(secDur),"PORT","HHZ");
S = differentiateWaveforms(S);
S = detrendWaveforms(S);
S = taperWaveforms(S,tw);
S = filterWaveforms(S,lfc,hfc);
S = taperWaveforms(S,tw);
S = intWaveforms(S);
S = taperWaveforms(S,tw);
S = resampleWaveforms(S,newFs);

%%
S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BBIL","BHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BBIL = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BMAS","BHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BMAS = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BPAT","BHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BPAT = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BRTU","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BRTU = double(pull(S));

%%
tw = 0.1; 
noiseWin = 5; 
secDur = 150; 
newFs = 8; 
lfc = 0.6;
hfc = 1.2;

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BBIL","BHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BBIL = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BMAS","BHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BMAS = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BPAT","BHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BPAT = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BRTU","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
BRTU = double(pull(S));

% 
% S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BRUN","BHZ");
% S = resampleWaveforms(S,newFs);
% S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
% BRUN = double(pull(S));
% 
% S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"BULB","BHZ");
% S = resampleWaveforms(S,newFs);
% S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
% BULB = double(pull(S));
% 
% S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"CHSH","HHZ");
% S = resampleWaveforms(S,newFs);
% S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
% CHSH = double(pull(S));
% 
% S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"COHC","HHZ");
% S = resampleWaveforms(S,newFs);
% S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
% COHC = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"JSCH","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
JSCH = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"PIAT","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
PIAT = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"PKYU","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
PKYU = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"PORT","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
PORT = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"PUYO","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
PUYO = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"TAIS","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
TAIS = double(pull(S));

S = extractWaveforms(t(tI)-seconds(noiseWin),seconds(secDur),"TAMH","HHZ");
S = resampleWaveforms(S,newFs);
S = intWaveforms(detrendWaveforms(taperWaveforms(filterWaveforms(differentiateWaveforms(S),lfc,hfc),tw)));
TAMH = double(pull(S));
