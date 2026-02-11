function [d,t_] = nxcorrPreprocess(S,lfc,hfc,newFs,npoles,zeroPhaseFlag,tw,diffFlag,transferFlag)
if nargin < 9
    transferFlag = false;
end

%%
if diffFlag
    S = differentiateWaveforms(S);
end

%%
S = detrendWaveforms(S);

%%
if transferFlag
    S = taperWaveforms(S,tw);
    S = transferWaveforms(S,lfc,hfc,npoles,newFs,'vel');
else
    S = filterWaveforms(S,lfc,hfc,npoles,tw,zeroPhaseFlag);
    S = resampleWaveforms(S,newFs);
end

%%
S = syncWaveforms(S);

%%
d = pull(S);
d = sign(d);

%%
t_ = getTimeVec(S(1));