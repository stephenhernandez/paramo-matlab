function S = noisePreprocess(S,newFs,lfc,hfc,npoles,tw,zeroPhaseFlag,transferFlag,smoother,expo)
%
% noisePreprocess routine to apply time domain normalization
%
% S = noisePreprocess(S,newFs,lfc,hfc,npoles,tw,zeroPhaseFlag,smoother)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 9; smoother = 1/lfc; end
if nargin < 10; expo = 2; end

%%
%S = padWaveforms(S);
S = differentiateWaveforms(S);
S = detrendWaveforms(S);
S = resampleWaveforms(S,newFs);

%% basic filtering
if transferFlag
    S = transferWaveforms(S,lfc,hfc,npoles,[],'vel',true);
else
    S = filterWaveforms(S,lfc,hfc,npoles,[],zeroPhaseFlag);
end

S = taperWaveforms(S,tw);
S = intWaveforms(S);

%% time-domain normalizing
% S = oneBitNorm(S);
% S = padWaveforms(S);
% S = fdWhiten(S,lfc,hfc,dW,newFs);
S = tdNorm(S,smoother,expo,newFs);
S = padWaveforms(S);
