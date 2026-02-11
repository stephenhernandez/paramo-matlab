function S = padWaveforms(S)
%
% S = padWaveforms(S)
%
% padWaveforms return structure with synced waveform
%
% Syncs traces so that they all have same begin and end
% Some traces will be longer than they originally were (zero-padded)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
sizeS = size(S);
S = S(:);
%lS = length(S);
if isscalar(S) %length(lS) == 1
    return
end

%%
ref = pull(S,'ref');
goodI = ~isnat(ref);

%%
b = pull(S,'b');
e = pull(S,'e');

startTime = min(ref+b);     % will potentially lengthen traces by zero-padding
endTime = max(ref+e);       % will potentially lengthen traces by zero-padding

%%
S_ = cutWaveforms(S(goodI),startTime,0,endTime-startTime,false,false);
S(goodI) = S_;

%%
S = reshape(S,sizeS);