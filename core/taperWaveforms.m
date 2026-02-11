function S = taperWaveforms(S,tw)
%
% taperWaveforms return structure with tapered waveforms
%
% S = taperWaveforms(S,tw,N)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019
% Updated: 19 AUG 2022

%%
if nargin < 2
    tw = 0.001;
end

if isempty(tw)
    return;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    S_ = S(i);
    d = S_.d;

    %%
    d = taper(d,tw);
    S_ = dealHeader(S_,d);
    S(i) = S_;
end
S = reshape(S,sizeS);
