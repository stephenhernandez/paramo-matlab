function S = signWaveforms(S)
%
% signWaveforms return structure with bit-normalized waveforms
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Oct 25, 2021

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    S_ = S(i);
    if isnat(S_.ref)
        fprintf(1,'\n');
        S(i) = populateWaveforms(1);
        continue;
    end

    %%
    d = double(S_.d);
    d = sign(d);

    S_ = dealHeader(S_,d);
    S(i) = S_;
end
S = reshape(S,sizeS);