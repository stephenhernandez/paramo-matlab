function S = intWaveforms(S,ncum,verboseFlag)
%
% intWaveforms return structure with integrated waveforms
%
% S = intWaveforms(S,ncum)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    ncum = 1;
end

if nargin < 3
    verboseFlag= false;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);

if verboseFlag
    if ncum == 1
        disp(['Integrating ',num2str(lS),' trace(s) in time domain once']);
    else
        disp(['Integrating ',num2str(lS),' trace(s) in time domain ',num2str(ncum),' times']);
    end
end

%%
for i = 1:lS
    S_ = S(i);
    if isnat(S_.ref)
        continue;
    end
    
    %%
    dt = S_.delta;
    d = S_.d;
    
    goodI = isfinite(d);
    n = 1;
    while n <= ncum
        d(goodI) = cumsum(d(goodI),"omitnan")*dt;
        n = n+1;
    end
    S_ = dealHeader(S_,d);
    S(i) = S_;
end
S = reshape(S,sizeS);
