function S = populatePhaseStruct(sizeS)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Saturday, Jul 27, 2019

%%
if nargin < 1
    sizeS = 1;
end

%%
S.stnm = "";
S.ntwk = "";
S.chan = "";
S.locid = "";

S.dist = NaN;       % distance in km
S.azimuth = NaN;    % back-azimuth
S.t = NaT;          % arrival time
S.res = NaN;        % residual
S.wt = NaN;         % weight (quality)
S.polarity = NaN;   % polarity
S.tt = NaN;         % travel-time (in milliseconds)
S.toa = NaN;        % take-off angle

%%
if ~isscalar(sizeS)
    S = repmat(S,sizeS);
    return;
end
S = repmat(S,sizeS,1); %make column vector