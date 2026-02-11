function M = populateMagnitudeStruct(n)

%%
if nargin < 1
    n = 1;
end

%%
M.stnm = "";
M.ntwk = "";
M.chan = "";
M.type = "";
M.dist = NaN;
M.azimuth = NaN;
M.value = NaN;
M.res = NaN;
M.amp = NaN;

%%
if n > 1
    M = repmat(M,n,1);
end