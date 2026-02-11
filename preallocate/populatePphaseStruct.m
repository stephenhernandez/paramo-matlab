function P = populatePphaseStruct(n)

%%
if nargin < 1
    n = 1;
end

%%
P.stnm = "";
P.ntwk = "";
P.chan = "";
P.dist = NaN;
P.azimuth = NaN;
P.t = NaT;
P.res = NaN;
P.wt = NaN;
P.polarity = NaN;

%%
if n > 1
    P = repmat(P,n,1);
end
