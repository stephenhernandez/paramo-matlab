function E = populateSeisCompStructure(N)

if nargin < 1
    N = 1; 
end

%%
E.t = NaN;
E.lat = NaN;
E.lon = NaN;
E.depth = NaN;
E.mag = NaN;
E.timerr = NaN;
E.laterr = NaN;
E.lonerr = NaN;
E.deptherr = NaN;
E.magerr = NaN;
E.magerr2 = NaN;
E.rms = NaN;
E.azgap = NaN;
E.usedPhases = NaN;
E.nmag = NaN;

E.id = "";
E.methodID = "";
E.earthModel = "";
E.evMode = "";
E.evStatus = "";
E.agencyID = "";
E.magtype = "";
E.magmethod = "";
E.evDescription = "";
E.evType = "";
E.authorID = "";

E.creationTime = NaN;
E.nPphases = NaN;
E.Pphases = NaN;
E.nSphases = NaN;
E.Sphases = NaN;
E.nMLv = NaN;
E.MLv = NaN;
E.nMjma = NaN;
E.Mjma = NaN;
E.nML = NaN;
E.ML = NaN ;
E.nMsBB = NaN;
E.MsBB = NaN;
E.nMwp = NaN;
E.Mwp = NaN;
E.nmB = NaN;
E.mB = NaN;
E.nmb = NaN;
E.mb = NaN;

%%
if N > 1
    E = repmat(E,N,1);
end
