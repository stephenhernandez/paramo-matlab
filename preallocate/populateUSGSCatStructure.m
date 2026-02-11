function U = populateUSGSCatStructure(N)
if nargin < 1; N = 1; end
U.mag = NaN;
U.lat = NaN;
U.lon = NaN;
U.depth = NaN;
U.place = '';
U.t = NaT;
U.updated = NaT;
U.tz = NaN;
U.url = '';
U.detail = '';
U.felt = NaN;
U.cdi = NaN;
U.mmi = NaN;
U.alert = '';
U.status = '';
U.tsunami = NaN;
U.sig = NaN;
U.net = '';
U.code = '';
U.ids = '';
U.sources = '';
U.types = '';
U.nst = NaN;
U.dmin = NaN;
U.rms = NaN;
U.gap = NaN;
U.magType = '';
U.type = '';
U.title = '';

%%
if N > 1
    U = repmat(U,N,1);
end