function S = populateWaveforms(sizeS)
if nargin < 1; sizeS = 1; end
S.ref = NaT;
S.b = seconds(0);
S.e = seconds(0);
S.npts = NaN;
S.delta = NaN;
S.d = [];

%%
S.kstnm = "";
S.kcmpnm = "";
S.knetwk = "";
S.khole = "";
S.evid = "";

%%
S.evla = NaN;
S.evlo = NaN;
S.evdp = NaN;
S.eqmag = NaN;
S.stla = NaN;
S.stlo = NaN;
S.stel = NaN;

%%
S.depmin = NaN;
S.depmax = NaN;
S.depmen = NaN;
S.az = NaN;
S.baz = NaN;
S.gcarc = NaN;
S.dist = NaN;

%%
S.gapFlag = false;
S.gapInfo = [];

%%
S.cmpinc = NaN;
S.cmpaz = NaN;
%S.timeCorrection = NaN;

%%
S.user0 = NaN;
S.user1 = NaN;
S.user2 = NaN;
S.user3 = NaN;
S.user4 = NaN;
S.user5 = NaN;
S.user6 = NaN;
S.user7 = NaN;
S.user8 = NaN;
S.user9 = NaN;

%% experimental filter and response information
S.lfc = -inf;
S.hfc = -inf;
S.npoles = NaN;
S.zeroPhaseFlag = false;
S.constant = NaN;
S.gain = NaN;
S.sensitivity = NaN;
S.A0 = NaN;
S.zeroes = [];
S.poles = [];

%%
if isscalar(sizeS)
    S = repmat(S,sizeS,1); %make column vector
else
    S = repmat(S,sizeS);
end