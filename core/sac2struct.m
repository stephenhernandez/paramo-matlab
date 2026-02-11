function S = sac2struct(fname)
%
% sac2struct read sac file and populate into waveform struct
%
% S = sac2struct(fname)
% We assume the file fname contains the path to a single sac files to
% be read and converted to a MATLAB structure.
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
S = populateWaveforms();

% Read individual files one by one
c = 1;
[d_,h,k] = readsacfile(fname);

yyyy_ = h(29,c);
doy_ = h(30,c);
hour_ = h(31,c);
minute_ = h(32,c);
second_ = h(33,c) + h(34,c)/1000;
delta = h(1,c);

%%
ref = datetime(yyyy_,01,doy_)+hours(hour_)+minutes(minute_)+seconds(second_);
b = seconds(h(4,c));

%%
kstnm = k(1,c)'; %kstnm = kstnm{1};
khole = k(3,c)'; %khole = khole{1};
kcmpnm = k(2,c)'; %kcmpnm = kcmpnm{1};
knetwk = k(4,c)'; %knetwk = knetwk{1};

stla = h(17,c)';    stlo = h(18,c)';
evla = h(21,c)';    evlo = h(22,c)';    evdp = h(24,c)';
dist = h(25,c)';    az = h(26,c)';
baz = h(27,c)';     gcarc = h(28,c)';
cmpaz = h(15,c)';   cmpinc = h(16,c)';

%%
delta = delta(c);
if b == 0
    S.ref = ref;
else
    S.ref = ref + b;
end

%%
S.b = seconds(0);
S.kstnm = upper(kstnm);
S.khole = upper(khole);
S.kcmpnm = upper(kcmpnm);
S.knetwk = upper(knetwk);
S.stla = stla;
S.stlo = stlo;
S.evla = evla;
S.evlo = evlo;
S.evdp = evdp;
S.dist = dist;
S.az = az;
S.baz = baz;
S.gcarc = gcarc;
S.cmpaz = cmpaz;
S.cmpinc = cmpinc;
S.delta = delta;
S.d = d_;
S.npts = length(d_);
S.e = seconds((S.npts - 1)*S.delta);

%%
[minVals,maxVals,meanVals] = minmaxmean(d_);
S.depmin = minVals;
S.depmax = maxVals;
S.depmen = meanVals;
