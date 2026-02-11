clear; close all;

load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');

[stlat,stlon,stelev] = metaDataFromStationList(["SN11";"SN12"]);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(stlat,stlon,eqlat,eqlon,refEllipse)*1e-3;


%%
dI = t >= datetime(2018,04,01) & d_ < 10000 & d_ >= 100 & t < datetime(2019,04,01);
d2 = d_(dI);
S = extractWaveforms(t(dI)+seconds(d2/3.4),hours(2),"SN11","HHZ","9D","",true,true,'~/data/iguana/BROADBAND/');
T = extractWaveforms(t(dI),hours(2),"SN12","HHZ","9D","",true,true,'~/data/iguana/BROADBAND/');
d2 = d_(dI);
refs = pull(S,'ref');
rI = isnat(refs);
S(rI) = [];
d2(rI) = [];
refs(rI) = [];
T(rI) = [];
rI = isnat(pull(T,'ref'));
S(rI) = [];
d2(rI) = [];
T(rI) = [];
refs(rI) = [];