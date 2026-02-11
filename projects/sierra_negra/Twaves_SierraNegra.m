clear; close all; clc;

load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');

kstnm = "SN11";
[stlat,stlon,stelev] = metaDataFromStationList(kstnm); %["SN11";"SN12"]);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(stlat(1),stlon(1),eqlat,eqlon,refEllipse)*1e-3;

dI = t >= datetime(2018,04,10) & d_ < 10000 & d_ >= 100 & t < datetime(2019,04,01);
d2 = d_(dI);
t2 = t(dI);

ll = sum(dI);
S = populateWaveforms(ll);
T = S;
U = S;
%for i = 1:ll
%     tStart = t2(i)+seconds(d2(i)/3.4);
%     tEnd = t2(i)+seconds(d2(i)/1.5) + minutes(10);
%     disp(tEnd-tStart)
 
    tStart = t2+seconds(d2/1.5)-seconds(5);
    tDur = minutes(10); %t2+seconds(d2/1.5) + minutes(10);
%     disp(tEnd-tStart)
    %%
%    tic;
%     S(i) = extractWaveforms(tStart,tEnd-tStart,"SN11","HHZ","9D","",true,true,'~/data/iguana/BROADBAND/');
%     toc;
%     T(i) = extractWaveforms(tStart,tEnd-tStart,"SN12","HHZ","9D","",true,true,'~/data/iguana/BROADBAND/');
%     toc;
%     U(i) = extractWaveforms(tStart,tEnd-tStart,"VCH1","HHZ","9D","",true,true,'~/data/iguana/BROADBAND/');
%    toc;
%end
% 

U = extractWaveforms(tStart,tDur,kstnm,["HHZ";"HHN";"HHE"],"9D","",true,true,'~/data/iguana/BROADBAND/');

% %%
% refs = pull(S,'ref');
% rI = isnat(refs);
% S(rI) = [];
% d2(rI) = [];
% refs(rI) = [];
% T(rI) = [];
% 
% %%
% rI = isnat(pull(T,'ref'));
% S(rI) = [];
% d2(rI) = [];
% T(rI) = [];
% refs(rI) = [];