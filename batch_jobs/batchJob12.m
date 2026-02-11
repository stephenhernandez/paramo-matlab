clear; close all; clc;

load ~/research/now/reventador/ReventadorSubspaceDetectorResults_v10.mat
[stla,stlo] = metaDataFromStationList(["CASC";"BONI";"ANTS";"ANTG"]);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(stla,stlo,-0.080850,-77.657995,refEllipse)*1e-3;
CASC = 1e9*median(z2p(:,1:3),2,"omitnan")/3.141950e+08;

bI = tabs >= datetime(2020,01,338);
BONI = 1e9*median(z2p(:,4:6),2,"omitnan");
BONI(bI) = BONI(bI)/5.037350e+08;
BONI(~bI) = BONI(~bI)/2.014940e+09;

bI = tabs >= datetime(2014,01,226);
ANTS = 1e9*median(z2p(:,7:9),2,"omitnan");
ANTS(bI) = ANTS(bI)/5.037350e+08;
ANTS(~bI) = ANTS(~bI)/3.141950e+08;

bI = tabs >= datetime(2017,01,105);
ANTG = 1e9*median(z2p(:,10:12),2,"omitnan");
ANTG(bI) = ANTG(bI)/3.141950e+08;
ANTG(~bI) = ANTG(~bI)/3.141950e+08;

amps = [CASC BONI ANTS ANTG];
fourI = sum(isfinite(amps),2) == 4;
fourI = fourI & CASC >= 100;
fourI = fourI & BONI >= 40;
fourI = fourI & ANTS >= 80;
fourI = fourI & ANTG >= 40;

tabs2 = tabs(fourI); 
z2p2 = CASC(fourI,1);  
[tabs3,z2p3] = filterCatalog(tabs2,z2p2(:,1),30);

%
load reventadorExcelData.mat;
refI = refStation == "REVS" & amp >= 2000 & tsheet >= datetime(2013,01,01);
tsheet2 = tsheet(refI);
amp2 = amp(refI);

[tsheet3,amp3] = filterCatalog(tsheet2,amp2,30);

[t1_commonI,t2_commonI,just_t1I,just_t2I] = synchronizeCatalog(tsheet3,tabs3,30,true,false);

figure(); 
semilogy(tsheet3(t1_commonI),amp3(t1_commonI),'.'); zoom on; grid on; hold on; 
semilogy(tabs3(t2_commonI),z2p3(t2_commonI),'.');

%%
sampleIndex = sort(randsample(length(t2_commonI),200)); 
tsample = tabs3(sampleIndex);

ampsCorrected = NaN(length(tsample),10);

lfc = 0.375; hfc = 0.75;
A_ = extractWaveforms(tsample,seconds(60),["CASC"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,true,true,false]); 
ampsCorrected(:,1) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["CASC"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,false,true,false]); 
ampsCorrected(:,2) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["BONI"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,true,true,false]); 
ampsCorrected(:,3) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["BONI"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,false,true,false]); 
ampsCorrected(:,4) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["ANTS"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,true,true,false]); 
ampsCorrected(:,5) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["ANTS"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,false,true,false]); 
ampsCorrected(:,6) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["ANTG"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,true,true,false]);
ampsCorrected(:,7) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["ANTG"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,false,true,false]);
ampsCorrected(:,8) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["REVS"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,true,true,false]); 
ampsCorrected(:,9) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

A_ = extractWaveforms(tsample,seconds(60),["REVS"],"HHZ","EC",[""],true,true,1,true,[lfc,hfc,false,true,false]);
ampsCorrected(:,10) = (pull(A_,'depmax') - pull(A_,'depmin'))*0.5;

figure(); 
semilogy(tsample,ampsCorrected,'.'); zoom on;

clearvars -except tsample ampsCorrected
save('~/research/now/reventador/CASC_BONI_ANTS_ANTG_REVS_ampsCorrected');
