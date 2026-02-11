clear
close all
cd ~/research/now/reventador/
load CASC.mat
tCASC = trefs;
clear trefs tGood
CZ = ACZ; CN = ACN; CE = ACE;
clear A*
load BONI.mat
tBONI = trefs;
clear trefs tGood
refTime = min([min(tBONI) min(tCASC)]); [lia,locb] = ismembertol(seconds(tCASC-refTime),seconds(tBONI-refTime),2,'DataScale',1);
%
tCASC = tCASC(lia);
tBONI = tBONI(locb(lia));
CZ = CZ(lia); CN = CN(lia); CE = CE(lia);
BZ = BZ(locb(lia)); BN = BN(locb(lia)); BE = BE(locb(lia));
tw = 0.04; dOrig = cumsum(zpkFilter([taper(detrend(diff(double(pull(CZ)))),tw); taper(detrend(diff(double(pull(CN)))),tw); taper(detrend(diff(double(pull(CE)))),tw); taper(detrend(diff(double(pull(BZ)))),tw); taper(detrend(diff(double(pull(BN)))),tw); taper(detrend(diff(double(pull(BE)))),tw)],3/8,3/4,10,4,0));
close all; d = dOrig; figure(); semilogy(tCASC,0.5*peak2peak(d)','.'); zoom on; grid on;
figure(); semilogy(tCASC,0.5*peak2rms(d)','.'); zoom on; grid on; tmpStackAligned = plot_family(d(:,1:100),1:size(d(:,1:100),2),10,10);

d = [];
tw = 0.04;
lfc = 3/8;
hfc = 3/4;
npoles = 4;
zeroPhaseFlag = false;
newFs = 10;
d = [d; normalizeWaveforms(zpkFilter(taper(detrend(double(pull(CZ))),tw),lfc,hfc,newFs,npoles,zeroPhaseFlag))]; 
d = [d; normalizeWaveforms(zpkFilter(taper(detrend(double(pull(CN))),tw),lfc,hfc,newFs,npoles,zeroPhaseFlag))]; 
d = [d; normalizeWaveforms(zpkFilter(taper(detrend(double(pull(CE))),tw),lfc,hfc,newFs,npoles,zeroPhaseFlag))]; 
d = [d; normalizeWaveforms(zpkFilter(taper(detrend(double(pull(BZ))),tw),lfc,hfc,newFs,npoles,zeroPhaseFlag))]; 
d = [d; normalizeWaveforms(zpkFilter(taper(detrend(double(pull(BN))),tw),lfc,hfc,newFs,npoles,zeroPhaseFlag))]; 
d = [d; normalizeWaveforms(zpkFilter(taper(detrend(double(pull(BE))),tw),lfc,hfc,newFs,npoles,zeroPhaseFlag))]; 

d = normalizeWaveforms(d); 
tmpStackAligned = plot_family(d(:,1:100),1:size(d(:,1:100),2),10,10);

%%
clearvars -except d
[maxccp,plags,maxccn,nlags] = doccFreqCircShift(d,true);

clearvars -except max* *lags d
save('CASC_BONI_docc2');