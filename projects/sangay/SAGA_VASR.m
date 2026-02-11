% sangay VASR


% % % clear
% % % [tabs,z2p,NCC,Neff,p2rms,kurt] = filterUniqueEvents('~/research/now/sangay/sangaySubspaceDetectorSAGA_v2');
% % % uniqDays = unique(dateshift(tabs,'start','day'));
% % % 
% % % n = 1;
% % % z2p2 = NaN(length(tabs),1); t2 = NaT(length(tabs),1);
% % % for i = 1:length(uniqDays)
% % % day_ = uniqDays(i);
% % % disp(datestr(day_));
% % % tI = tabs >= day_ & tabs < day_+1;
% % % S = extractWaveforms(tabs(tI)+seconds(30),seconds(120),"SAGA","HHZ","EC","",true,false);
% % % refs = pull(S,'ref');
% % % rI = isnat(refs);
% % % refs(rI) = [];
% % % S(rI) = [];
% % % nref = length(refs);
% % % if nref
% % % d2 = cumsum(detrend(zpkFilter(taper(detrend(demean(diff(double(pull(resampleWaveforms(S,100)))))),0.04),0.5,4,100,6)));
% % % z2p_ = max(abs(d2),[],1)';
% % % t2(n:n+nref-1) = refs;
% % % z2p2(n:n+nref-1) = z2p_;
% % % n = n + nref;
% % % end
% % % end
% % % z2p2 = z2p2(1:n-1);
% % % t2 = t2(1:n-1);
% % % 
% % % n = 1;
% % % z2pBDF = NaN(length(tabs),1); tBDF = NaT(length(tabs),1);
% % % for i = 1:length(uniqDays)
% % % day_ = uniqDays(i);
% % % disp(datestr(day_));
% % % tI = tabs >= day_ & tabs < day_+1;
% % % S = extractWaveforms(tabs(tI)+seconds(30),seconds(120),"SAGA","BDF","EC","01",true,false);
% % % refs = pull(S,'ref');
% % % rI = isnat(refs);
% % % refs(rI) = [];
% % % S(rI) = [];
% % % nref = length(refs);
% % % if nref
% % % d2 = cumsum(detrend(zpkFilter(taper(detrend(demean(diff(double(pull(resampleWaveforms(S,100)))))),0.04),0.6,2.4,100,6)));
% % % z2p_ = max(abs(d2),[],1)';
% % % tBDF(n:n+nref-1) = refs;
% % % z2pBDF(n:n+nref-1) = z2p_;
% % % n = n + nref;
% % % end
% % % end
% % % z2pBDF = z2pBDF(1:n-1);
% % % tBDF = tBDF(1:n-1);
% % % figure(); semilogy(t2,z2p2,'.'); zoom on; grid on;
% % % clearvars -except t2 z2p2 tBDF z2pBDF
% % % save('SAGA_VASR');

clear; 

cd ~/research/now/sangay/

%%
load SAGA_VASR.mat
z2p2(t2 <= datetime(2018,11,11)) = z2p2(t2 <= datetime(2018,11,11))*4;
close all; zI = z2p2 >= 1e4; figure(); plot(t2(zI),1:sum(zI),'.'); zoom on;
t3 = t2(zI);
z2p3 = z2p2(zI);
tI = ismembertol(seconds(t3-min([min(t3) min(tBDF)])),seconds(tBDF-min([min(t3) min(tBDF)])),1,'DataScale',1);
sum(tI)
t4 = t3(tI);
z2p4 = z2p3(tI);
tI = ismembertol(seconds(tBDF-min([min(t4) min(tBDF)])),seconds(t4-min([min(t4) min(tBDF)])),1,'DataScale',1);
sum(tI)
z2pBDF2 = z2pBDF(tI);
VASR = z2pBDF2./z2p4;
figure(); semilogy(t4,VASR,'.'); zoom on; grid on;
