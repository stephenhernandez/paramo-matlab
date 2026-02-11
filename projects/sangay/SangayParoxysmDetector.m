% SangayParoxysmDetector
clear; close all;

%[z2p,tcut] = rsam2magnitude(15,datetime(2020,01,01),datetime(2021,10,20)); % 15 minutes == 30 samples
[z2p,tcut] = rsam2magnitude(90,datetime(2021,01,01),datetime(2021,10,20)); % 30 minutes == 60 samples
z2p = z2p(1:length(tcut),:);
[Msangay,Err_sangay,mlv,wa,dmlv] = sangayMagnitudeCalculation(tcut,z2p,false);

close all; 
figure(); ax_(1) = subplot(211); 
plot(tcut,Msangay,'.'); zoom on; 
ax_(2) = subplot(212); 
semilogy(tcut,Err_sangay,'.'); 
linkaxes(ax_,'x');

Neff = 9 - sum(~isfinite(z2p(1:length(tcut),:)),2);
goodI = Err_sangay <= 0.11 & Neff >= 4;

% 
%%
% figure(); histogram(Msangay(goodI),200); zoom on; grid on;
% goodI = Err_sangay <= 0.11 & Neff >= 4 & Msangay >= 0.8 & Msangay <= 2;
% 
% tthresh = tcut(goodI);
% rateGood = t2r(tthresh,minutes(10));
% figure(); plot(tthresh,rateGood,'.'); zoom on; grid on;
% figure(); ax_(1) = subplot(211); plot(tcut(goodI),Msangay(goodI),'.'); zoom on; ax_(2) = subplot(212); semilogy(tcut(goodI),Err_sangay(goodI),'.'); linkaxes(ax_,'x');
% figure(); semilogy(tcut(goodI),z2p(goodI,:),'.'); zoom on; grid on;
% [Mcut,startIndex,endIndex,badFlag,nwindows] = cutWindows(Msangay(goodI),30,29,false);
% tcutgood = dn2dt(cutWindows(datenum(tcut(goodI)),30,29,false));
% tcutgood = seconds(tcutgood - tcutgood(1,:));
% bb = NaN(size(Mcut,2),1);
% for i = 1:size(tcutgood,2)
%     bb_ = robustfit(tcutgood(:,i),Mcut(:,i));
%     disp(i);
%     bb(i) = 1e6*bb_(2);
% end
% tgood = tcut(goodI);
% tEnd = tgood(endIndex);
% figure(); plot(tEnd(bb>= 50),bb(bb>50),'.'); zoom on; grid on;