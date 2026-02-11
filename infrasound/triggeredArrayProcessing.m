clear; close all;
lfc = 0.25;
hfc = 1;
S = loadWaveforms(datetime(2022,02,08),1,"SAG1","BDF","EC",["01";"02";"03";"04";"05"]);
Sf = resampleWaveforms(detrendWaveforms(intWaveforms(detrendWaveforms(filterWaveforms(...
    taperWaveforms(detrendWaveforms(differentiateWaveforms(S)),0.0004),lfc,hfc)))),50);
Scut = cutWaveforms(Sf,dateshift(Sf(1).ref,'start','day')+hours(00),0,hours(24));
% % close all; [tMain,ampMain,velMain,bazMain,medCCsMain,shiftsMain,waveformMain] = ...
% %     infrasoundArrayProcessing(Scut,-inf,-inf);

%%
close all;
tic;
secDur = 10;
maxRange = 2;
sta = 5;
lta = 5;
mph = 1.4;
hFlag = false;
plotFlag = false;
envFiltFlag = false;
hfc = -inf;
verboseFlag = true;

Scut = syncWaveforms(Scut);
locsMain = [];
for i = 1:length(Scut)
    locs_ = stalta(Scut(i),sta,lta,mph,hFlag,plotFlag,envFiltFlag,hfc,verboseFlag);
    locsMain = [locsMain; locs_];
end

%
lS = length(Scut);
%plotWaveforms(Scut);
t = getTimeVec(Scut);
tTrial = sort(t(locsMain));
% figure();
% plot(tTrial,(0:length(tTrial)-1)','.'); zoom on;
% 
rate = t2r(tTrial(:,1),seconds(maxRange));
% figure(); plot(tTrial(:,1),rate,'.'); zoom on; grid on;

%
fI = rate >= 3;
tPot = tTrial(fI);

delI = false(size(tPot));
diff_tPot = seconds(diff(tPot));
fI2 = find(diff_tPot < maxRange);

delI(fI2+1) = true;
tPot(delI) = [];

% figure();
% plot(tPot,(0:size(tPot,1)-1)','.'); zoom on; grid on;
% 
% t2 = tPot;
% [~,ax] = plotWaveforms(Scut);
% ax = ax(end);
% hold(ax,'on');
% ylim_ = ax.YLim;
% for i = 1:length(t2)
%     t_ = t2(i) - seconds(2);
%     plot(ax,[t_ t_],ylim_,'k','linewidth',2);
%     aa = area(ax,[t_; t_+seconds(secDur)],max(ylim_)*ones(2,1),min(ylim_));
%     aa.FaceAlpha = 0.5;
% end
toc;
