clear; close all; clc;
cd ~/research/now/pichincha/pichincha_nxcorr/
load ggp_swarm_analysis_2008_2020.mat

%%
lfc = 2.5;
hfc = 10;
indiv_events = cumsum(zpkFilter(taper(diff(indiv_events),0.04),lfc,hfc,100));
% indiv_events = zpkFilter(taper(diff(indiv_events),0.04),lfc,hfc,100);
figure(); semilogy(tabs,rms(indiv_events)','o'); zoom on;

%%
stack = normalizeWaveforms(nanmean(normalizeWaveforms(indiv_events),2));
norms = sqrt(size(indiv_events,1))*rms(indiv_events)';
figure(); plot(tabs,norms,'o');
estimates = repmat(norms',size(indiv_events,1),1).*repmat(stack,1,size(indiv_events,2));
noise = indiv_events - estimates;

cc_all_vs_stack_2 = sum(normalizeWaveforms(repmat(stack,1,size(indiv_events,2))).*normalizeWaveforms(indiv_events))';
figure(); plot(cc_all_vs_stack_2,'o'); zoom on;

%%
figure()
plot(tabs,cc_all_vs_stack,'o'); zoom on;

pI = cc_all_vs_stack < 0.3 & cc_all_vs_stack > 0.2;
tmp_stack_ = plot_family(indiv_events(:,pI),1:sum(pI),10,100);
% pI = cc_all_vs_stack >= 0.9;
% tmp_stack_ = plot_family(shifted_data(:,pI),1:sum(pI),10,100);

pI = cc_all_vs_stack > -1;
tabsOrig = tabs;
maxAmpRMSOrig = maxAmpRMS;

tabs = tabs(pI);
maxAmpRMS = maxAmpRMS(pI);

%%
dateshift(min(tabs),'start','hour');
iet = seconds(diff(tabs))/60;
figure();
a(1) = subplot(131);
semilogy(tabs(2:end),(iet),'.'); zoom on;

a(2) = subplot(132);
plot(tabs,1:length(tabs),'.'); zoom on;

a(3) = subplot(133);
semilogy(tabs,maxAmpRMS,'.'); zoom on;
linkaxes(a,'x');

%%
close all;
figure();
subplot(121);
semilogx(iet,maxAmpRMS(1:end-1),'.');
subplot(122);
semilogx(iet,maxAmpRMS(2:end),'.');
zoom on;

figure();
N = 251;
zeroPhaseFlag = true;
[iet,sI] = sort(iet);

maxAmpRMS_ = maxAmpRMS(1:end-1);
ax(1) = subplot(121);
loglog(iet,maxAmpRMS_(sI),'.');
hold on;
loglog(ax(1),medfiltSH(iet,N,zeroPhaseFlag),medfiltSH(maxAmpRMS_(sI),N,zeroPhaseFlag),'k-','linewidth',3);
ax(2) = subplot(122);
maxAmpRMS_ = maxAmpRMS(2:end);
loglog(iet,maxAmpRMS_(sI),'.');
hold on;
loglog(ax(2),medfiltSH(iet,N,zeroPhaseFlag),medfiltSH(maxAmpRMS_(sI),N,zeroPhaseFlag),'k-','linewidth',3);
zoom on;
linkaxes(ax,'xy');

%%
dayFraction = seconds(tabs - dateshift(tabs,'start','day'))/86400;
figure();
subplot(121);
histogram(dayFraction,100); zoom on;
subplot(122);
plot(sort(dayFraction),(1:length(dayFraction))'/length(dayFraction),'.'); zoom on;

%%
pseudoMag = log10(maxAmpRMSOrig);
Mw = (2/3)*pseudoMag + 1.15;
m0 = mw2m0(Mw);

figure();
subplot(211);
pseudoMag = log10(maxAmpRMSOrig);
Mw = (2/3)*pseudoMag + 1.15;
m0 = mw2m0(Mw);
semilogy(tabsOrig,m0,'.'); zoom on;
hold on;
pI = cc_all_vs_stack < 0.4;
tabs = tabsOrig(pI);
maxAmpRMS = maxAmpRMSOrig(pI);
pseudoMag = log10(maxAmpRMS);
Mw = (2/3)*pseudoMag + 1.15;
m0 = mw2m0(Mw);
semilogy(tabs,m0,'.'); zoom on;

subplot(212);
pseudoMag = log10(maxAmpRMSOrig);
Mw = (2/3)*pseudoMag + 1.15;
m0 = mw2m0(Mw);
semilogy(tabsOrig,cumsum(m0),'.'); zoom on;
hold on;
pI = cc_all_vs_stack < 0.4;
tabs = tabsOrig(pI);
maxAmpRMS = maxAmpRMSOrig(pI);
pseudoMag = log10(maxAmpRMS);
Mw = (2/3)*pseudoMag + 1.15;
m0 = mw2m0(Mw);
semilogy(tabs,cumsum(m0),'.'); zoom on;
