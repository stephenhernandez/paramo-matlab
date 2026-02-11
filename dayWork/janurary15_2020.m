clear; close all; clc;
cd ~/research/now/pichincha_nxcorr/
load('ggp_2008_2020_PLD');

%%
pI = cc_all_vs_stack_2 >= 0.4;
tabs = tabs(pI);
p2p = p2p(pI);
p2rms = p2rms(pI);
indiv_events = indiv_events(:,pI);
maxAmpRMS = maxAmpRMS(pI);
iet = minutes(diff(tabs));

%%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tabs,1:length(tabs),'.');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
a(1) = subplot(121);
semilogy(a(1),tabs,p2p,'.');
a(2) = subplot(122);
plot(a(2),tabs,cc_all_vs_stack_2(pI),'.');
zoom on;

%%
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(121);
% semilogx(p2p(1:end-1),iet,'.');
% subplot(122);
% semilogx(iet,p2p(2:end),'.');
% zoom on;

%%
N = 301;
zeroPhaseFlag = true;
[p2p_,sI] = sort(p2p(1:end-1));

figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(121);
semilogy(log10(p2p_),iet(sI),'.');
hold on;
semilogy(ax(1),log10(medfiltSH(p2p_,N,zeroPhaseFlag)),medfiltSH(iet(sI),N,zeroPhaseFlag),'k-','linewidth',3);
title('time-predictable model');

[iet,sI] = sort(iet);
ax(2) = subplot(122);
maxAmpRMS_ = p2p(2:end);
semilogx(iet,log10(maxAmpRMS_(sI)),'.');
hold on;
semilogx(ax(2),medfiltSH(iet,N,zeroPhaseFlag),log10(medfiltSH(maxAmpRMS_(sI),N,zeroPhaseFlag)),'k-','linewidth',3);
zoom on;
title('slip-predictable model');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tabs,(1:length(tabs))'/length(tabs),'.');
hold on;
plot(tabs,cumsum(log10(p2p))/sum(log10(p2p)),'.'); zoom on;

%%
tI = tabs >= datetime(2014,05,03) & tabs <= datetime(2014,05,06);
t_ = tabs(tI);
S = extractWaveforms(t_-seconds(60),seconds(120),"PINO","SHZ");
d = double(pull(S));
tmp_stack_ = plot_family(d,1:size(d,2),40,100);