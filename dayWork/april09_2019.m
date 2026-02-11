%april 09, 2019
%clear;
clear;
cd ~/research/now/sn_eruption/
load VCH1_templates_from_James_picks.mat
%ieF = pull(filterSacData(intSacData(detrendSacData(differentiateSacData(S))),3,12));
ieF = pull(filterSacData(detrendSacData(differentiateSacData(S)),3,12));
ieF = ieF(300:1499,:);
t = pull(S,'ref');

% ieF = indiv_events;
% t = tabs;

[maxccp,plags,maxccn,nlags] = doccFreqDom(ieF,true);
%cd ~/research/now/sn_eruption/
%load vch1_igepn_event_cuts
close all;
%threshes = [0.8 0.8 0.8 0.75 0.75 0.75 0.7 0.7 0.7]; 
threshes = [0.9 0.8 0.7 0.9 0.8 0.7 0.9 0.8 0.7]; %95 0.95 0.95 0.90 0.90 0.90 0.85 0.85 0.85 0.80 0.80 0.80 ...
%    0.75 0.75 0.75 0.70 0.70 0.70]; %0.65 0.65 0.65 0.65 0.60 0.60];
method = 'weighted';

%%
i = 1;
ieF = constrictFamilies(maxccp,ieF,threshes(i),method);
[maxccp,plags,maxccn,nlags] = doccFreqDom(ieF,true);
for i = 2:length(threshes)
    ieF = constrictFamilies(maxccp,ieF,threshes(i),method);
    [maxccp,plags,maxccn,nlags] = doccFreqDom(ieF,true);
end

%%
maxccp = squareform(maxccp);
figure('units','normalized','outerposition',[0 0 1 1]);
h = imagesc(maxccp);
title('Similarity Matrix');
ylabel('Event Number');
xlabel('Event Number');
set(h, 'alphadata', maxccp >= -1);
axis square;
colorbar;
caxis([0 1]);
zoom on;

%%
figure(); plot(squareform(maxccp),'o'); zoom on;

%%
I1 = 237; I2 = 238; d1 = ieF(:,I1); d2 = ieF(:,I2);
[cc,lags] = xcorr(d1,d2,'coeff');
figure(); plot(lags/100,cc); zoom on;
figure(); plot(apply_vdcc(normalize_traces([d1 -d2]))); zoom on;

%%
function [new,singletons,family,l_uniq_indices] = constrictFamilies(maxccp,ieF,thresh,method)
[family,l_uniq_indices] = generateFamilies(maxccp,thresh,method);
new = ieF;
for ii = 1:sum(l_uniq_indices>1)
    new_trace = normalize_traces(sum(apply_vdcc(normalize_traces(ieF(:,family{ii}))),2));
    new(:,ii) = new_trace;
end

%%
singletons = cat(1,family{ii+1:end});
new(:,((ii+1):length(family))) = normalize_traces(ieF(:,singletons));
new = new(:,1:length(family));
end