%april15_2019
clear; close all; clc;
cd ~/research/now/reve_silvi/
load infrasound_waveform_similarity.mat
t = tabs;

% cd ~/research/now/sn_eruption/
% load VCH1_templates_from_James_picks.mat
% indiv_events = pull(filterSacData(detrendSacData(differentiateSacData(S)),3,12));
% indiv_events = indiv_events(300:1499,:);
% t = pull(S,'ref');
% [maxccp,plags,maxccn,nlags] = doccFreqDom(indiv_events,true);

%%
threshes = 0.9; %flipud((0.7:0.01:0.99)'); %[0.95 0.9 0.85 0.8 0.75 0.7 0.7 0.7 0.6 0.6 0.5 0.5 0.5]; % 0.6 0.5 0.5];
%method = ["single"]; %, "average", "weighted"];
method = ["average"]; %,"average","single"];

for i = 1:length(method)
    method_ = method(i);
    for j = 1:length(threshes)
        disp([method_,', ',num2str(threshes(j))]);
        [indiv_events,ii,mccpCopy] = constrictFamilies(maxccp,indiv_events,threshes(j),method_);
        if ii
            maxccp_ = doccFreqDom(indiv_events,true,ii);
            maxccp = [maxccp_; squareform(mccpCopy)'];
            clear maxccp_ mccpCopy
        end
    end
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
I1 = 56;
I2 = 62;
d1 = indiv_events(:,I1); d2 = indiv_events(:,I2);
[cc,lags] = xcorr(d1,d2,'coeff');
figure(); plot(lags/100,cc); zoom on;
figure(); plot(apply_vdcc(normalize_traces([d1 d2]))); zoom on;

%%
function [indiv_events,numConstrictedFams,mccpCopy] = constrictFamilies(maxccp,indiv_events,thresh,method)
[family,l_uniq_indices] = generateFamilies(maxccp,thresh,method);

new = indiv_events;
numConstrictedFams = sum(l_uniq_indices>1);
if numConstrictedFams
    disp(' ');
    disp(' ');
    disp(' ');
    disp(' ');
    disp(['Number of constricted families: ',num2str(numConstrictedFams)]);
    disp(' ');
    disp(' ');
    disp(' ');
    disp(' ');
    for ii = 1:sum(l_uniq_indices>1)
        new_trace = normalize_traces(sum(apply_vdcc(normalize_traces(indiv_events(:,family{ii}))),2));
        new(:,ii) = new_trace;
    end
    
    %%
    singletons = cat(1,family{ii+1:end});
    new = new(:,1:ii);
    mccpCopy = maxccp;
    mccpCopy = squareform(mccpCopy);
    mccpCopy = mccpCopy(:,singletons);
    mccpCopy = mccpCopy';
    mccpCopy = mccpCopy(:,singletons);
    singletons = indiv_events(:,singletons);
    indiv_events = [new singletons];
else
    mccpCopy = [];
    return;
end
end
