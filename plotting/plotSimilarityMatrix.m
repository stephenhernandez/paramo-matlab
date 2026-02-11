function plotSimilarityMatrix(indiv_events,thresh,Fs,lfc,hfc,zeroPhaseFlag)
if nargin < 2
    thresh = 0.7;
end
if nargin < 3
    Fs = 100;
end
if nargin < 4
    lfc = -inf;
end
if nargin < 5
    hfc = -inf;
end
if nargin < 6
    zeroPhaseFlag = false;
end

%%
npoles = 4;
if any(isfinite([lfc hfc]))
    disp('filtering data');
    indiv_events = zpkFilter(indiv_events,lfc,hfc,Fs,npoles,zeroPhaseFlag);
else
    disp('no filtering requested');
end

%%
indiv_events = normalizeWaveforms(indiv_events);
[maxccp,~,maxccn] = doccFreqCircShift(indiv_events,true);
mI = maxccn > maxccp;
maxccp(mI) = maxccn(mI);
maxccp = squareform(maxccp);
mI = maxccp == 0;
maxccp(mI) = NaN;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
h = imagesc(maxccp);
title('Similarity Matrix');
ylabel('Event Number');
xlabel('Event Number');
set(h, 'alphadata', maxccp >= thresh);
axis square;
colorbar;
clim([thresh 1]);
zoom on;