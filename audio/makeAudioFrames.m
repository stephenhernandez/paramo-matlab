function [S,h,hp,hp2,ax] = makeAudioFrames(varargin)

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Aug 07, 2019

%%
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2015,08,13),...
    60,...
    120,...
    0.25,...
    0.75,...
    12,...
    "BREF",...
    "BHZ",...
    "EC",...
    "",...
    '~/rawdata/'};

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[tStart,dur,nFrames,clipLevel,lfc,hfc,stnm,chan,net,locID,rawDataDir] = deal(optsToUse{:});

%%
tEnd = tStart + seconds(dur)*nFrames;
tvec = tStart:seconds(dur):tEnd;
tvec = tvec(1:nFrames)';

%%
S = extractWaveforms(tStart,tEnd,stnm,chan,net,locID,false,false,rawDataDir);

%%
cornersfin = isfinite([lfc hfc]);
if any(cornersfin)
    S = filterWaveforms(S,lfc,hfc);
end

%%
S = normalizeWaveforms(S,true,true);
S = clipWaveforms(S,clipLevel);
S = normalizeWaveforms(S,true,true);

%%
blur = cutWaveforms(S,repmat(tStart,nFrames,1),zeros(nFrames,1),tvec-tStart);
bold = cutWaveforms(S,repmat(tEnd,nFrames,1),tvec - tEnd,zeros(nFrames,1));

%%
PWD = pwd;
cd ~/animations/
dirName = strcat(stnm,".",locID,".",net,".",chan,".",datestr(tStart,30));
if ~exist(dirName,'dir')
    mkdir(dirName);
end
cd(dirName);

%%
close all;
[ax,hp] = plotWaveforms(S);
hp.LineWidth = 1;
frameName = strcat('frame_',num2str(1));
ylim([-1 1]);
title(datestr(S.ref));
print(frameName,'-djpeg','-r100');
h = gcf;
h.Visible = 'off';
set(h,'defaultLegendAutoUpdate','off');
delete(hp);
ax.ColorOrderIndex = 1;

%%
for i = 2:nFrames
    disp(i);
    frameName = strcat('frame_',num2str(i));
    [~,hp] = plotWaveforms(ax,blur(i),[],[],'k-');
    hp.LineWidth = 0.1;
    [~,hp2] = plotWaveforms(ax,bold(i));
    hp2.LineWidth = 1;
    ylim([-1 1]);
    title(datestr(bold(i).ref));
    print(frameName,'-djpeg','-r100');
    delete(hp);
    delete(hp2);
    ax.ColorOrderIndex = 1;
end


%%
[~,hp] = plotWaveforms(ax,S,[],[],'k-');
hp.LineWidth = 0.1;
ylim([-1 1]);
title(datestr(S.ref+S.e));
frameName = strcat('frame_',num2str(nFrames+1));
print(frameName,'-djpeg','-r100');

%%
cd(PWD);
