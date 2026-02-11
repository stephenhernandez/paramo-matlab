clear; close all; clc;
tStart = datetime(2020,06,01);
tEnd = datetime(2020,06,30);
S = nanGapWaveforms(rmsGather(tStart,tEnd,60,2,4,"SAGA","HHZ",true,false,"EC",""),0); % true-false is mean of abs value
W4 = S;
S = nanGapWaveforms(rmsGather(tStart,tEnd,60,1/4,1,"SAGA","HHZ",true,false,"EC",""),0); % true-false is mean of abs value
S(2,1) = W4;
S = syncWaveforms(S);

%%
N = 30; 

%%
S2 = nanGapWaveforms(convWaveforms(nanGapWaveforms(medfiltWaveforms(S,1,true),0),N));
d1 = S2(1).d; %low band
d2 = S2(2).d; %high band
t = getTimeVec(S2);
inCorr = nanmedian(d2./d1);

close all; 
figure(); scatter(t,d2./d1./inCorr,[],log10(d1),'filled'); 
Cbar = colorbar; zoom on; grid on;
ax = gca; ax.YScale = 'log';
caxis([2 4]);

[~,ax] = plotWaveforms(flipud(S2)); for i = 1:length(ax); ax(i).YScale = 'log'; grid(ax(i),'on'); end; linkaxes(ax,'xy')
figure(); 
histogram(d2./d1./inCorr,logspace(floor(log10(min(d2./d1./inCorr))),ceil(log10(max(d2./d1./inCorr))),1001)); 
zoom on; grid on; ax = gca; ax.XScale = 'log';

figure(); scatter(t,d1,[],1./(d2./d1./inCorr),'filled'); 
Cbar = colorbar; zoom on; grid on;
axC = gca; axC.YScale = 'log';
caxis([0.5 2]);

figure(); scatter(t,d2,[],1./(d2./d1./inCorr),'filled'); 
Cbar = colorbar; zoom on; grid on;
axC(2,1) = gca; axC(2,1).YScale = 'log';
caxis([0.5 2]); linkaxes(axC,'x');