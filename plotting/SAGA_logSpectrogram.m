function [h,aa,c2] = SAGA_logSpectrogram(S)
S = S(1);
if isnat(S.ref)
    disp('no data');
    h = [];
    return;
end

%%
S = differentiateWaveforms(S);
S = demeanWaveforms(S);
S = detrendWaveforms(S);
S = interpolateWaveforms(S);

tw = 0.02;
S = taperWaveforms(S,tw);
S = filterWaveforms(S,1/20);
S = intWaveforms(S);
S = detrendWaveforms(S);

d = S.d;
Fs = round(1./S.delta);
t = getTimeVec(S);

nMinutes = 10;
winlen = nMinutes*60*Fs;

overlapFactor = 1-(1/8); %0.75; 
noverlap = round(overlapFactor*winlen);

detrendFlag = true;
[dcut,~,endIndex] = cutWindows(d,winlen,noverlap,detrendFlag);
dcut = detrend(demean(detrend(dcut)));

tdown = t(endIndex);

npow2 = nextpow2(winlen);
nfft = 2^npow2;
maxNFFT = 2^14;
nfft = min([nfft maxNFFT]);

%%
[pxx,fxx] = pmtm(dcut.*parzenwin(winlen),{4,'trace'},nfft,Fs);
%[pxx,fxx] = pwelch(dcut,parzenwin(winlen),[],nfft,Fs);
%[pxx,fxx] = pwelch(dcut,[],[],nfft,Fs);
%pxx = bsxfun(@rdivide,pxx,max(abs(pxx)));
pxx = 10*log10(abs(pxx).^2);

%%
fxxI = fxx > 0;
fxx = fxx(fxxI);
pxx = pxx(fxxI,:);

%%
minTdown = min(tdown);
sss = seconds(tdown-minTdown);

%%
nSubs = 6;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
%%
%fig.Visible = 'off';
aa(1) = subplot(nSubs,1,1);
%t = seconds(t - min(t));
%d = d(winlen/2:end);

t = t(1:length(d));
t = datenum(t);

minT = min(t);
maxT = max(t);
plot(t,d,'k-','linewidth',2); zoom on;
c = colorbar;
c.Visible = 'off';
title(strcat(S(1).kstnm,',',S(1).kcmpnm));
set(aa(1), 'YLimSpec', 'Tight');
aa(1).XTickLabel = '';
xlim([minT maxT]);

%%
aa(2) = subplot(nSubs,1,(2:nSubs));
fxx2 = logspace(-2,log10(Fs/2),length(fxx)/4);
pxx = interp2(sss,fxx,pxx,sss,fxx2,'spline');

%%
fxxI = fxx2 > 1/10 & fxx2 <= Fs/2; %min([100 Fs/2]);
fxx2 = fxx2(fxxI);
pxx = pxx(fxxI,:);
disp(size(pxx));

%%
tdumb = datenum(seconds(sss)+minTdown);
%     disp(dn2dt(tdumb(1:4)));
%     disp(dn2dt(t(1:4)));

tdumb = tdumb - seconds((nMinutes*60*(1-overlapFactor))*2);
h = pcolor(datenum(tdumb),fxx2,pxx);
axis xy;
%ax = gca;
c2 = colorbar;
maxmax = max(max(pxx));
minmin = max([maxmax - 200, min(min(pxx))]);
caxis([minmin maxmax]);
h.EdgeColor = 'none';
aa(2).YScale = 'log';

%xlabel(['Hours since: ',datestr(minTdown)]);
%ylabel('Period [sec.]');
%xlabel(['seconds since: ',datestr(minTdown)]);
ylabel('frequency [hz.]');
%set(h, 'alphadata', pxx >= 20, 'facealpha', 'flat');
zoom on;
alpha 0.8;
colormap('jet');
xlabel('UTC'); %['seconds since: ',datestr(minTdown)]);
linkaxes(aa,'x');
aa(2).XTickLabel = datestr(aa(2).XTick);
aa(2).XTickLabelRotation = 20;
xlim([minT maxT]);
caxis([0 200]);
