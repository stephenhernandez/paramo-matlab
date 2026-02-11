clear; close all; clc; 

tic;
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = loadRepeaterCatalog("CHL1");
toc;

badI = ismember(templateNumber,[34; 94; 98; 181; 191; 71; 127; 217; 134; 125; 102; 182; 65; 77; 166; 177; 187; 173; 18; 140; 127]);

tMain = tMain(~badI);
ccMain = ccMain(~badI);
ampMain = ampMain(~badI);
dMag = dMag(~badI);
magMain = magMain(~badI);
templateNumber = templateNumber(~badI);
madMain = madMain(~badI);
nUsedMain = nUsedMain(~badI);

[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    filterCatalog(tMain,ccMain,10,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain);

ccI = (ccMain >= 0.4 | madMain >= 12) & nUsedMain >= 3 & ampMain >= 5; %4e2 & tMain >= datetime(2023,03,08);% & tMain <= datetime(2023,03,11);

%%
ampThresh = 3e2;
ccI2 = find(ccI & ampMain >= ampThresh & tMain >= datetime(2023,02,16));
%ccI2 = find(ccI & ampMain >= ampThresh & tMain >= datetime(2023,03,01) & tMain <= datetime(2023,06,01));
newFs = 100;
CHL1 = extractWaveforms(tMain(ccI2)-seconds(0),seconds(10),"CHL1","HHZ","EC","",true,true,...
    1,true,[1/20,-inf,false,false,false,newFs]);
d = pull(CHL1);
d(~isfinite(d)) = 0;
badI = logical(sum(~isfinite(d))') | rms(d)' == 0;
d(:,badI) = [];
ccI2(badI) = [];
t2 = tMain(ccI2);
a2 = ampMain(ccI2);
tNumber = templateNumber(ccI2);

tic;
fprintf('computing spectra for %d events, may take a while...\n',size(d,2));
%[pxxN,fxxN] = pmtm(normalizeWaveforms(detrend((d))),2,2^12,newFs);
[pxxN,fxxN] = pwelch(normalizeWaveforms(detrend(d)),2^7,[],2^12,newFs);
fI1 = fxxN >= 0.5 & fxxN < 2;
fI2 = fxxN >= 1 & fxxN < 4;
fI3 = fxxN >= 4 & fxxN <= 16;
lb = sum(pxxN(fI1,:))';
mb = sum(pxxN(fI2,:))';
hb = sum(pxxN(fI3,:))';
fIndex1 = lb./hb;
fIndex2 = mb./hb;
toc;

%%
close all;
figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(211); plot(tMain(ccI),t2r(tMain(ccI),hours(1)),'.'); zoom on; grid on;

%
figure(1); yyaxis right; plot(tMain(ccI),t2r(tMain(ccI),hours(24),[],false),'.'); zoom on; grid on; 
ylabel("Daily Rate"); legend("hourly","daily");
yyaxis left;
ax(1).YLabel.String = "Hourly Rate";
ax(2) = subplot(212); 
semilogy(tMain(ccI),ampMain(ccI),'.'); zoom on; grid on; hold on; 
semilogy(tMain(ccI),medfiltSH(ampMain(ccI),101),'.'); 
ax(2).YLabel.String = "Amplitude [counts]";
legend("","101-pt. moving median");
linkaxes(ax,'x');
figure('units','normalized','outerposition',[0 0 0.7 1]); 
semilogy(t2,fIndex1,'.'); zoom on; grid on; hold on; 
plot([min(t2) max(t2)],[0.5 0.5],'linewidth',4);
ylabel("frequency index"); 
title("$A_{0.5 \le f < 2 Hz.} / A_{4 \le f < 16 Hz.}$"); axis tight;

lpI = fIndex1 >= 0.3;
nlpI = find(~lpI);
lpI = find(lpI);

figure('units','normalized','outerposition',[0 0 0.7 1]); 
plot(t2(lpI),1:length(lpI),'.'); zoom on; grid on;
ylabel("cumulative number"); axis tight;
title("high ($\ge$0.5) FI events");

figure('units','normalized','outerposition',[0 0 0.7 1]); 
ll1 = loglog(fxxN,pxxN(:,nlpI(randsample(length(nlpI),500))),'linewidth',0.3,'color',[0.5 0.5 0.5]); 
zoom on; grid on; hold on; ll2 = loglog(fxxN,pxxN(:,lpI(randsample(length(lpI),500))),'linewidth',0.3,'color','k');
legend([ll1(1) ll2(1)],"$FI < 0.5$","$FI \ge 0.5$"); axis tight;
xlabel("frequency [hz.]");
ax = gca;
ax.YTickLabel = "";
xlim([0.5 20]); ylabel("normalized amplitude");
title("black: possible LP; gray: likely not LP");

clear ax; 
figure('units','normalized','outerposition',[0 0 0.8 1]); 
ax(1) = subplot(211); plot(t2,t2r(t2,hours(1),[],false),'.'); zoom on; grid on;
yyaxis right; 
plot(t2,t2r(t2,hours(24),[],false),'.'); zoom on; grid on; ylabel("Daily Rate");
legend("hourly","daily")
yyaxis left;
ax(1).YLabel.String = "Hourly Rate";

ax(2) = subplot(212); 
semilogy(t2,a2,'.'); zoom on; grid on; hold on; 
semilogy(t2,medfiltSH(a2,101),'.'); 
legend("","101-pt. moving median");
ax(2).YLabel.String = "Amplitude [counts]";
linkaxes(ax,'x'); axis tight;


vlpI = fIndex1 >= 10;
figure('units','normalized','outerposition',[0 0 0.7 1]); 
semilogy(t2(lpI),a2(lpI),'.'); zoom on; grid on; hold on; 
semilogy(t2(vlpI),a2(vlpI),'o','linewidth',2);
ylabel("amplitud [cuentas]"); 
legend("$0.5\le FI <10$","$FI \ge 10$"); 
title("Blue = Possible LP; Red = Possible VLP"); axis tight;

cd ~/research/now/chiles/
load("CVCCN_CHL1_March2023Templates_N233");
Torig = T;
Torig = Torig';

tLP = t2(lpI);
figure('units','normalized','outerposition',[0 0 1 1]); 
histogram(tLP,((dateshift(min(tLP),'start','day')):(dateshift(max(tLP),'start','day')))'+1); 
zoom on; grid on; title("Daily Number Possible LPs");

fI4 = fxxN >= 0.5 & fxxN < 20;
pxx4 = pxxN(fI4,:);
[~,maxI] = max(pxx4);
fxx4 = fxxN(fI4);
newLPI = (fxx4(maxI') > 0.513 & fxx4(maxI') < 2);

fxx4 = fxxN(fI4);
ydum = interp1(datenum(t2),t2r(t2,days(1)),datenum(tLP));

figure('units','normalized','outerposition',[0 0 1 1]); 
ax(1) = subplot(211); 
plot(tLP,t2r(tLP,hours(24)),'.'); 
zoom on; grid on; title("Daily Rate of Possible LPs");
cbar1 = colorbar;
cbar1.Visible = 'off';

ax(2) = subplot(212);
SS = scatter(ax(2),tLP,t2r(tLP,days(1))./ydum,5*exp(log10(ydum)),ydum,'o','filled'); zoom on; grid on; 
colorbar; LPax = gca; set(LPax,'ColorScale','log'); colormap turbo;
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5;
title(sprintf("Fraction of LPs Over All CVCCN Seismicity (Amp. $\\ge$ %d)",ampThresh)); 
linkaxes(ax,'x');

figure('units','normalized','outerposition',[0 0 1 1]);
SS = scatter(t2,fxx4(maxI'),0.01+exp(log10(rms(d)'))*5,fIndex1,'o','filled');
zoom on; grid on; ax = gca; ax.YScale = 'log';
ylim([0.5 20]); c = colorbar; SS.MarkerFaceAlpha = 0.5; clim([3e-3 3e1]);
c.Label.String = "Frequency Index"; c.Label.Interpreter = 'latex';
ylabel("Peak Frequency [hz.]"); set(ax,"ColorScale",'log');
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.2; colormap turbo; 

cd ~/research/now/chiles/
for i = 1:9
    fName = sprintf("Figure_%02d",i); 
    figure(i); print('-djpeg',fName); pause(1);
end
