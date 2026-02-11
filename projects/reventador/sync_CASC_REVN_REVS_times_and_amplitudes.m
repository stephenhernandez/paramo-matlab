clear; close all;

cd ~/research/now/reventador/
load LizGoodData_03DEC2021_2018Amps.mat

refTime = min([min(casc_t) min(revn_t) min(revs_t)]);
otherT = revn_t;
otherAmp = revn_amps;
otherkstnm = "REVN";

[lia,locb] = ismembertol(seconds(casc_t-refTime),seconds(otherT-refTime),2,'DataScale',1);

%close all;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(sort(casc_t(lia)),(1:length(casc_t(lia)))'/length(casc_t(lia)),'.');
zoom on; hold on;
plot(sort(otherT(locb(lia))),(1:length(otherT(locb(lia))))'/length(otherT(locb(lia))),'.');

[cascAmpSort,sI] = sort(1e9*casc_amps(lia)); 
otherAmpSort = 1e9*otherAmp(locb(lia)); 
otherAmpSort = otherAmpSort(sI);

%sI = revnAmpSort >= 1e4 & 
sI = cascAmpSort <= 1e4;
otherAmpSort = otherAmpSort(sI);
cascAmpSort = cascAmpSort(sI);

figure(); 
loglog(cascAmpSort,otherAmpSort,'.'); zoom on; grid on;

b_ = robustfit(log10(cascAmpSort),log10(otherAmpSort));
b_ = flipud(b_);
yq = polyval(b_,log10(cascAmpSort));

figure(4); 
hold on; 
ll = loglog(cascAmpSort,10.^yq,'linewidth',2); ll.Color(4) = 0.6; title(['slope = ',num2str(b_(1))]);

xlabel('CASC Amplitude [nm]');
ylabel([otherkstnm,' Amplitude [nm]']);
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';