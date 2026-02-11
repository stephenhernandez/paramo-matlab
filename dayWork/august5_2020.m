clear; close all; clc;

%%load reventador excel sheet
%[tsheet,type,Tsp,amp,mag,codaDuration,period,refStation] = loadVolcanicEvents('reve',6);
cd ~/research/now/reventador/
load reventadorExcelData.mat

%% get filter
expI = strcmp(type,"EXP") & (strcmp(refStation,"REVS") | strcmp(refStation,"REVN")) & amp >= 2e3;

%% filter data
texp = tsheet(expI);
ampexp = amp(expI);

%% get another filter and filter data again
iet = seconds(diff(texp));
diffampexp = diff(ampexp);
amp1 = ampexp(1:end-1);
amp2 = ampexp(2:end);

dI = iet > 10 & iet <= 1e5;
iet = iet(dI);
diffampexp = diffampexp(dI);
amp1 = amp1(dI);
amp2 = amp2(dI);

%% generate plots
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(texp,ampexp,'.'); zoom on;

figure('units','normalized','outerposition',[0 0 1 1]); %slip-predictable?
loglog(iet,amp2,'.'); zoom on;
xlabel('interevent time');
ylabel('amplitude of second event');
title('test slip-predictable');

figure('units','normalized','outerposition',[0 0 1 1]); %time-predictable
loglog(amp1,iet,'.'); zoom on;
ylabel('interevent time');
xlabel('amplitude of first event');
title('test time-predictable');

% slip predictable?
b = robustfit(log10(iet),log10(amp2));
figure('units','normalized','outerposition',[0 0 1 1]);
plot(log10(iet),log10(amp2),'.'); zoom on; hold on; 
ysynth = b(1) + b(2)*log10(iet); 
plot(log10(iet),ysynth,'.');
xlabel('log10(interevent time)');
ylabel('log10(amplitude of first event)');
title({'fitting slip-predictability';['slope:',num2str(b(2))]});

% time predictable?
b = robustfit(log10(amp1),log10(iet));
figure('units','normalized','outerposition',[0 0 1 1]);
plot(log10(amp1),log10(iet),'.'); zoom on; hold on; 
ysynth = b(1) + b(2)*log10(amp1);
plot(log10(amp1),ysynth,'.');
ylabel('log10(interevent time)');
xlabel('log10(amplitude of first event)');
title({'fitting time-predictability';['slope:',num2str(b(2))]});

%%
% figure('units','normalized','outerposition',[0 0 1 1]);
% loglog(iet(dI),diffampexp(dI),'.'); zoom on;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% loglog(iet(dI),diffampexp(dI)./amp1(dI),'.'); zoom on;
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% loglog(iet(dI),diffampexp(dI)./amp2(dI),'.'); zoom on;