clear; close all; clc;

%
cd ~/research/now/sangay/
% load SAGA_WoodAnderson.mat
% load sangayWoodAndersonAmplitudes.mat
load sangayRegionalDisplacementsNM

%%
%kstnms = ["SAGA";"PUYO";"TAIS";"TAMH";"PORT";"PKYU";"BMAS";"BULB";"BPAT";"BRUN"];
kstnms = ["SAGA"; "BPAT"; "BMAS"; "BULB"; "BRUN"; "PUYO"; "TAMH"; "PORT"; "PKYU"; "TAIS"];
lK = length(kstnms);

%%
refEllipse = referenceEllipsoid('wgs84');
[stla,stlo,stel] = metaDataFromStationList(kstnms);
dist = distance(stla,stlo,-2.00535,-78.341294,refEllipse)*1e-3;
[dist,sI] = sort(dist);
kstnms = kstnms(sI);

%%
%allAmps = [z2pWA z2pWA_PUYO z2pWA_TAIS z2pWA_TAMH z2pWA_PORT z2pWA_PKYU...
%    z2pWA_BMAS z2pWA_BULB z2pWA_BPAT z2pWA_BRUN];

%%
experimentalFlag = true;
if experimentalFlag
    load('raw2wa_sangayRegionalSensors_2','b','sagab'); %version 2 little bit better?
    %load('raw2wa_sangayRegionalSensors','b','sagab'); %version 1
    allAmps = [z2pWA_SAGA_raw z2pWA_BPAT_raw z2pWA_BMAS_raw z2pWA_BULB_raw...
        z2pWA_BRUN_raw z2pWA_PUYO_raw z2pWA_TAMH_raw z2pWA_PORT_raw...
        z2pWA_PKYU_raw z2pWA_TAIS_raw];
    for i = 1:lK
        if i == 1
            raw = allAmps(:,i);
            rI = raw >= 0 & isfinite(raw) & refs <= datetime(2018,11,11);
            raw(rI) = 10.^(sagab(1,2) + sagab(2,2)*log10(raw(rI)*4));
            allAmps(rI,i) = raw(rI);
            % raw(~rI) = 10.^(sagab(1,2) + sagab(2,2)*log10(raw(~rI)));
            raw(~rI) = 10.^(sagab(1,1) + sagab(2,1)*log10(raw(~rI)));   % version 2, using nanometers
            allAmps(~rI,i) = raw(~rI);
        else
            raw = allAmps(:,i);
            rI = raw >= 0 & isfinite(raw);
            raw(rI) = 10.^(b(1,i-1) + b(2,i-1)*log10(raw(rI)));
            allAmps(:,i) = raw;
        end
    end
else
    allAmps = [z2pWA_SAGA z2pWA_BPAT z2pWA_BMAS z2pWA_BULB z2pWA_BRUN z2pWA_PUYO...
        z2pWA_TAMH z2pWA_PORT z2pWA_PKYU z2pWA_TAIS];
end

allAmps = allAmps(:,sI);
goodI = allAmps > 1e-3;
allAmps(~goodI) = NaN;

%%
mlvRichter = NaN(size(allAmps));

%%
figure();
for i = 1:lK
    a = allAmps(:,i);
    if strcmp("BULB",kstnms(i))
        bI = refs >= datetime(2016,01,32) & refs <= datetime(2017,01,91);
        a(bI) = a(bI)/4;
        allAmps(:,i) = a;
    end
    plot(refs,a,'.','DisplayName',kstnms(i)); zoom on; grid on; hold on;
    mlvRichter(:,i) = log10(a) + 1.11*log10(dist(i)) + 0.00189*dist(i) + 0.591;
end
legend('show','location','northeastoutside');

%%
ax = gca;
ax.YScale = 'log';

%%
G = [];
d = [];
NeffRicter = sum(isfinite(mlvRichter),2);

for i = 1:lK-1
    thisAmp = allAmps(:,i);
    minNeff = 5;
    if i > 1
        kI1 = isfinite(thisAmp) & thisAmp > 1e-3 & thisAmp <= 2 & NeffRicter >= minNeff;
    else
        kI1 = isfinite(thisAmp) & thisAmp <= 60 & NeffRicter >= minNeff;
    end
    
    % version 2
    %     if i > 1
    %         kI1 = isfinite(thisAmp) & thisAmp > 1e0 & thisAmp <= 5e2 & NeffRicter >= minNeff;
    %     else
    %         kI1 = isfinite(thisAmp) & thisAmp <= 6e4 & NeffRicter >= minNeff;
    %     end
    
    %%
    for j = i+1:lK
        nextAmp = allAmps(:,j);
        kI = kI1 & nextAmp > 1e-3 & nextAmp <= 2 & isfinite(nextAmp);
        %kI = kI1 & nextAmp > 1e0 & nextAmp <= 5e2 & isfinite(nextAmp); %version 2
        sumNext = sum(kI);
        if sumNext
            disp(['Found ',num2str(sumNext),' amp. ratios, (',char(kstnms(i)),'-',char(kstnms(j)),')']);
            Gsc = zeros(1,lK);
            Gsc(i) = 1;
            Gsc(j) = -1;
            G_ = repmat([log10(dist(j)/dist(i)) dist(j)-dist(i) Gsc],sumNext,1);
            G = [G; G_];
            d_ = log10(thisAmp(kI)./nextAmp(kI));
            d = [d; d_];
        else
            disp(' ');
            disp(['Couldnt find matches for pair: ',char(kstnms(i)),'-',char(kstnms(j))]);
            disp(' ');
        end
    end
end

%%
G = [G; 0 0 ones(1,lK)];
d = [d; 0];

%%
m1 = ((G'*G)^-1)*G'*d;
disp(m1);

%%
m = pinv(G)*d;
disp(m);

%%
rsynth = (1:1000)';
logA0_richter = 1.11*log10(rsynth) + 0.00189*rsynth + 0.591;
logA01 = m1(1)*log10(rsynth) + m1(2)*rsynth;
logA0 = m(1)*log10(rsynth) + m(2)*rsynth;

%%
fix17 = true;
if fix17
    rI = find(rsynth == 17);
    gamma1 = 2 - logA01(rI);
    gamma = 2 - logA0(rI);
else
    rI = find(rsynth == 100);
    gamma1 = 3 - logA01(rI);
    gamma = 3 - logA0(rI);
end

%%
figure();
semilogx(rsynth,logA0_richter,'.-','DisplayName','richter (california)');
zoom on;
grid on;
hold on;
semilogx(rsynth,logA01+gamma1,'s-','DisplayName','generalized inverse');
semilogx(rsynth,logA0+gamma,'o-','DisplayName','min. norm. (this study)');
legend('show','location','northwest');

%% gen new ml (generalized inverse)
newmlv1 = NaN(size(allAmps));
for i = 1:lK
    a = allAmps(:,i);
    newmlv1(:,i) = log10(a) + m1(1)*log10(dist(i)) + m1(2)*dist(i) + gamma1 - m1(i+2);
end
Neff1 = sum(isfinite(newmlv1),2);

figure('units','normalized','outerposition',[0 0 1 1]);
ax_(1) = subplot(2,2,1);
plot(refs,newmlv1,'.'); zoom on; grid on; hold on;
title('individual magnitudes, generalized inverse')

%%
errorPerEvent1 = nanstd(newmlv1,0,2);
%errorPerEvent1 = mad(newmlv1,1,2);
ax_(2) = subplot(2,2,3);
plot(refs,errorPerEvent1,'o'); zoom on; grid on;
title('error for generalized inverse');
sum(errorPerEvent1);

%% gen new ml
newmlv = NaN(size(allAmps));
for i = 1:lK
    a = allAmps(:,i);
    newmlv(:,i) = log10(a) + m(1)*log10(dist(i)) + m(2)*dist(i) + gamma - m(i+2);
end
Neff = sum(isfinite(newmlv),2);

ax_(3) = subplot(2,2,2);
plot(refs,newmlv,'.'); zoom on; grid on; hold on;
title('individual magnitudes, min. norm. inversion');

%%
errorPerEvent = nanstd(newmlv,0,2);
Mlv = nanmedian(newmlv,2);
%errorPerEvent = mad(newmlv,1,2);
ax_(4) = subplot(2,2,4);
plot(refs,errorPerEvent,'o'); zoom on; grid on;
title('error for min. norm. inversion');
sum(errorPerEvent)

%%
linkaxes(ax_([1,3]),'y');
% linkaxes(ax_([2,4]),'y');
linkaxes(ax_,'x');
linkaxes(ax_([1,3]),'y');

%%
figure();
%pp = errorbar(datenum(refs(Mlv <= 2.75)),Mlv(Mlv <= 2.75),errorPerEvent(Mlv <= 2.75),'o');
pp = errorbar(datenum(refs),Mlv,errorPerEvent,'o');
zoom on; grid on;
datetick('x');

%figure(); plot(errorPerEvent(Neff>=5),newmlv(Neff>=5),'.'); zoom on; grid on;

%% this section applies mag parameters to WA amplitudes that are converted from raw counts data
load('raw2wa_sangayRegionalSensors_2','b','sagab');
allAmps = [z2pWA_SAGA_raw z2pWA_BPAT_raw z2pWA_BMAS_raw z2pWA_BULB_raw...
    z2pWA_BRUN_raw z2pWA_PUYO_raw z2pWA_TAMH_raw z2pWA_PORT_raw...
    z2pWA_PKYU_raw z2pWA_TAIS_raw];

%% convert to wood-anderson
for i = 1:lK
    if i == 1
        raw = allAmps(:,i);
        rI = raw >= 0 & isfinite(raw) & refs <= datetime(2018,11,11);
        raw(rI) = 10.^(sagab(1,2) + sagab(2,2)*log10(raw(rI)*4));
        
        allAmps(rI,i) = raw(rI);
        %raw(~rI) = 10.^(sagab(1,2) + sagab(2,2)*log10(raw(~rI)));
        raw(~rI) = 10.^(sagab(1,1) + sagab(2,1)*log10(raw(~rI)));       % version 2
        allAmps(~rI,i) = raw(~rI); % WA here
    else
        raw = allAmps(:,i);
        if strcmp("BULB",kstnms(i))
            raw = allAmps(:,i);
            bI = refs >= datetime(2016,01,32) & refs <= datetime(2017,01,91);
            raw(bI) = raw(bI)/4;
            allAmps(:,i) = raw;
        end
        
        rI = raw >= 0 & isfinite(raw);
        raw(rI) = 10.^(b(1,i-1) + b(2,i-1)*log10(raw(rI)));
        allAmps(:,i) = raw; % WA here
    end
end

%%
newmlv2 = NaN(size(allAmps));
for i = 1:lK
    a = allAmps(:,i); %PSEUDO-WOODANDERSON
    aI = a > 0 & isfinite(a); %<-- this step is important!
    newmlv2(aI,i) = log10(a(aI)) + m(1)*log10(dist(i)) + m(2)*dist(i) + gamma - m(i+2);
end

%% fig. 5
figure();
plot(refs,newmlv2,'.');
zoom on; grid on;
title('all individual mag. estimates')

%%
errorPerEvent2 = nanstd(newmlv2,0,2);
Mlv2 = nanmedian(newmlv2,2);

figure(4); hold on;
pp = errorbar(datenum(refs),Mlv2,errorPerEvent2,'s');
zoom on; grid on;
datetick('x');
legend('true WA','synthetic WA');

%%
figure();
plot(refs,Mlv-Mlv2,'o');
zoom on; grid on;
title('Difference between Min. Norm. parameters applied to true WA and Synthetic WA')

%%
figure();
subplot(211);
plot(refs,errorPerEvent,'o');
zoom on; grid on; hold on;
plot(refs,errorPerEvent2,'s');
title('Errors associated with Min. Norm. parameters applied to: true WA vs. Synthetic WA')
legend('true WA','synthetic WA','Locatiom','NorthWest');

subplot(212);
plot(refs,errorPerEvent-errorPerEvent2,'.');
title('difference in error (error($M_{trueWA}$) - error($M_{synthWA}$))');

%%
figure();
plot(refs,newmlv2(:,1),'.'); zoom on; grid on; hold on;
plot(refs,nanmedian(newmlv2(:,2:end),2),'o');
title('SAGA-only vs. Regional-only (min. norm.)');
legend('SAGA-only','Regional-only');


