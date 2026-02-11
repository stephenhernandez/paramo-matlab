%function ax = SierraNegraSumOfMoments()
clear; close all; clc;
cd '~/research/now/sierra_negra'
load('SierraNegraCatalog30Aug2019');
d = d*1e-3;

experimentFlag = true;
if experimentFlag
    azI = d <= 6.5;
    az = az(azI);
    eqlat = eqlat(azI);
    eqlon = eqlon(azI);
    eqdepth = eqdepth(azI);
    eqmag = eqmag(azI);
    t = t(azI);
    
    [az,azI] = sort(az);
    eqlat = eqlat(azI);
    eqlon = eqlon(azI);
    eqdepth = eqdepth(azI);
    eqmag = eqmag(azI);
    t = t(azI);
    
    N = 51;
    detrendFlag = false;
    az = cutWindows(az,N,N-1,detrendFlag);
    eqlat = cutWindows(eqlat,N,N-1,detrendFlag);
    eqlon = cutWindows(eqlon,N,N-1,detrendFlag);
    eqmag = cutWindows(eqmag,N,N-1,detrendFlag);
    eqdepth = cutWindows(eqdepth,N,N-1,detrendFlag);
    t = dn2dt(cutWindows(datenum(t),N,N-1,detrendFlag));
    
    for i = 1:size(eqlat,2)
        d_ = distance(eqlat(:,i),eqlon(:,i),mean(eqlat(:,i)),mean(eqlon(:,i)),refEllipse)*1e-3;
        meanD(i) = median(d_);
        stdD(i) = mad(d_,1); %mad(d_,1);
    end
    
    figure();
    S = scatter(median(eqlon),median(eqlat),2.5*exp(10*stdD),meanD,'o');
    colorbar;
    zoom on;
    hold on;
    plot(referenceLon,referenceLat,'k.');
    axis equal;
    
    plot(eqlon(:),eqlat(:),'k.');
    ax = plotEcuadorianCatalog(median(t)',median(eqlat)',median(eqlon)',median(eqdepth)',median(eqmag)',[-91.3 -91.0 -0.95 -0.7]);
end

%%
Mw = (2/3)*eqmag + 1.15;
m0 = mw2m0(Mw);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
S = scatter(eqlon,eqlat,2*exp(eqmag),az,'filled'); zoom on; colorbar;
S.MarkerEdgeColor = 'k';
S.LineWidth = 0.1;
S.MarkerFaceAlpha = 0.5;
axis equal;
hold on;
ax = gca;
plot(referenceLon,referenceLat,'k.','MarkerSize',20);
[circLat,circLon] = reckon(referenceLat,referenceLon,6500,(0:359)',refEllipse);
hold on; plot(circLon,circLat,'k.'); zoom on; axis equal;

%%
%tRefs = [datetime(2018,04,22); datetime(2018,06,26); datetime(2018,07,01); datetime(2018,08,17)];
%tRefs = [datetime(2018,04,20); datetime(2018,06,26,09,00,00); datetime(2018,06,26,20,00,00); datetime(2018,06,30); datetime(2018,08,17)];
%tRefs = [datetime(2018,04,22); datetime(2018,06,26,20,00,00); datetime(2018,06,30); datetime(2018,08,17)];
tRefs = [datetime(2018,04,22); datetime(2018,06,26,09,00,00); datetime(2018,06,26,20,00,00); datetime(2018,06,30); datetime(2018,08,17)];

magFact = 10;
azSearch = (0:24:360)';
symbols = ['o','d','v','s'];

%%
totalMoment = sum(m0(d <= 6.5));
disp(['total moment: ',num2str(totalMoment)]);
totalDays = days(tRefs(end) - tRefs(1));
meanMom = zeros(length(azSearch)-1,length(tRefs)-1);
meanAz = azSearch(1:end-1) + diff(azSearch)/2;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(azSearch)-1
    for j = 1:length(tRefs)-1
        sI = az >= azSearch(i) & az < azSearch(i+1) & d <= 6.5 & t >= tRefs(j) & t < tRefs(j+1);
        if sum(sI)
            hold on;
            meanMomentForBin = sum(m0(sI))/(days(tRefs(j+1)-tRefs(j))+1);
            S = scatter(mean(eqlon(sI)),mean(eqlat(sI)),magFact*exp(m02mw(meanMomentForBin)),mean(az(sI)),symbols(j),'filled');
            zoom on;
            colorbar;
            caxis([0 360]);
            S.MarkerEdgeColor = 'k';
            S.MarkerFaceAlpha = 0.5;
            S.LineWidth = 1;
            clear S
            meanMom(i,j) = meanMomentForBin;
        else
            disp(azSearch(i))
        end
    end
end

%%
%subplot(121); hold on; plot(circLon,circLat,'k.'); plot(referenceLon,referenceLat,'k.','MarkerSize',20); zoom on; axis equal;
%subplot(122);
c = colorbar;
caxis([0 360]);
c.Label.Interpreter = 'latex';
c.Label.String = 'azimuth';
hold on; plot(circLon,circLat,'k.'); plot(referenceLon,referenceLat,'k.','MarkerSize',20); zoom on; axis equal;
legend([datestr(tRefs(1)),' - ',datestr(tRefs(2))],...
    [datestr(tRefs(2)),' - ',datestr(tRefs(3))],...
    [datestr(tRefs(3)),' - ',datestr(tRefs(4))],...
    [datestr(tRefs(4)),' - ',datestr(tRefs(5))]);

%%
meanMomNorm = meanMom ./ meanMom(:,1);
figure('units','normalized','outerposition',[0 0 1 1]);
%subplot(211)
bar(meanAz,log10(meanMomNorm)); %,'stacked');
zoom on;
legend('1','2','3','4');
disp(nanmean(meanMomNorm));

% %%
% meanMomNorm = meanMom ./ sum(meanMom(:,1));
% meanMomNorm(:,1) = 1;
% %figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(212)
% bar(meanAz,log10(meanMomNorm)); %,'stacked');
% zoom on;
% legend('1','2','3','4');
% disp(nanmean(meanMomNorm));
