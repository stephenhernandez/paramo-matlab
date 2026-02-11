%june17_2019
clear;
close all;
clc;

%%
[t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,nPhases,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,eqType,creationTime] = readCat1();
eqmag = round(eqmag*10)/10;
winSize = (1:1:24)';%[(1:9)'; (10:10:90)';  (100:100:400)']; %hours
refTime = datetime(2019,05,26,07,41,15);
minMag = 0;
magFact = 5;
alphaValue = 0.4;
minDist = 5;
maxMedianDist = 10;
nearestN = 40;

%%
tI = t >= refTime - days(max(winSize)) & eqmag >= minMag;
t = t(tI);
eqlat = eqlat(tI);
eqlon = eqlon(tI);
eqdepth = eqdepth(tI);
eqmag = eqmag(tI);
refEllipse = referenceEllipsoid('wgs84');
load('~/igdata/denseEcuadorGridNonInsular');

%%
load ~/igdata/ec_boundaries.mat
N1 = 98092; N2 = 234003; %these are the boundaries for the largest conterminous area of continental
lonEC = lonEC(N1:N2);
latEC = latEC(N1:N2);

tic;
[in,on] = inpolygon(eqlon,eqlat,lonEC,latEC);
toc;

dNear = qlon;
rateBefore = qlon;
rateAfter = qlon;
minD = qlon;

%%
close all;
lw = length(winSize);
for i = 1:lw
    winSize_ = days(winSize(i));
    bI = t >= refTime - winSize_ & t <= refTime & in;
    aI = t >= refTime & t <= refTime + winSize_ & in;
    aorb = aI|bI;
    
    if i == lw
        eqlat_ = eqlat(aorb);
        eqlon_ = eqlon(aorb);
        eqlonBefore = eqlon(bI);
        eqlonAfter = eqlon(aI);
        eqlatBefore = eqlat(bI);
        eqlatAfter = eqlat(aI);
        parfor j = 1:length(qlon)
            disp(j)
            d_ = distance(qlat(j),qlon(j),eqlat_,eqlon_,refEllipse)*1e-3;
            [dSort,dI] = sort(d_);
            dNear(j) = min([maxMedianDist median(dSort(1:nearestN))]);
            minD(j) = dSort(1);
            dBefore = distance(qlat(j),qlon(j),eqlatBefore,eqlonBefore,refEllipse)*1e-3;
            dAfter = distance(qlat(j),qlon(j),eqlatAfter,eqlonAfter,refEllipse)*1e-3;
            rateBefore(j) = sum(dBefore<=dNear(j));
            rateAfter(j) = sum(dAfter<=dNear(j));
        end
    end
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    ax(1) = subplot(3,2,[1 3]);
    S = scatter(eqlon(bI),eqlat(bI),magFact*exp(eqmag(bI)),eqmag(bI),'filled');
    S.MarkerFaceAlpha = alphaValue;
    S.MarkerEdgeColor = 'k';
    hold on;
    plot(lonEC,latEC,'k','linewidth',1);
    axis equal;
    caxis([minMag 5]);
    title('Before');
    zoom on;
    
    ax(2) = subplot(3,2,[2 4]);
    S = scatter(eqlon(aI),eqlat(aI),magFact*exp(eqmag(aI)),eqmag(aI),'filled');
    S.MarkerFaceAlpha = alphaValue;
    S.MarkerEdgeColor = 'k';
    hold on;
    plot(lonEC,latEC,'k','linewidth',1);
    axis equal;
    caxis([minMag 5]);
    title('After');
    zoom on;
    linkaxes([ax(1) ax(2)],'xy');
    
    ax(3) = subplot(3,2,[5 6]);
    maxN = sum(aI|bI);
    plot(t(aI|bI),1:maxN,'k.-');
    if sum(bI) > 0 && sum(aI) > 0
        title(['Events Before: ', num2str(sum(bI)),'; Events After: ',num2str(sum(aI)), '; Ratio: ',num2str(sum(aI)/sum(bI))]);
    else
        title(['Events Before: ', num2str(sum(bI)),'; Events After: ',num2str(sum(aI))]);
    end
    
    xlim([refTime-winSize_ refTime+winSize_]);
    yyaxis(ax(3),'right');
    stem(t(bI),eqmag(bI),'s-');
    hold on;
    stem(t(aI),eqmag(aI),'d-');
    yyaxis(ax(3),'left');
    hold on;
    h2 = plot([refTime refTime],[0 max(ax(3).YLim)],'k--','linewidth',2);
    zoom on;
    suptitle([num2str(winSize(i)),' day(s)']);
    pause(1);
end

%%
% close(2);
% close(3);
figure('units','normalized','outerposition',[0 0 1 1]);
axR(1) = subplot(1,2,1);
S = scatter(qlon,qlat,[],rateBefore,'.');
S.MarkerFaceAlpha = alphaValue;
hold on;
plot(lonEC,latEC,'k','linewidth',1);
axis equal;
caxis([0 40]);
title('Before');
zoom on;

axR(2) = subplot(1,2,2);
S = scatter(qlon,qlat,[],rateAfter,'.');
S.MarkerFaceAlpha = alphaValue;
hold on;
plot(lonEC,latEC,'k','linewidth',1);
axis equal;
caxis([0 40]);
title('After');
zoom on;
linkaxes([axR(1) axR(2)],'xy');

%
figure('units','normalized','outerposition',[0 0 1 1]);
minN = 4;
rI = rateAfter >= minN & rateBefore >= minN & minD <= minDist;
qlonD = qlon(rI);
qlatD = qlat(rI);
rc = rateAfter(rI)./rateBefore(rI);
ii = rc >= 1.005;
S = scatter(qlonD(ii),qlatD(ii),2*magFact*rc(ii),rc(ii),'filled');
S.MarkerFaceAlpha = alphaValue;
hold on;
plot(lonEC,latEC,'k','linewidth',1);

rc = -1./rc;
ii = rc <= -1.005;
S = scatter(qlonD(ii),qlatD(ii),2*magFact*abs(rc(ii)),rc(ii),'filled');
S.MarkerFaceAlpha = alphaValue;
hold on;
axis equal;
caxis(3*[-1 1]);
colorbar;
zoom on;
%figure(3);
%colormap jet;
