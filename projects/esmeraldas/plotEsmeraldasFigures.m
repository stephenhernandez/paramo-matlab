clear; close all; clc;
cd ~/research/now/esmeraldas/hypodd/
!\rm -rf hypoDD.inp; \rm -rf ph2dt_esmeraldas.inp;
!ln -s hypoDD_M3d_2.inp hypoDD.inp
%!ln -s hypoDD_MSERGIO.inp hypoDD.inp
%!ln -s hypoDD_SERGIO.inp hypoDD.inp
%!ln -s hypoDD_CENTRAL.inp hypoDD.inp
%!ln -s hypoDD_MACQUET.2.inp hypoDD.inp
%!ln -s hypoDD_MFONT.inp hypoDD.inp
%!ln -s hypoDD_SOTO.inp hypoDD.inp
!ln -s ph2dt_sandro.inp ph2dt_esmeraldas.inp

%
regionName = "esmeraldas";
minLon = -79.9; maxLon = -79.61; minLat = 0.7; maxLat = 1.02;
%minLon = -80.1; maxLon = -79.25; minLat = 0.66; maxLat = 1.7;
%minLon = -80.15; maxLon = -79.35; minLat = 0.6; maxLat = 1.25;
tStart = datetime(2022,03,26);
tEnd = datetime(2022,05,12);
boundaryBox = [minLon; maxLon; minLat; maxLat];
depthCorrection = 0; %positive number for elevation above sea-level
diasFlag = false;
minMag = 0;
maxDepth = 50; depthRange = [21 16];
[E,t,eqlat,eqlon,eqdepth,eqmag,id,tOrig,eqlatOrig,eqlonOrig,eqdepthOrig,eqmagOrig,idOrig] = ...
    genHypoDDFiles(tStart,tEnd,regionName,boundaryBox,depthCorrection,diasFlag,minMag,maxDepth,depthRange);

figure('units','normalized','outerposition',[0 0 1 1]);
SS = scatter3(111.19*(eqlon-mean(eqlon)),111.19*(eqlat-mean(eqlat)),-eqdepth,5*exp(eqmag),-eqdepth,'filled');
Cbar = colorbar; zoom on; grid on; axis equal; zlabel('Depth'); ylabel('Northing [km]'); xlabel('Easting [km.]'); grid on;

figure('units','normalized','outerposition',[0 0 1 1]);
spax(1) = subplot(211); SS1 = scatter(spax(1),eqlat,-eqdepth,8*exp(eqmag),datenum(t),'filled'); colorbar;
zoom on; xlabel('Latitude');SS1.MarkerFaceAlpha = 0.5; SS1.MarkerEdgeColor = 'k'; SS1.MarkerEdgeAlpha = 0.25; grid on;
xlim([boundaryBox(3) boundaryBox(4)]); clim(datenum([tStart tEnd])); hold on;
spax(2) = subplot(212);
SS2 = scatter(spax(2),eqlon,-eqdepth,8*exp(eqmag),datenum(t),'filled'); colorbar;
zoom on; xlabel('Longitud'); SS2.MarkerFaceAlpha = 0.5; SS2.MarkerEdgeColor = 'k'; SS2.MarkerEdgeAlpha = 0.25; grid on;
xlim([boundaryBox(1) boundaryBox(2)]); clim(datenum([tStart tEnd])); hold on;
linkaxes(spax,'y');

%
ecuador_xsection();

%
msI = eqmag == max(eqmag);
pmsI = t > t(msI);
d_ = distance(eqlat(msI),eqlon(msI),eqlat,eqlon,refEllipse)*1e-3;
figure();
SS = scatter(seconds(t(pmsI)-t(msI)),d_(pmsI),10*exp(eqmag(pmsI)),eqdepth(msI)-eqdepth(pmsI),'o','filled'); zoom on; grid on;
c1 = colorbar; ax = gca; ax.YScale = 'log'; ax.XScale = 'log'; SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.8; SS.MarkerFaceAlpha = 0.5;
clim([-3 3]);
ylabel('Distance from Mainshock [km.]'); xlabel('Time since Mainshock [sec.]');
c1.Label.String = 'Depth Relative to Mainshock [km.]'; c1.Label.Interpreter = 'latex';
axis tight;

%
figure();
loglog(seconds(t(pmsI)-t(msI)),(1:sum(pmsI))','.'); zoom on; grid on;
id = char(string(id));
idIndex = str2double(string(id(:,5:end)));
newID = idOrig(idIndex);
noRelocateID = true(length(eqmagOrig),1);
noRelocateID(idIndex) = false;

%
cd ~/research/now/esmeraldas/
load EsmeraldasTemplates_2.mat

[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    loadRepeaterCatalog2("esmeraldas");
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    filterCatalog(tMain,ccMain,5,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain);

%tRef = datetime(2022,01,01);
tRef = datetime(2022,03,27,04,28,00);
ccI = ccMain >= 0.12 & madMain >= 8 & nUsedMain >= 9 & magMain >= -1 & tMain >= tRef & templateNumber ~= 33;
tMain(~ccI) = [];
ccMain(~ccI) = [];
ampMain(~ccI) = [];
dMag(~ccI) = [];
magMain(~ccI) = [];
templateNumber(~ccI) = [];
madMain(~ccI) = [];
nUsedMain(~ccI) = [];

%%
eqlatTS_Orig = pull(T(1,:),'evla')';
eqlonTS_Orig = pull(T(1,:),'evlo')';
eqdepthTS_Orig = pull(T(1,:),'evdp')';
idTS_Orig = pull(T(1,:),'evid',"")';

idMain = idTS_Orig(templateNumber);
[lia,locb] = ismember(idMain,newID);
idMain = idMain(lia);
latMain = eqlat(locb(lia));
lonMain = eqlon(locb(lia));
depthMain = eqdepth(locb(lia));
magMain2 = eqmag(locb(lia));

%%
ccI = lia;
tMain(~ccI) = [];
ccMain(~ccI) = [];
ampMain(~ccI) = [];
dMag(~ccI) = [];
magMain(~ccI) = [];
templateNumber(~ccI) = [];
madMain(~ccI) = [];
nUsedMain(~ccI) = [];
magMain3 = magMain2+dMag;

figure(); plot(tMain,1:length(tMain),'.'); zoom on; grid on;
figure(); semilogy(tMain,ampMain,'.'); zoom on; grid on;
figure(); plot(tMain,templateNumber,'.'); zoom on; grid on;
figure(); loglog(ccMain,madMain,'.'); zoom on; grid on;
figure(); plot(tMain,dMag,'.'); zoom on;
figure(); plot(tMain,[magMain magMain3],'.'); zoom on; grid on;

%
d2_ = distance(latMain(2:end),lonMain(2:end),latMain(1),lonMain(1),refEllipse)*1e-3;
figure();
scatter(seconds(tMain(2:end)-tMain(1)),d2_,'k.'); hold on;
% S1 = scatter(seconds(tMain(2:end)-tMain(1)),d2_,[],[1 1 1],'o'); zoom on; hold on;
% S1.MarkerFaceColor = 'k';
% S1.MarkerEdgeColor = [0.5 0.5 0.5];
% S1.MarkerEdgeAlpha = 0.5;

msI = eqmag == max(eqmag);
pmsI = t > t(msI);
d_ = distance(eqlat(msI),eqlon(msI),eqlat,eqlon,refEllipse)*1e-3;

SS = scatter(seconds(t(pmsI)-t(msI)),d_(pmsI),10*exp(eqmag(pmsI)),eqdepth(msI)-eqdepth(pmsI),'o','filled'); zoom on; grid on;
c1 = colorbar; ax = gca; ax.YScale = 'log'; ax.XScale = 'log'; SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.8; SS.MarkerFaceAlpha = 0.5;
clim([-3 3]);
ylabel('Distance from Mainshock [km.]'); xlabel('Time since Mainshock [sec.]');
c1.Label.String = 'Depth Relative to Mainshock [km.]'; c1.Label.Interpreter = 'latex';
axis tight;

%
figure();
loglog(seconds(tMain-min(tMain)),(1:length(tMain))','.'); zoom on; grid on;
table(unique(idMain),groupcounts(idMain))
%sum(groupcounts(idOrig(templateNumber(ccI))))


% omori stuff
tI2 = tMain > tMain(1) & tMain <= tMain(564);
t2 = tMain(tI2);
t2 = seconds(t2 - tMain(1));

%figure(); plot(t2,1:length(t2),'.'); zoom on;
% mc = 1.8;
% tAfter = t(t > t(3) & eqmag >= mc);
% length(tAfter)
% tRef = t(3);
% tAfter = seconds(tAfter - tRef);
% tAfter(1:3)
% magAfter = eqmag(t > t(3) & eqmag >= mc);

S = 0;
TE = max(t2);
N = length(t2);
c_ = 1;
p_ = 3.5;
A_ = Acp(c_,p_,S,TE);
D_ = Dcp(c_,p_,S,TE);
K_ = N/A_;
dLdc = @(c) -p_*sum(1/(t2+c)) - K_*((TE+c)^-p_ - (S+c)^-p_);

%%
lastSection = ~false;
if lastSection
    close all;
    % [tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    %     loadRepeaterCatalog2("Esmeraldas");
    [tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
        loadRepeaterCatalog2("esmeraldas");
    [tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
        filterCatalog(tMain,ccMain,3,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain);

    fig = figure();
    SS = scatter(nUsedMain.*madMain,ccMain,50,nUsedMain,'.'); zoom on; grid on; colorbar; ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
    colormap hot
    ax = gca; ax.Color = 0.8*[1 1 1];

    clear medMAD medCC;
    uniqN = unique(nUsedMain);
    n = 0;
    for i = uniqN'
        n = n + 1;
        ii = nUsedMain == i;
        medMAD(n,1) = median(i*madMain(ii));
        medCC(n,1) = median(ccMain(ii));
        figure(fig);
        hold on;
        plot(medMAD(n,1),medCC(n,1),'ko');
    end
    medCC
    b_ = robustfit(log10(medMAD),log10(medCC));
    b_ = flipud(b_)
    yfit = polyval(b_,log10(logspace(0,4,201))');
    figure(fig); hold on; plot(logspace(0,4,201)',10.^yfit,'.');
    b_2 = b_; b_2(2) = -0; yfit2 = polyval(b_2,log10(logspace(0,4,201))');

    figure(fig); hold on;
    plot(logspace(0,4,201)',10.^yfit2,'.');

    ccI = ccMain >= 10.^(polyval(b_2,log10(nUsedMain.*madMain))) & tMain < datetime(2022,05,13) & templateNumber ~= 33;
    clear ax; figure(); ax(1) = subplot(311); semilogy(tMain(ccI),ccMain(ccI),'.'); zoom on; grid on;
    hold on; semilogy(tMain(~ccI),ccMain(~ccI),'.'); zoom on; grid on; ax(2) = subplot(312); semilogy(tMain(ccI),madMain(ccI),'.'); zoom on; grid on; hold on; semilogy(tMain(~ccI),madMain(~ccI),'.'); ax(3) = subplot(313); semilogy(tMain(ccI),ampMain(ccI),'.'); zoom on; grid on; hold on; semilogy(tMain(~ccI),ampMain(~ccI),'.'); linkaxes(ax,'x');

    figure(); plot(tMain(ccI),nUsedMain(ccI),'.'); zoom on;

    fig = figure(); ax(1) = subplot(211); plot(tMain(ccI),1:sum(ccI),'.'); zoom on; grid on; ylabel("cumulative number");
    ax(2) = subplot(212); semilogy(tMain(ccI),ampMain(ccI),'.'); zoom on; grid on; ylabel("amplitude"); linkaxes(ax,'x');
    figure(fig); pause(2);

    figure();
    axM(1) = subplot(211);
    plot(tMain(ccI),dMag(ccI),'.'); zoom on; grid on;
    axM(2) = subplot(212);
    semilogy(sort(dMag(ccI)),1-((0:sum(ccI)-1)'/sum(ccI)),'.'); zoom on; grid on;

    figure();
    tMainDum = tMain(ccI);
    semilogy(tMainDum(1:end-1),seconds(diff(tMainDum)),'.'); zoom on; grid on;

    %
    %tRef = datetime(2022,01,01);
    ccI = ccI & tMain >= tRef;
    tMain(~ccI) = [];
    ccMain(~ccI) = [];
    ampMain(~ccI) = [];
    dMag(~ccI) = [];
    magMain(~ccI) = [];
    templateNumber(~ccI) = [];
    madMain(~ccI) = [];
    nUsedMain(~ccI) = [];

    eqlatTS_Orig = pull(T(1,:),'evla')';
    eqlonTS_Orig = pull(T(1,:),'evlo')';
    eqdepthTS_Orig = pull(T(1,:),'evdp')';
    idTS_Orig = pull(T(1,:),'evid',"")';

    idMain = idTS_Orig(templateNumber);
    [lia,locb] = ismember(idMain,idOrig); %newID);
    idMain = idMain(lia);
    latMain = eqlatOrig(locb(lia));
    lonMain = eqlonOrig(locb(lia));
    depthMain = eqdepthOrig(locb(lia));
    magMain2 = eqmagOrig(locb(lia));

    %
    ccIOrig = ccI;
    ccI = lia;
    tMain(~ccI) = [];
    ccMain(~ccI) = [];
    ampMain(~ccI) = [];
    dMag(~ccI) = [];
    magMain(~ccI) = [];
    templateNumber(~ccI) = [];
    madMain(~ccI) = [];
    nUsedMain(~ccI) = [];
    magMain3 = magMain2+dMag;

    figure(); plot(tMain,1:length(tMain),'.'); zoom on; grid on;
    figure(); semilogy(tMain,ampMain,'.'); zoom on; grid on;
    figure(); plot(tMain,templateNumber,'.'); zoom on; grid on;
    figure(); loglog(ccMain,madMain,'.'); zoom on; grid on;
    figure(); plot(tMain,dMag,'.'); zoom on;
    figure(); plot(tMain,[magMain magMain3],'.'); zoom on; grid on;

    dM = 0.01;
    magMain4 = round(magMain3/dM)*dM;
    Nb = 250;
    [magCut,startIndex,endIndex,badFlag,nwindows] = cutWindows(magMain4,Nb,Nb-1,false);
    clear b;
    for i = 1:size(magCut,2)
        mag_ = magCut(:,i);
        magDiff = diff(mag_);
        mPositive = magDiff >= dM*2;
        mag2 = magDiff(mPositive);
        b(i,1) = log10(exp(1))./(mean(mag2)-dM+(dM/2));
    end

    figure();
    plot(tMain(startIndex),b,'.'); zoom on; grid on;
end
