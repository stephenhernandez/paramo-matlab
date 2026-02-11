clear; close all;

%
tic;
% [eqType,evStatus,t,eqlat,eqlon,eqdepth,~,ids,stderr,azgap,nPhases,nMLv,...
%     timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
%      locMethod,earthModel,~,~,evMode,scbullmag] = readCat1();

[evDescription,evStatus,t,eqlat,eqlon,eqdepth,~,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = ...
    readCat1(datetime(2000,01,01),datetime(2030,01,01),-10);

eqmag = scbullmag;
[eqmag,eqmagSort] = sort(eqmag,'descend');
evType = evType(eqmagSort);
evDescription = evDescription(eqmagSort);
evStatus = evStatus(eqmagSort);
t = t(eqmagSort);
eqlat = eqlat(eqmagSort);
eqlon = eqlon(eqmagSort);
eqdepth = eqdepth(eqmagSort);
ids = ids(eqmagSort);
stderr = stderr(eqmagSort);
azgap = azgap(eqmagSort);
phaseTot = phaseTot(eqmagSort);
nMLv = nMLv(eqmagSort);
timerr = timerr(eqmagSort);
eqlaterr = eqlaterr(eqmagSort);
eqlonerr = eqlonerr(eqmagSort);
eqdeptherr = eqdeptherr(eqmagSort);
eqmagerr = eqmagerr(eqmagSort);
magType = lower(magType(eqmagSort));
meanMethod = meanMethod(eqmagSort);
scbullmag = scbullmag(eqmagSort);
toc;

%
totErr = sqrt(eqdeptherr.^2 + eqlonerr.^2 + eqlaterr.^2);
depthConstraint = totErr./abs(eqdepth);

%loose
% minNmlv = 2;
% minPhases = 10;
% magerrThresh = 1;
% azGapThresh = 180;
% spatialErrThresh = 5;
% timeErrThresh = 3;
% rmsThresh = 0.5;
% minDepth = 5;
% minMag= 1;
% maxMag = 7;
% minDeptherr = 0; %0;
% maxAmp = 550;

% strict
minSPhases = 6;
minNmlv = 8; %minSPhases;            %2;
minPhases = 12;
magerrThresh = 1;
azGapThresh = 330;
spatialErrThresh = 5;
timeErrThresh = 3;
rmsThresh = 1.2;
resThresh = 1; %mag res. thresh == 1 == one order of magnitude!
minDepth = 0;
minMag= 0;
maxMag = 6.5;
minDeptherr = 0.1; %0;
maxAmp = 315;
minAmp = 3e-2;
errorDepthRatio1 = 0.3;
errorDepthRatio2 = 2;
maxDepth = 30;

tI = t >= datetime(2000,01,01) & t <= datetime(2026,01,01) & nMLv >= minNmlv & ...
    phaseTot >= minPhases & eqmag <= maxMag & eqmag >= minMag & azgap <= azGapThresh & ...
    abs(eqmagerr) <= magerrThresh & strcmp(magType,"mlv") & ...
    timerr > 0.01 & timerr <= timeErrThresh & stderr <= rmsThresh & ...
    eqdepth > minDepth & eqdepth <= maxDepth & eqdeptherr >= minDeptherr & eqlonerr >= minDeptherr & eqlaterr >= minDeptherr & ...
    (totErr < spatialErrThresh | depthConstraint <= errorDepthRatio1) & depthConstraint <= errorDepthRatio2 &...
    strcmp(evType,"earthquake");% & meanMethod~="mean"; 
sum(tI)

%%
tIOrig = tI;
tI = find(tI);
toc;

%
minObsPerStationForInversion = 10;
maxGood = 2e2;
lE = length(tI);
refEllipse = referenceEllipsoid('wgs84');
dist = [];
amp = [];
residuals = [];
magMain = [];
magerr2 = NaN(lE,1);
eqmag2 = magerr2;
minDist = magerr2;
maxDist = magerr2;
nS = magerr2;
goodI = tI;
n = 0;
uniqKstnms = [];
allKstnms = [];
idMain = [];

%%
importAll = true;
if importAll
    goodCmpList = ["HHZ";"BHZ";"ENZ";"BLZ"]; %["HHZ";"BHZ";"BLZ";"ENZ";"HNZ";"SHZ"];
    load('~/igdata/ecuadorSensorDataTable10.mat','allSNCLs','kcmpnm','knetwk','kstnm');
    lia = ismember(kcmpnm,goodCmpList) & ismember(knetwk,...
        ["EC";"8G";"4B";"CM";"CN";"IU";"OP";"PE";"TU";"XE";"XF";"XX";"Y2";"Z3";"9D"]);
    goodStations = unique(strcat(knetwk(lia),kstnm(lia),kcmpnm(lia))); %IGNORE Khole
    clear allSNCLs kcmpnm knetwk kstnm
else
    goodStations = unique(["CMBBACHHZ";"CMFLO2HHZ";"ECAATCHNZ";"ECANGUHHZ";...
        "ECANTGHHZ";"ECANTMHHZ";"ECANTSHHZ";"ECARDOHHZ";"ECARNLHHZ";...
        "ECBBILBHZ";"ECBMASBHZ";"ECBMASHHZ";"ECBMORBHZ";"ECBNASBHZ";...
        "ECBNASHHZ";"ECBONIHHZ";"ECBOSCHHZ";"ECBPATBHZ";"ECBPATHHZ";"ECBREFBHZ";...
        "ECBRTUHHZ";"ECBRUNBHZ";"ECBTAMBHZ";"ECBTERHHZ";"ECBULBBHZ";...
        "ECBULBHHZ";"ECBV15HHZ";"ECBVC2BHZ";"ECCAB1HHZ";"ECCABPHHZ";...
        "ECCASCHHZ";"ECCAYAHHZ";"ECCHL1HHZ";"ECCHL2HHZ";"ECCHMAHHZ";...
        "ECCHSHHHZ";"ECCO1VHHZ";"ECCOHCBLZ";"ECCOHCHHZ";"ECCOSEHHZ";...
        "ECCOTAHHZ";"ECCUICHHZ";"ECCUSEHHZ";"ECCUSWHHZ";"ECESM1HHZ";...
        "ECFLF1HHZ";"ECGGPCHHZ";"ECGGPTHHZ";"ECGONZHHZ";"ECILLIHHZ";...
        "ECIMBAHHZ";"ECISPGHHZ";"ECISPTHHZ";"ECJIPIHHZ";"ECJSCHHHZ";...
        "ECLAMOHHZ";"ECLNGLHHZ";"ECMCRAHHZ";"ECMCRABLZ";"ECMILOHHZ";...
        "ECMONBHHZ";"ECMORRHHZ";"ECNINAHHZ";"ECPAC1HHZ";"ECPAS1HHZ";...
        "ECPCRAHHZ";"ECPIATBLZ";"ECPIATHHZ";"ECPIS1HHZ";"ECPKYUHHZ";...
        "ECPONDHHZ";"ECPORTHHZ";"ECPPLPBLZ";"ECPPLPHHZ";"ECPTGLHHZ";...
        "ECPUEMHHZ";"ECPULUHHZ";"ECPUYOHHZ";"ECREVNHHZ";"ECREVSHHZ";...
        "ECRUNAHHZ";"ECRVRDHHZ";"ECSAG1HHZ";"ECSAGAHHZ";"ECSAGOHHZ";...
        "ECSALIHHZ";"ECSLORHHZ";"ECSNLRBLZ";"ECSNLRHHZ";"ECSUSEHHZ";...
        "ECTAISHHZ";"ECTAMHHHZ";"ECTULMHHZ";"ECTUYUHHZ";"ECURCUHHZ";...
        "ECVCESHHZ";"ECYAHUHHZ";"ECZUMBHHZ";"IUOTAVBHZ";"PENIEVBHZ";...
        "PETBM0BHZ";"PESNIGBHZ"]);
    badStations = ["ECALAUHHZ";...
        "ECBIBLHHZ";...
        "ECBVC2HHZ";...
        "ECCHIBHHZ";...
        "ECCSOLHHZ";...
        "ECELARHHZ";...
        "ECFLFRHHZ";...
        "ECHSPRHHZ";...
        "ECPARUHHZ";...
        "ECPDNSHHZ";...
        "ECROQEHHZ";...
        "ECLGCBHHZ";...
        "XFHB16HHZ";...
        "XFHB18HHZ";...
        "XFHB12HHZ"];
end
badStations = ["8GEC01HHZ";...
    "8GEC02HHZ";...
    "8GEC05HHZ";...
    "8GEC07HHZ";...
    "8GEC09HHZ";...
    "8GEC19HHZ";...
    "8GEC20HHZ";...
    "ECALAUHHZ";...
    "ECANTCHHZ";...
    "ECBIBLHHZ";...
    "ECELARHHZ";...
    "ECFLFRHHZ";...
    "ECPARUHHZ";...
    "ECPDNSHHZ";...
    "ECPIKAHHZ";...
    "ECTST1HHZ";...
    "XEBOCAHHZ";...
    "XEBUCEHHZ";...
    "XECHIBHHZ";...
    "XELOLAHHZ";...
    "XEQINDHHZ";...
    "XFHB17HHZ"];
[lia,locb] = ismember(badStations,goodStations);
goodStations(locb(lia)) = [];

%%
gCount = 1;
gFlag = true;
if gFlag
    gMax = 1e6;
    Gmain = zeros(gMax,length(goodStations)+2);
    dMain = NaN(gMax,1);
end
cumG = 0;
allDepths = [];

%
potentialStations = [];
for i = 1:lE
    id_ = ids(tI(i));
    E = readSCBulletin(id_);
    nS_ = E.nSphases;
    if nS_ < minSPhases
        fprintf("skipping %s, too few S phases\n",id_);
        continue;
    end
    MLv = E.MLv;
    distsTmp = pull(MLv,'dist');
    [~,dsI] = sort(distsTmp);
    MLv = MLv(dsI);
    kstnm = pull(MLv,'stnm');
    knetwk = pull(MLv,'ntwk');
    kcmpnm = pull(MLv,'chan');
    magEstimate = pull(MLv,'value');
    % kI = kcmpnm == "HNZ" | kcmpnm == "BLZ";
    % if sum(kI)
    %     disp(kstnm(kI))
    % end
    amp_ = pull(MLv,'amp');
    res_ = pull(MLv,'res');
    %[stla,~] = metaDataFromStationList([kstnm knetwk kcmpnm]);
    %disp([kstnm knetwk kcmpnm]);

    nPerSNC = groupcounts(strcat(knetwk,kstnm,kcmpnm));
    if any(nPerSNC > 1)
        [~,b] = unique(strcat(knetwk,kstnm,kcmpnm)); %find and remove duplicates (ignoring khole)
        b = sort(b);
        keepI = false(length(strcat(knetwk,kstnm,kcmpnm)),1); %ignore khole
        keepI(b) = true;
        MLv(~keepI) = [];
        kstnm(~keepI) = [];
        knetwk(~keepI) = [];
        kcmpnm(~keepI) = [];
        amp_(~keepI) = [];
        res_(~keepI) = [];
        magEstimate(~keepI) = [];
    end

    [stla,~] = metaDataFromStationList(kstnm,knetwk,kcmpnm); %ignore khole
    stla(isnan(stla)) = [];
    sncls_ = strcat(knetwk,kstnm,kcmpnm); %ignore khole
    potentialStations = [potentialStations; sncls_];
    potentialStations = unique(potentialStations);
    try
        good = amp_ <= maxAmp & amp_ > minAmp & magEstimate >= minMag & magEstimate <= maxMag & ...
            isfinite(stla) & abs(res_) <= resThresh;
    catch
        fprintf(2,'something went wrong\n')
        continue;
    end

    if gFlag
        good = good & ismember(sncls_,goodStations);
    end
    lgood = sum(good);

    if lgood < minNmlv
        fprintf("skipping %s, not enough observations (%d/%d)\n",id_,lgood,minNmlv);
        continue;
    end
    n = n+1;

    goodI(n) = tI(i);
    good = find(good);
    if lgood > maxGood
        good = good(1:maxGood);
    end
    lgood = maxGood;
    MLv = MLv(good);
    kstnm = pull(MLv,'stnm');
    knetwk = pull(MLv,'ntwk');
    kcmpnm = pull(MLv,'chan');
    mlvs = pull(MLv,'value');
    amp_ = pull(MLv,'amp');
    res_ = pull(MLv,'res');

    lat_ = E.lat;
    lon_ = E.lon;
    eqmagerr2 = E.magerr2;
    depth_ = E.depth;
    sncls = strcat(knetwk,kstnm,kcmpnm);

    [stla,stlo,stel] = metaDataFromStationList(kstnm,knetwk,kcmpnm);
    stel = stel*1e-3;

    if ~isempty(uniqKstnms)
        [lia,locb] = ismember(sncls,uniqKstnms);
        lia = ~lia;
        if sum(lia)
            uniqKstnms = sort([uniqKstnms; sncls(lia)]);
        end
    else
        uniqKstnms = sort([uniqKstnms; sncls]);
    end
    allKstnms = [allKstnms; [kstnm knetwk kcmpnm]];
    horDist = distance(lat_,lon_,stla,stlo,refEllipse)*1e-3;

    nAmps = length(amp_);
    mag_ = median(mlvs);
    eqmag2(n) = mag_;
    magerr2(n) = mad(mlvs,1);
    nS(n) = nS_;
    r_ = sqrt((depth_+stel).^2 + horDist.^2);
    minDist(n) = min(horDist);
    maxDist(n) = max(horDist);
    residuals = [residuals; mag_ - mlvs];
    allDepths = [allDepths; repmat(depth_,nAmps,1)];
    dist = [dist; r_];
    amp = [amp; amp_];
    magMain = [magMain; repmat(mag_,nAmps,1)];
    idMain = [idMain; repmat(id_,nAmps,1)];

    %
    [~,goodSnclsI] = ismember(sncls,goodStations);
    diffCount = nAmps*(nAmps - 1)*0.5;
    if gFlag
        Gtmp = full(Gvdcc(nAmps));
        Gtmp = Gtmp(1:end-1,:);
        difflog = getDD(log10(r_));
        difflin = getDD(r_);
        difflogamp = -getDD(log10(amp_));
        Gmain(gCount:gCount+diffCount-1,1) = difflog;
        Gmain(gCount:gCount+diffCount-1,2) = difflin;
        dMain(gCount:gCount+diffCount-1) = difflogamp;
        for j = 1:nAmps
            Gmain(gCount:gCount+diffCount-1,goodSnclsI(j)+2) = Gtmp(:,j);
        end
    end
    gCount = gCount + diffCount;
    fprintf("%d/%d %d %d %d\n",i,lE,n,length(magMain),gCount-1);
end

if gFlag
    dMain = dMain(1:gCount-1);
    Gmain = Gmain(1:gCount-1,:);
end

goodI = goodI(1:n);
eqmag2 = eqmag2(1:n);
magerr2 = magerr2(1:n);
minDist = minDist(1:n);
maxDist = maxDist(1:n);
nS = nS(1:n);

%%
tI = goodI;
t2 = t(tI);
nMLv2 = nMLv(tI);
nPhases2 = phaseTot(tI);
azgap2 = azgap(tI);
totErr2 = totErr(tI);
stderr2 = stderr(tI);
Gorig = Gmain;
dOrig = dMain;
toc;

%%
close all;
dataDir = '~/igdata';
load(fullfile(dataDir,'ec_boundaries.mat'))
S = shaperead(fullfile(dataDir,'volcanes/volcanes_2.shp'));

% Figure 1
figure('units','normalized','outerposition',[0 0 1 1]);
SS = scatter(dist,magMain - log10(amp),exp(magMain),magMain,'filled'); zoom on; colorbar; clim([3 6]);
ax = gca;
ax.XScale = 'log';
SS.MarkerFaceAlpha = 0.6;
SS.MarkerEdgeColor = 'k';
SS.MarkerEdgeAlpha = 0.5;
grid on; %ylim([1 5]); xlim([8 500]);

% Figure 2
figure('units','normalized','outerposition',[0 0 1 1]);
SS = scatter(dist,log10(amp),exp(magMain),magMain,'filled'); zoom on; colorbar; clim([3 6]);
ax = gca;
ax.XScale = 'log';
SS.MarkerFaceAlpha = 0.6;
SS.MarkerEdgeColor = 'k';
SS.MarkerEdgeAlpha = 0.5;
grid on;

% Figure 3
figure('units','normalized','outerposition',[0 0 1 1]);
SS = scatter3(111.19*(eqlon(tI)-mean(eqlon(tI))),111.19*(eqlat(tI)-mean(eqlat(tI))),-eqdepth(tI),5*exp(eqmag2),(eqdepth(tI)),'filled');
Cbar = colorbar; zoom on; grid on; axis equal; zlabel('Depth');
ylabel('Northing [km]'); xlabel('Easting [km.]'); title('Original Data');
SS.MarkerFaceAlpha = 0.6;
SS.MarkerEdgeColor = 'k';
SS.MarkerEdgeAlpha = 0.5;
ax = gca; set(ax,'ColorScale','log');
colormap turbo;
hold on;
plot(111.19*(lonEC(lonEC >= -84) - mean(eqlon(tI))),111.19*(latEC(lonEC >= -84) - mean(eqlat(tI))),'k-','linewidth',1);

%
kstnms2 = [allKstnms(:,2) allKstnms(:,1) allKstnms(:,3)];
kstnms3 = unique(kstnms2,'rows');
kstnms4 = strcat(kstnms2(:,1),kstnms2(:,2),kstnms2(:,3));
[stla,stlo,stel] = metaDataFromStationList(kstnms3(:,2),kstnms3(:,1),kstnms3(:,3));

% Figure 4
figure('units','normalized','outerposition',[0 0 1 1]);
nSubplots = length(uniqKstnms);
Nrows = ceil(sqrt(nSubplots));
Ncols = ceil(nSubplots/Nrows);
ax = gobjects(Nrows,Ncols);
tiledlayout(Nrows,Ncols, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:nSubplots
    ax(Nrows,Ncols,i) = nexttile;
    lia = strcmp(uniqKstnms(i),kstnms4);
    plot(residuals(lia),'.'); zoom on;
    title(uniqKstnms(i));
end

% figure 5
[fig,boundaryBox,Zlimits] = loadBasemap([-84 -77 -6 2],'winter',false,false);
fig.Visible = 'on';
ax = gca;
axis(ax(1),'equal');
ax(2) = axes;
hold(ax(2),'on');
SS = scatter(ax(2),eqlon(tI),eqlat(tI),exp(eqmag2),(eqdepth(tI)),'filled'); zoom on; grid on; axis equal; colorbar;
SS.MarkerFaceAlpha = 0.6;
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5;
plot(ax(2),[S.X],[S.Y],'k--','linewidth',1);
set(ax(2),'ColorScale','log');
axis(ax(2),'equal');
colormap(ax(2),'turbo');
linkaxes(ax,'xy');
axis(ax,[-83 -74 -6 2]);
ax(2).Visible = 'off';

%
cd ~/igdata/
load ecTrench.mat
figure(5); hold on; plot(lonTrench,latTrench,'r-','linewidth',5);
% uncomment for usgs solutions
% load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');
% figure(5); hold on; tI = eqmag > 6 & eqdepth >= 50 & t >= datetime(1970,01,01);
%SS = scatter(eqlon(tI),eqlat(tI),200,(eqdepth(tI)),'kp'); zoom on; grid on; axis equal; colorbar;
ax = gca; set(ax,'ColorScale','log');
figure(5); hold on; axis equal;
xlim([-83 -75]);
ylim([-5 2]);
ylabel('Latitude');
xlabel('Longitude');
clim([50 170])
load ecSlabModelUSGS.mat
figure(5); hold on; zlevs = 1000*(50:10:170);
ax = gca;
for iii = 1:length(zlevs)
    Ctmp = contour(xq,yq,-1000*newDepth,[zlevs(iii) zlevs(iii)],'k.-');
    [xend,yend] = plotContourMatrix(ax,Ctmp);
    text(xend(1),yend(1),[num2str(zlevs(iii)/1000),' km.'],'FontSize',15);
end
fig.OuterPosition = [0 0 1 1];

% Figure 6
figure();
[gcAll,gcI] = sort(groupcounts(kstnms2),'descend');
plot(gcAll,'.'); zoom on; grid on; hold on;
text((1:length(gcAll))',gcAll,uniqKstnms(gcI,1),'FontSize',10)

% Figure 7
figure('units','normalized','outerposition',[0 0 1 1]);
hold on; plot(lonEC,latEC,'k-','linewidth',2);
SS = scatter(stlo,stla,40*exp(log10(groupcounts(kstnms2))),(groupcounts(kstnms2)),'filled'); zoom on; grid on; axis equal; colorbar;
SS.MarkerFaceAlpha = 0.6;
hold on; plot(lonEC,latEC,'k-','linewidth',2);
text(stlo,stla,kstnms3(:,2),'FontSize',12); axis equal;
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5;
cmap = parula(512);
colormap((cmap));
ax = gca; set(ax,'ColorScale','log');
xlim([-83 -75]);
ylim([-5.5 2.5]);

%
if gFlag
    kstnms5 = unique(kstnms4);
    inversionI = groupcounts(kstnms4) >= minObsPerStationForInversion;
    lia = ismember(goodStations,kstnms5(inversionI));
    badI = find(~lia)+2;
    goodI = true(length(lia)+2,1);
    goodI(badI) = false;
    dMain = dOrig;
    Gmain = Gorig(:,goodI); %(:,badI) = []; %delete columsn where not enough observations collected...

    %
    badI = sum(sum(abs(Gmain(:,3:end)),2)==0);
    if sum(badI)
        Gmain(badI,:) = [];
        dMain(badI) = [];
    end
    Gmain = [Gmain; [0 0 ones(1,sum(inversionI))]];
    dMain = [dMain; 0];
    solution = Gmain\dMain;
    rdum = logspace(log10(1),log10(2000),701);
    att1 = solution(1)*log10(rdum) + solution(2)*rdum;  %hernandez
    att2 = 1.11*log10(rdum) + 0.00189*rdum + 0.591;     %uhrhammer

    % Figure 8
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on; plot(lonEC,latEC,'k-','linewidth',2);
    SS = scatter(stlo(inversionI),stla(inversionI),150*exp(abs(solution(3:end))),solution(3:end),'filled'); zoom on; grid on; axis equal; colorbar;
    SS.MarkerFaceAlpha = 0.7;
    hold on; plot(lonEC,latEC,'k-','linewidth',2);
    text(stlo,stla,kstnms3(:,2),'FontSize',12); axis equal;
    SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5;
    cmap = parula(512);
    colormap((cmap));
    clim([-0.5 0.5]);
    xlim([-83 -75]);
    ylim([-5.5 2.5]);

    % Figure 9
    figure('units','normalized','outerposition',[0 0 1 1]);
    semilogx(rdum,att1,'.'); zoom on; grid on; hold on;
    semilogx(rdum,att2,'.');
    ylim([0 6]);

    gamma17 = -interp1(rdum,att1,17) + 2;
    gamma100 = -interp1(rdum,att1,100) + 3;

    attGamma17 = att1 + gamma17;
    hold on;
    semilogx(rdum,attGamma17,'.','color',[0.5 0.5 0.5]); grid on;

    attGamma100 = att1 + gamma100;
    hold on;
    semilogx(rdum,attGamma100,'.','color',[0.3 0.3 0.3]); grid on;

    legend('new','uhrhammer','fix17','fix100','location','northwest');

    figure(1); hold on;
    semilogx(rdum,att2,'k.');
    semilogx(rdum,attGamma100,'.','color',[0.5 0.5 0.5]); grid on;

    M100tmp = log10(amp) + solution(1)*log10(dist) + solution(2)*dist + gamma100;
    M17tmp = log10(amp) + solution(1)*log10(dist) + solution(2)*dist + gamma17;
    M100 = eqmag2;
    M17 = M100;
    M2 = magMain;
    R2 = residuals;
    A2 = amp;
    solutionKstnms = kstnms5(inversionI);
    for i = 1:length(tI)
        idTmp = ids(tI(i));
        lia = ismember(idMain,idTmp);
        [lia2,goodSnclsI2] = ismember(kstnms4(lia),solutionKstnms); %goodStations);
        stationCorrections = zeros(sum(lia),1);
        stationCorrections(lia2) = solution(goodSnclsI2(lia2)+2);
        M100tmp(lia) = M100tmp(lia) + stationCorrections;
        M100(i) = mean(M100tmp(lia));

        M17tmp(lia) = M17tmp(lia)+stationCorrections;
        M17(i) = mean(M17tmp(lia));
        M2(lia) = M100(i);
        A2(lia) = log10(amp(lia)) + stationCorrections;
        R2(lia) = M100(i) - M100tmp(lia);
    end

    % Figure 10
    figure('units','normalized','outerposition',[0 0 1 1]);
    SS = scatter(dist,M2-A2,exp(M2),(allDepths),'filled'); zoom on; grid on; colorbar; %caxis([2 6]);
    SS.MarkerFaceAlpha = 0.6;
    SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5;
    ax = gca; set(ax,'ColorScale','log');
    ax.XScale = 'log';

    figure(10); hold on; semilogx(rdum,attGamma100,'-','color',[0.5 0.5 0.5],'linewidth',4); grid on;
    ylim([-1 7])

    % Figure 11
    figure('units','normalized','outerposition',[0 0 1 1]);
    ax11(1) = subplot(121);
    plot(R2,'.'); zoom on;
    xlabel('residuals');
    ax11(2) = subplot(122);
    histogram(R2,1001); zoom on; grid on;
    xlabel('residual');

    % Figure 12
    figure('units','normalized','outerposition',[0 0 1 1]);
    semilogx(dist,R2,'.'); zoom on; grid on;
    xlabel('dist');
    ylabel('residuals');

    % Figure 13
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(sort(R2),(0:length(R2)-1)'/length(R2),'.'); zoom on; grid on; hold on;
    plot(sort(residuals),(0:length(residuals)-1)'/length(residuals),'.');
    legend('New Residuals','SC3 Residuals','Location','NorthWest');

    % Figure 14
    figure('units','normalized','outerposition',[0 0 0.8 1]);
    SS = scatter(eqmag2,M100,100*exp(abs(eqmag2-M100)),(eqdepth(tI)),'filled'); zoom on; grid on;
    cbar1 = colorbar; %caxis([2 6]);
    SS.MarkerFaceAlpha = 0.6;
    SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5;
    ax = gca; set(ax,'ColorScale','log');
    hold on; plot([1.1 6.9],[1.1 6.9],'k--','linewidth',3)
    xlabel('SC5 Magnitude'); ylabel('This Study');
    axis square;
    cbar1.Label.String = 'depth [km]';

    % Figure 15
    figure('units','normalized','outerposition',[0 0 1 1]);
    nSubplots = sum(inversionI);
    Nrows = ceil(sqrt(nSubplots));
    Ncols = ceil(nSubplots/Nrows);
    ax = gobjects(Nrows,Ncols);
    tiledlayout(Nrows,Ncols, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:nSubplots
        ax(Nrows,Ncols,i) = nexttile;
        thisKstnm_ = solutionKstnms(i);
        lia = strcmp(thisKstnm_,kstnms4);
        plot(residuals(lia)-solution(2+i),'.'); zoom on;
        title(thisKstnm_);
    end

    % Figure 16
    figure('units','normalized','outerposition',[0 0 1 1]);
    SS = scatter3(111.19*(eqlon(tI)-mean(eqlon(tI))),111.19*(eqlat(tI)-mean(eqlat(tI))),-eqdepth(tI),5*exp(M100),M100-eqmag2,'filled');
    Cbar = colorbar; zoom on; grid on; axis equal; zlabel('Depth');
    ylabel('Northing [km]'); xlabel('Easting [km.]'); title('Original Data');
    SS.MarkerFaceAlpha = 0.6; %ax = gca; set(ax,'ColorScale','log');
    clim([-0.5 0.5]);
    colormap turbo
    hold on;
    plot(111.19*(lonEC(lonEC >= -84) - mean(eqlon(tI))),111.19*(latEC(lonEC >= -84) - mean(eqlat(tI))),'k-','linewidth',1);

    % Figure 17
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(sort(M100-eqmag2),(0:length(M100)-1)'/length(M100),'.');
    zoom on; grid on;

    % Figure 18
    figure('units','normalized','outerposition',[0 0 1 1]);
    semilogx(eqdepth(tI),M100-eqmag2,'.'); zoom on; grid on;
    xlabel('depth [km]'); ylabel('mag. residual');

    % Figure 19
    figure('units','normalized','outerposition',[0 0 1 1]);
    semilogx(stel(inversionI),solution(3:end),'.');
    zoom on; grid on;
    xlabel('station elevation [masl]');
    ylabel('station correction');
end
