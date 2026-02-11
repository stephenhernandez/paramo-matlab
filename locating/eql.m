function [Mbest,minRMS,minIter,mag,NobsP,origModel,minRMSOrig] = eql...
    (eventIDs,ttp,dTdR,dTdZ,tts,dTsdR,dTsdZ,Xorig,Yorig,plotFlag,diasFlag) %,fixedDepth)
disp(eventIDs);

minRMSOrig = -999;
doDDFlag = false;
SFlag = true;
weightFlag = true;
refEllipse = referenceEllipsoid('wgs84');
maxIter = 9;
minTol = 1e-10;
maxDist = 500;
SWW = 5;

%%
eventIDs = string(eventIDs);
E = readSCBulletin(eventIDs,diasFlag);
origTime = E.t;
referenceLat = -1/4+round(E.lat*1e5)/1e5;
referenceLon = -1/4+round(E.lon*1e5)/1e5;
[tt,stla,stlo,stel,SFlag] = getPhaseDataAndMetaData(E,SFlag);

%% convert lon,lat to X,Y at stations and at initial model (units are in meters due to refEllipse)
[Kdists,Kaz] = distance(referenceLat,referenceLon,stla,stlo,refEllipse);
Y = Kdists.*cosd(Kaz); %longitude
X = Kdists.*sind(Kaz); %latitude

%% convert SC3 location to initial model location (x,y)
[Kdists,Kaz] = distance(referenceLat,referenceLon,E.lat,E.lon,refEllipse);
Ym0 = Kdists.*cosd(Kaz); %longitude
Xm0 = Kdists.*sind(Kaz); %latitude

%% convert from meters to km
Y = Y/1000;
X = X/1000;
Ym0 = Ym0/1000;
Xm0 = Xm0/1000;

%% filter according to distance threshold
dsP = sqrt((Xm0-X).^2 + (Ym0-Y).^2);
dI = dsP <= maxDist;
stlaF = stla(dI);
stloF = stlo(dI);
stelF = stel(dI)/1000;
XF = X(dI);
YF = Y(dI);
ttF = tt(dI,1);
ttF = seconds(ttF-origTime); %travel time since estimated origin time [sec.]
nSphases = 0;
if SFlag
    ttsF = tt(dI,2);
    ttsFI = isfinite(ttsF);
    ttsF = seconds(ttsF(ttsFI)-origTime);
    ttsFOrig = ttsF;

    stelFS = stel(ttsFI)/1000;
    XFS = XF(ttsFI);
    YFS = YF(ttsFI);
    nSphases = sum(ttsFI);
    Sweights = ones(nSphases,1);
end

if nSphases < 2
    SFlag = false;
end

%% shift reference elevation
% stelF = stelF - median(stelF);
% if SFlag
%     stelFS = stelFS - median(stelF);
% end

%% get observables
nPphases = sum(dI);
Mbest = -999;
minRMS = -999;
minIter = -999;
mag = -999;
NobsP = -999;
origModel = -999;
minRMSOrig = -999;

%%
if nPphases < 3
    fprintf('not enough data to continue. aborting.\n');
    return;
end

Pweights = ones(nPphases,1);
disp(['Number of usable P phases: ',num2str(nPphases)]);
NobsP = 0.5*nPphases*(nPphases-1);
obsDiffP = getDD(ttF);

if SFlag
    disp(['Number of usable S phases: ',num2str(nSphases)]);
    NobsS = 0.5*nSphases*(nSphases-1);
    obsDiffS = getDD(ttsF);
end

origDepth = E.depth; %-3; %depth relative to sea level
if origDepth < 0
    origDepth = 0;
end

Zm0 = origDepth;
if Zm0 > 200
    Zm0 = Zm0/10;
end

Tm0 = 0;
ttFOrig = ttF;
Mcurrent = [Xm0; Ym0; Zm0; Tm0];
numberOfParameters = length(Mcurrent);
origModel = Mcurrent;
dM = zeros(numberOfParameters,1);

Mall = zeros(numberOfParameters,maxIter);
%sigmaM = Mall;
rms = NaN(maxIter,1);
partialsP = ones(nPphases,numberOfParameters);
disp(['Number of P double differences: ',num2str(NobsP)]);
if SFlag
    partialsS = ones(nSphases,numberOfParameters);
    disp(['Number of S double differences: ',num2str(NobsS)]);
end

tol = 1;
ii = 0;

while tol > minTol && ii < maxIter
    %% do forward modeling
    ii = ii+1;

    disp(['iteration number: ',num2str(ii)]);
    disp(['tol: ',num2str(tol)]);
    Mcurrent = Mcurrent + dM;
    Mall(:,ii) = Mcurrent;
    disp(['Current Model: ', num2str(Mcurrent')]);
    Xm0 = Mcurrent(1);
    Ym0 = Mcurrent(2);
    Zm0 = Mcurrent(3);
    Tm0 = Mcurrent(4);

    dX = Xm0-XF;
    dY = Ym0-YF;
    dsP = sqrt(dX.^2 + dY.^2); %distance to stations <-- flat-earth approximation
    %Pweights = 1*dsP(1)./dsP; %Pweights./dsP; %1*dsP(1)./dsP;
    obsDiffPWeights = abs(getDD(dsP));
    obsDiffPWeights = obsDiffPWeights(1)./obsDiffPWeights;
    pddWI = obsDiffPWeights>1;
    obsDiffPWeights(pddWI) = 0.75;
    drdx = dX./dsP; %these are vectors, use in loop later
    drdy = dY./dsP;
    correctedDepthP = Zm0*ones(nPphases,1);%+stelF;
    correctedDepthVectP = correctedDepthP; %*ones(Ns,1);
    synth = interp2(Xorig,Yorig,ttp,correctedDepthVectP,dsP,'spline'); %<-- get synthetic travel time
    synth = synth+Tm0; %arrival time

    if SFlag
        dXs = Xm0-XFS;
        dYs = Ym0-YFS;
        dsS = sqrt(dXs.^2 + dYs.^2); %distance to stations
        Sweights = SWW*dsS(1)./dsS; %Sweights./dsS;
        obsDiffSWeights = abs(getDD(dsS));
        obsDiffSWeights = obsDiffSWeights(1)./obsDiffSWeights;
        %         pddWI = obsDiffSWeights>1;
        %         obsDiffSWeights(pddWI) = 0.5;

        drsdx = dXs./dsS; %these are vectors, use in loop later
        drsdy = dYs./dsS;
        correctedDepthS = Zm0*ones(nSphases,1);%+stelFS;
        correctedDepthVectS = correctedDepthS;
        synthS = interp2(Xorig,Yorig,tts,correctedDepthVectS,dsS,'spline'); %<-- get synthetic travel time
        synthS = synthS+Tm0; %arrival time
    end

    %% get synthetic data and intial model residuals
    synthDiffP = getDD(synth);
    Psd = ttF - synth; %observed arrival time minus synthetic arrival time
    Pdd = obsDiffP - synthDiffP; %

    if SFlag
        synthDiffS = getDD(synthS);
        Ssd = ttsF - synthS;
        Sdd = obsDiffS - synthDiffS;
    end

    %% get residual vector
    if doDDFlag
        if SFlag
            res = [Psd;Ssd;Pdd;Sdd];
            if weightFlag
                weights = [Pweights; Sweights; obsDiffPWeights; obsDiffSWeights];
                weights = weights/sum(weights);
            else
                weights = ones(size(res));
                weights = weights/sum(weights);
            end
        else
            res = [Psd;Pdd];
            if weightFlag
                weights = [Pweights;obsDiffPWeights];
                weights = weights/sum(weights);
            else
                weights = ones(size(res));
                weights = weights/sum(weights);
            end
        end
    else
        if SFlag
            res = [Psd;Ssd];
            if weightFlag
                weights = [Pweights;Sweights];
                weights = weights/sum(weights);
            else
                weights = ones(size(res));
                weights = weights/sum(weights);
            end
        else
            res = Psd;
            if weightFlag
                weights = Pweights;
                weights = weights/sum(weights);
            else
                weights = ones(size(res));
                weights = weights/sum(weights);
            end
        end
    end
    rms(ii) = sqrt(sum(weights.*(res.^2)));
    if ii > 1
        tol = abs(rms(ii)-rms(ii-1));
        disp(['New tol: ',num2str(tol)])
    end
    disp(['RMS: ',num2str(rms(ii))]);

    %% get partial derivatives evaluated at each station (for each phase)
    partials_dTdR_P = interp2(Xorig,Yorig,dTdR,correctedDepthVectP,dsP,'spline');
    partialsP(:,1) = partials_dTdR_P.*drdx; %dTdx
    partialsP(:,2) = partials_dTdR_P.*drdy; %dTdy
    partialsP(:,3) = interp2(Xorig,Yorig,dTdZ,correctedDepthVectP,dsP,'spline'); %dTdz
    if SFlag
        partials_dTdR_S = interp2(Xorig,Yorig,dTsdR,correctedDepthVectS,dsS,'spline');
        partialsS(:,1) = partials_dTdR_S.*drsdx; %dTdx
        partialsS(:,2) = partials_dTdR_S.*drsdy; %dTdy
        partialsS(:,3) = interp2(Xorig,Yorig,dTsdZ,correctedDepthVectS,dsS,'spline'); %dTdz
    end

    %% get G matrix
    GP = getDD(partialsP);
    if SFlag
        GS = getDD(partialsS);
    end

    if doDDFlag
        if SFlag
            G = [partialsP;partialsS;GP;GS];
        else
            G = [partialsP;GP];
        end
    else
        if SFlag
            G = [partialsP;partialsS];
        else
            G = partialsP;
        end
    end

    %% perform inversion
    %sigmaM_ = sigmaConst2*GtGprime;
    %sigmaM(:,ii) = sqrt(diag(sigmaM_));
    if weightFlag
        %dM = lscov(G,res,1*weights);
        w2 = sqrt(weights);
        G = repmat(w2,1,size(G,2)).*G;
        res = w2.*res;
    end
    %dM = pinv(G)*res;
    %dM = lsqnonneg(G,res);
    %dM = pinv(G)*res
    %dM = lsqnonneg(G,res)
    dM = lscov(G,res)
    dM = (res\G)';

    %GtGprime = (G'*G)^-1;
    %GtGprimeGT = GtGprime*G';
    %dM = GtGprimeGT*res;
    disp(['dM: ',num2str(dM')])
    disp(' ')
end
rms = rms(1:ii);
[~,minIter] = nanmin(rms);
Mbest = Mall(:,minIter);
disp('Mbest')
disp(Mbest)

Xm0 = Mbest(1);
Ym0 = Mbest(2);
Zm0 = Mbest(3);
Tm0 = Mbest(4);
vertical_off = -(Zm0 - E.depth); %-3;

if SFlag
    ttF = ttFOrig; %observed arrival time (relative to original origin time)
    ttsF = ttsFOrig;
    dX = -XF;
    dY = -YF;
    dsP = sqrt(dX.^2 + dY.^2); %distance to stations
    correctedDepthP = origDepth*ones(nPphases,1); %+stelF;
    correctedDepthVectP = correctedDepthP; %*ones(Ns,1);
    dXs = -XFS;
    dYs = -YFS;
    dsS = sqrt(dXs.^2 + dYs.^2); %distance to stations
    correctedDepthS = origDepth*ones(nSphases,1); %+stelFS;
    correctedDepthVectS = correctedDepthS;
    synth = interp2(Xorig,Yorig,ttp,correctedDepthVectP,dsP,'spline'); %synthetic arrival time (relative to original origin time)
    synthOrig = synth+0; % corrected synthetic arrival time
    synthS = interp2(Xorig,Yorig,tts,correctedDepthVectS,dsS,'spline');
    synthSOrig = synthS+0;
    origRes = [ttF - synthOrig; ttsF - synthSOrig];
    resMaxOrig = ceil(max(abs(origRes)));
    minRMSOrig = rms(1);

    ttF = ttFOrig; %observed arrival time (relative to original origin time)
    ttsF = ttsFOrig;
    dX = Xm0-XF;
    dY = Ym0-YF;
    dsP = sqrt(dX.^2 + dY.^2); %distance to stations
    correctedDepthP = Zm0*ones(nPphases,1);%+stelF;
    correctedDepthVectP = correctedDepthP; %*ones(Ns,1);
    dXs = Xm0-XFS;
    dYs = Ym0-YFS;
    dsS = sqrt(dXs.^2 + dYs.^2); %distance to stations
    correctedDepthS = Zm0*ones(nSphases,1); %+stelFS;
    correctedDepthVectS = correctedDepthS;
    synth = interp2(Xorig,Yorig,ttp,correctedDepthVectP,dsP,'spline'); %synthetic arrival time (relative to original origin time)
    synth = synth+Tm0; % corrected synthetic arrival time
    synthS = interp2(Xorig,Yorig,tts,correctedDepthVectS,dsS,'spline');
    synthS = synthS+Tm0;
    actualRes = [ttF - synth; ttsF - synthS];

    resMax = ceil(max(abs(actualRes)));
    minRMS = rms(minIter); %sqrt(mean(actualRes.^2));

    disp('minRMS')
    disp(minRMS)

    [theta,rho] = cart2pol(Mbest(1)*1000,Mbest(2)*1000);
    theta = 90 - theta*180/pi;
    [Mbest(2),Mbest(1)] = reckon(referenceLat,referenceLon,rho,theta,refEllipse);
    horiz_off = distance(E.lat,E.lon,Mbest(2),Mbest(1),refEllipse)*1e-3;

    origModel(1) = E.lon;
    origModel(2) = E.lat;

    disp('Mbest')
    disp(Mbest)

    disp('Morig')
    disp(origModel)
    mag = E.mag;

    Kdists = distance(Mbest(2),Mbest(1),stlaF,stloF,refEllipse);
    Kdists = Kdists/1000;

    if plotFlag
        load ~/igdata/soam_noec.mat
        load ~/igdata/ec_boundaries.mat
        load ~/igdata/ecTrench
        load ~/igdata/ecuador_slab_model.mat
        LineWidth = 2;

        fig(1) = figure('units','normalized','outerposition',[0 0 1 1]);
        plot(rms,'s-','linewidth',LineWidth);
        xlabel('Iteration Number');
        ylabel('RMS')
        grid on; zoom on;

        fig(2) = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2,2,[1 2]);
        plot([0 maxDist],[minRMSOrig minRMSOrig],'m--','linewidth',LineWidth); hold on;
        plot(Kdists,origRes(1:length(XF)),'d');
        plot(Kdists(ttsFI),origRes(length(XF)+1:end),'^','linewidth',LineWidth);
        legend('Residual Root Mean Square','P-wave Residual (SC3)','S-wave Residual (SC3)','Location','SouthEast');
        plot([0 maxDist],[0 0],'-','color',[0.5 0.5 0.5],'linewidth',LineWidth);
        ylabel('residual [sec.]');
        title('IASP91 Velocity Model, Original SC3 Solution, No Topography Correction')
        grid on; zoom on;

        subplot(2,2,3);
        plot(Kdists,ttF,'o','linewidth',LineWidth); hold on;
        plot(Kdists,synthOrig,'d','linewidth',LineWidth);
        legend('Observed P','SC3 Predicted P','Location','SouthEast')
        ylabel('Travel Time from Best Solution');
        xlabel('Distance from Best Solution Epicenter [km.]');
        title(['Number of used P phases: ',num2str(length(ttF))])
        grid on; zoom on;

        subplot(2,2,4);
        plot(Kdists(ttsFI),ttsF,'s','linewidth',LineWidth); hold on;
        plot(Kdists(ttsFI),synthSOrig,'d','linewidth',LineWidth);
        legend('Observed S','SC3 Predicted S','Location','SouthEast')
        ylabel('Travel Time from Best Solution');
        xlabel('Distance from Best Solution Epicenter [km.]');
        title(['Number of used S phases: ',num2str(sum(ttsFI))])
        grid on; zoom on;

        fig(3) = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2,2,[1 2]);
        plot([0 maxDist],[minRMS minRMS],'m--','linewidth',LineWidth); hold on;
        plot(Kdists,actualRes(1:length(XF)),'o','linewidth',LineWidth);
        plot(Kdists(ttsFI),actualRes(length(XF)+1:end),'s','linewidth',LineWidth);
        legend('Residual Root Mean Square','P-wave Residual (SH)','S-wave Residual (SH)','Location','SouthEast');
        plot([0 maxDist],[0 0],'-','color',[0.5 0.5 0.5]);
        ylabel('residual [sec.]');
        title('IASP91 Velocity Model, Best SH Solution, No Topography Correction')
        grid on; zoom on;
        %
        subplot(2,2,3);
        plot(Kdists,ttF,'o','linewidth',LineWidth); hold on;
        plot(Kdists,synth,'p','linewidth',LineWidth);
        legend('Observed P','SH Predicted P','Location','SouthEast')
        ylabel('Travel Time from Best Solution');
        xlabel('Distance from Best Solution Epicenter [km.]');
        title(['Number of used P phases: ',num2str(length(ttF))])
        grid on; zoom on;
        %
        subplot(2,2,4);
        plot(Kdists(ttsFI),ttsF,'s','linewidth',LineWidth); hold on;
        plot(Kdists(ttsFI),synthS,'^','linewidth',LineWidth);
        legend('Observed S','SH Predicted S','Location','SouthEast')
        ylabel('Travel Time from Best Solution');
        xlabel('Distance from Best Solution Epicenter [km.]');
        title(['Number of used S phases: ',num2str(sum(ttsFI))])
        grid on; zoom on;

        fig(4) = figure('units','normalized','outerposition',[0 0 1 1]);
        leg1 = plot(Mbest(1),Mbest(2),'p','linewidth',2,'markersize',20); hold on;
        leg2 = plot(origModel(1),origModel(2),'s','linewidth',2,'markersize',20);
        plot(lonEC,latEC,'k','linewidth',2);
        plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
        plot(lonTrench,latTrench,'k--','linewidth',2);
        geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
        grid on; zoom on;
        scatter(stloF,stlaF,100,actualRes(1:length(XF)),'filled');
        colorbar;
        caxis([-resMax resMax]);
        xlabel('Longitude'); ylabel('Latitude');
        zlevs = 1000*(10:5:40);
        ax = gca;
        for iii = 1:length(zlevs)
            Ctmp = contour(lonSlab,latSlab,slabDepth,[zlevs(iii) zlevs(iii)],'k.-');
            [xend,yend] = plotContourMatrix(ax,Ctmp);
            text(xend(1),yend(1),[num2str(zlevs(iii)/1000),' km.'],'FontSize',15);
        end
        axis equal;
        axis([-82 -77 -5 2]);
        title({['\makebox[6in][c]{SC3 $T_0$: ',datestr(origTime,'dd-mmm-yyyy HH:MM:SS.FFF'),...
            ', Depth: ',num2str(origModel(3)),' km., wRMS: ',num2str(rms(1)),' sec.}'],...
            ['\makebox[6in][c]{SH $T_0$: ',datestr(origTime+Mbest(4)/86400,'dd-mmm-yyyy HH:MM:SS.FFF'),...
            ', Depth: ',num2str(Mbest(3)),' km., wRMS: ',num2str(minRMS),' sec.}'],...
            ['\makebox[6in][c]{horizontal offset: ',num2str(horiz_off),', vertical offset: ',num2str(vertical_off),...
            ', timing offset: ',num2str(Tm0),'}']});
        legend([leg1, leg2],'SH Solution','SC3 Solution','Location','NorthWest');
        grid on; zoom on;

        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(211);
        scatter(Kdists,actualRes(1:length(XF)),100,stelF,'filled');
        colorbar; grid on; zoom on;

        subplot(212);
        scatter(actualRes(1:length(XF)),stelF-mean(stelF),'o','linewidth',LineWidth);
        ylabel('deviation from mean station elevation [km.]');
        xlabel('residual [sec.]');
        grid on; zoom on;
    end
else
    ttF = ttFOrig; % + Tm0;
    dX = Xm0-XF;
    dY = Ym0-YF;
    dsP = sqrt(dX.^2 + dY.^2); %distance to stations
    correctedDepthP = Zm0*ones(nPphases,1);%+stelF;
    correctedDepthVectP = correctedDepthP; %*ones(Ns,1);
    synth = interp2(Xorig,Yorig,ttp,correctedDepthVectP,dsP,'spline');
    synth = synth+Tm0;
    actualRes = ttF - synth;
    resMax = ceil(max(abs(actualRes)));
    minRMS = sqrt(mean(actualRes.^2));

    disp('minRMS')
    disp(minRMS)

    [theta,rho] = cart2pol(Mbest(1)*1000,Mbest(2)*1000);
    theta = 90 - theta*180/pi;
    [Mbest(2),Mbest(1)] = reckon(referenceLat,referenceLon,rho,theta,refEllipse);
    horiz_off = distance(E.lat,E.lon,Mbest(2),Mbest(1),refEllipse)*1e-3;

    origModel(1) = E.lon;
    origModel(2) = E.lat;

    disp('Mbest')
    disp(Mbest)

    disp('Morig')
    disp(origModel)
    mag = E.mag;


    if plotFlag
        load ~/igdata/soam_noec.mat
        load ~/igdata/ec_boundaries.mat
        load ~/igdata/ecTrench
        load ~/igdata/ecuador_slab_model.mat
        LineWidth = 2;
        figure('units','normalized','outerposition',[0 0 1 1]);
        leg1 = plot(Mbest(1),Mbest(2),'p','linewidth',LineWidth); hold on;
        leg2 = plot(origModel(1),origModel(2),'s','linewidth',LineWidth);
        plot(lonEC,latEC,'k','linewidth',2);
        plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
        plot(lonTrench,latTrench,'k--','linewidth',2);
        geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
        grid on; zoom on;

        scatter(stloF,stlaF,[],actualRes(1:length(XF)),'filled');
        colorbar; caxis([-resMax resMax]);
        xlabel('Lontiude'); ylabel('Latitude');
        axis equal
        axis([-82 -77 -5 2]);

        title({['\makebox[6in][c]{SC3 $T_0$: ',datestr(origTime,'dd-mmm-yyyy HH:MM:SS.FFF'),...
            ', Depth: ',num2str(origModel(3)),' km., wRMS: ',num2str(rms(1)),' sec.}'],...
            ['\makebox[6in][c]{SH $T_0$: ',datestr(origTime+Mbest(4)/86400,'dd-mmm-yyyy HH:MM:SS.FFF'),...
            ', Depth: ',num2str(Mbest(3)),' km., wRMS: ',num2str(minRMS),' sec.}'],...
            ['\makebox[6in][c]{horizontal offset: ',num2str(horiz_off),', vertical offset: ',num2str(vertical_off),...
            ', timing offset: ',num2str(Tm0),'}']});
        legend([leg1, leg2],'SH Solution','SC3 Solution','Location','NorthWest');

        zlevs = 1000*(10:5:40);
        ax = gca;
        for iii = 1:length(zlevs)
            Ctmp = contour(lonSlab,latSlab,slabDepth,[zlevs(iii) zlevs(iii)],'k.-');
            [xend,yend] = plotContourMatrix(ax,Ctmp);
            text(xend(1),yend(1),[num2str(zlevs(iii)/1000),' km.'],'FontSize',15);
        end
        
        figure(2);
        plot(rms,'s-','linewidth',LineWidth);

        figure(3);
        plot(abs(diff(rms)),'o','linewidth',LineWidth);
    end

    Kdists = distance(Mbest(2),Mbest(1),stlaF,stloF,refEllipse);
    Kdists = Kdists/1000;

    if plotFlag
        load ~/igdata/soam_noec.mat
        load ~/igdata/ec_boundaries.mat
        load ~/igdata/ecTrench
        load ~/igdata/ecuador_slab_model.mat
        LineWidth = 2;

        figure(4); plot(Kdists,actualRes,'o','linewidth',LineWidth);
        hold on; plot([0 maxDist],[0 0],'-','color',[0.5 0.5 0.5],'linewidth',LineWidth);
        plot([0 maxDist],[minRMS minRMS],'m--','linewidth',LineWidth);
        xlim([0 maxDist])
        ylabel('residual [sec.]');
        grid on; zoom on;

        %         figure(5); plot(Kdists,ttF,'o'); hold on; plot(Kdists,synth,'p');
        %         xlim([0 maxDist])

        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(211);
        plot([0 maxDist],[minRMS minRMS],'m--','linewidth',LineWidth); hold on;
        plot(Kdists,actualRes(1:length(XF)),'o','linewidth',LineWidth);
        legend('Residual Root Mean Square','P-wave Residual','Location','SouthEast');
        plot([0 maxDist],[0 0],'-','color',[0.5 0.5 0.5],'linewidth',LineWidth);
        ylabel('residual [sec.]');
        title('IASP91 Velocity Model, (No Topography Correction)')
        grid on; zoom on;

        subplot(212);
        plot(Kdists,ttF,'o','linewidth',LineWidth); hold on;
        plot(Kdists,synth,'p','linewidth',LineWidth);
        legend('Observed P','Predicted P','Location','SouthEast')
        ylabel('Travel Time from Best Solution');
        xlabel('Distance from Best Solution Epicenter [km.]');
        title(['Number of used P phases: ',num2str(length(ttF))]);
        grid on; zoom on;
    end
end
