clear; close all; %clc;

centralLat = -2.00504;
centralLon = -78.34124;
westLat = -2.00481;
westLon = -78.34327;
northLat = -2.00194;
northLon = -78.34257;

dayStart = datetime(2022,10,13);
dayEnd = datetime(2022,10,13);

%dayEnd = datetime(2022,08,28);
% dayStart = datetime(2022,08,27);
% dayEnd = datetime(2022,08,28);

dInc = 1;
dayVec = (dayStart:dInc:dayEnd)';
lDays = length(dayVec);

thresh = 0.9;
initialThresh = 0.6;
newFs = 50;
Fs2 = 1e3;
secDur = 12;                                % 12 seconds
nOverlap = 3/4;                             %  9 seconds
winlen = round(secDur*newFs);
stride = round(winlen*(1-nOverlap))/newFs;  %  4 seconds

lfc = 1/2;
hfc = 4;
sensitivity = 1/3600;
nStrides = 2;
ampThresh = 1e-2;

%
refEllipse = referenceEllipsoid('wgs84');
elementLats = [-2.193637; -2.193408; -2.193701; -2.193944; -2.193641];
elementLons = [-78.099766; -78.099465; -78.099223; -78.099517; -78.099506];
elementElevs = [1236; 1236; 1234; 1233; 1238];
nComps = length(elementLats);

maxN = 1e4;
badComponents = [];
shiftsMain = NaN(nComps-length(badComponents),maxN);
meanCCsMain = NaN(maxN,1);
medCCsMain = meanCCsMain;
ampMain = meanCCsMain;
medratioMain = meanCCsMain;
stdratioMain = meanCCsMain;
velMain = meanCCsMain;
bazMain = meanCCsMain;
belMain = meanCCsMain;
tMain = NaT(maxN,1);

elementLats(badComponents) = [];
elementLons(badComponents) = [];
elementElevs(badComponents) = [];

lS = length(elementLats);
waveformMain = NaN(winlen*lS,maxN);
dx = [];
dy = [];
dz = [];
for i = 1:lS-1
    stla = elementLats(i);
    stlo = elementLons(i);
    [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
    dx = [dx; d_.*cosd(90-az)];
    dy = [dy; d_.*sind(90-az)];
    dz = [dz; elementElevs(i+1:end) - elementElevs(i)];

    %
    disp([elementLats(i),elementLons(i)])
    disp([d_ az])
end
Gorig = [dx dy]; % dz];
[Grows,Gcolumns] = size(Gorig);
plagsMain = NaN(Grows,maxN);
slownessMain = NaN(Gcolumns,maxN);
kholelist = ["01";"02";"03";"04";"05"];
minNumStations = 3;

%
n = 0;
for jj = 1:lDays
    elementLats = [-2.193637; -2.193408; -2.193701; -2.193944; -2.193641];
    elementLons = [-78.099766; -78.099465; -78.099223; -78.099517; -78.099506];
    elementElevs = [1236; 1236; 1234; 1233; 1238];
    day_ = dayVec(jj);

    fprintf("loading data for %s\n",datestr(day_));
    S = populateWaveforms(nComps);
    for i = 1:nComps
        khole_ = kholelist(i);
        S_ = loadWaveforms(day_,dInc,"SAG1","BDF","EC",khole_); %this order cannot change
        S(i) = S_;
    end

    refs = pull(S,'ref');
    badComponents = isnat(refs);
    goodComponents = ~badComponents;
    lS = sum(goodComponents);
    if lS < minNumStations
        fprintf(2,'no data for: %s\n',datestr(day_));
        continue;
    end

    elementLats(badComponents) = [];
    elementLons(badComponents) = [];
    elementElevs(badComponents) = [];

    dx = [];
    dy = [];
    dz = [];
    for i = 1:lS-1
        stla = elementLats(i);
        stlo = elementLons(i);
        [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
        dx = [dx; d_.*cosd(90-az)];
        dy = [dy; d_.*sind(90-az)];
    end
    Gorig = [dx dy];
    S(badComponents) = [];

    %
    tw = 0.004;
    fprintf(1,'processing: %s\n',datestr(day_));
    Sf = detrendWaveforms(...
        intWaveforms(...
        detrendWaveforms(...
        filterWaveforms(...
        taperWaveforms(...
        detrendWaveforms(...
        differentiateWaveforms(S)),tw),lfc,hfc))));

    fprintf('done filtering...\n');
    Sf = syncWaveforms(detrendWaveforms(resampleWaveforms(Sf,newFs)));
    %Sf = cutWaveforms(Sf,dateshift(Sf(1).ref,'start','day')+hours(23),0,hours(1));
    dstack = [];

    Fs = round(1./median(pull(Sf,'delta')));
    if Fs ~= newFs
        fprintf('something went wrong with resampling step\n');
        continue;
    end

    for i = 1:lS
        d = double(pull(Sf(i)));
        [dcut,startIndex] = cutWindows(d,winlen,nOverlap,true);
        dcut = taper(dcut,0.2);
        dstack = [dstack; dcut];
    end

    allAmps = sensitivity*0.5*peak2peak(dstack)';
    winI = allAmps >= ampThresh;

    t_ = getTimeVec(Sf(1));
    tdumcut = t_(startIndex);
    clear t_;
    nWindows = sum(winI); %size(dstack,2);
    if ~nWindows
        continue;
    end
    allAmps(~winI) = [];
    dstack(:,~winI) = [];
    si = 1 + winlen*(0:lS-1)';
    ei = winlen*(1:lS)';

    shifts = NaN(lS,nWindows);
    meanCCs = NaN(nWindows,1);
    medCCs = meanCCs;
    medratio = meanCCs;
    stdratio = meanCCs;

    %
    dispIndex = [];
    prcWins = round(100*(1:nWindows)/nWindows);
    for i = [1 10 20 30 40 50 60 70 80 90]
        dispIndex = [dispIndex; find(prcWins>i,1)];
    end

    %
    Ndiff = 0.5*lS*(lS-1);
    plags = NaN(Ndiff,nWindows);
    slownii = NaN(Gcolumns,nWindows);
    baz = NaN(nWindows,1);
    vel = baz;
    if Gcolumns > 2
        bel = baz;
    end

    %
    fprintf("running through: %d potential windows\n",nWindows);
    for i = 1:nWindows
        if ismember(i,dispIndex)
            fprintf('%g\n',floor(100*i/nWindows));
        end

        %
        G = Gorig;
        dslice = dstack(:,i);
%         ampTmp = allAmps(i);
%         if ampTmp < ampThresh
%             %fprintf("amp too small, skipping...\n");
%             continue;
%         end

        dslice = reshape(dslice,[winlen,lS]);
        [~,maxccp_] = apply_vdcc(dslice,[],false,false,false);
        sqmaxccp = squareform(maxccp_);
        sqmaxccp(sqmaxccp==0) = NaN;
        indivAvgCC = median(sqmaxccp,"omitnan")';
        newGood = indivAvgCC >= initialThresh;

        if sum(newGood) < minNumStations
            %fprintf(2,"Window %d: Skipping, not enough good data\n",i);
            continue;
        end

        G1 = squareform(G(:,1));
        G2 = squareform(G(:,2));
        G1 = G1(:,newGood);
        G2 = G2(:,newGood);
        G1 = G1'; G2 = G2';
        G1 = squareform(G1(:,newGood));
        G2 = squareform(G2(:,newGood));
        G = [G1(:) G2(:)];

        dslice = detrend(dstack(:,i));
        dslice = reshape(dslice,[winlen,lS]);
        dsliceOrig = dslice;
        dslice = resample(dsliceOrig,Fs2,newFs);
        dsliceResampled = dslice(:,newGood);
        amp_ = sensitivity*0.5*median(peak2peak(dsliceResampled));
        allAmps(i) = amp_;

        [~,~,~,~,raw_shifts,meancc,medcc] = apply_vdcc(dsliceResampled,[],false,false,false);
        raw_shifts = raw_shifts-min(raw_shifts);
        shifts(newGood,i) = raw_shifts;
        dslice = dsliceOrig;

        meanCCs(i) = meancc;
        medCCs(i) = medcc;
        slow_ = lscov(G,(getDD(raw_shifts)/Fs2)); %slow_ = ((G'*G).^(-1))*G'*(getDD(raw_shifts)/Fs2);
        slownii(:,i) = slow_;

%         %plags_ = -plags_(1:end-1)/Fs2;
%         %slow_ = plags_\G;
%         %slow_ = slow_';

        if Gcolumns > 2
            [azimuth,elevation,r] = cart2sph(slow_(1),slow_(2),slow_(3));
            vel(i,1) = r; %rssq(slow_);
            baz(i,1) = 90 - rad2deg(azimuth);
            bel(i,1) = rad2deg(elevation);
            if bel(i) > 90
                bel(i) = bel(i) - 90;
            end
            fprintf('%g %g %g\n',vel(i),baz(i),bel(i));
        else
            baz(i,1) = 90 - atan2d(slow_(1),slow_(2)); %+360;
            vel(i,1) = 1./rssq(slow_);
            if baz(i) > 360
                baz(i) = baz(i)-360;
            end
        end

        if baz(i) > 180
            baz(i) = baz(i) - 360;
        end

        if baz(i) < -180
            baz(i) = baz(i) + 360;
        end

        % super dirty, must be a way to speed up
        medampratios = NaN(1,Ndiff);
        dslice = abs(dslice);
        refBlock = dslice;

        ri = 1;
        ln = size(refBlock,2)-1;
        for j = 1:lS-1
            d1 = dslice(:,j);
            refBlock = circshift(refBlock,-1,2);
            refBlock = refBlock(:,1:end-1);
            ln = size(refBlock,2);
            medampratios(ri:ri+ln-1) = median(log10(d1./refBlock),'omitnan');
            ri = ri + ln;
        end
        medratio(i) = median(medampratios,'omitnan');
        stdratio(i) = mad(medampratios,1,2); %'omitnan');
    end
    
    medCCsOrig = medCCs;
    goodWindowsI = medCCs >= thresh & allAmps >= ampThresh;
    summi = sum(goodWindowsI);
    if ~summi
        continue;
    end

    %
    t = getTimeVec(Sf);
    goodWindowsI = find(goodWindowsI);
    if summi > 1
        ccGood = medCCs(goodWindowsI);
        tGood = t(startIndex(goodWindowsI));
        tdiff = seconds(diff(tGood));
        lstride = tdiff < nStrides*stride;
        sumlstride = sum(lstride);
        while sumlstride
            lstride = find(lstride,1);
            i1 = lstride(1);
            i2 = i1 + 1;
            cc1 = ccGood(i1);
            cc2 = ccGood(i2);
            if cc1 > cc2
                %fprintf('deleting: (%s,%g), keeping (%s,%g)\n',datestr(tGood(i2)),cc2,datestr(tGood(i1)),cc1);
                goodWindowsI(i2) = [];
            else
                %fprintf('deleting: (%s,%g), keeping (%s,%g)\n',datestr(tGood(i1)),cc1,datestr(tGood(i2)),cc2);
                goodWindowsI(i1) = [];
            end
            ccGood = medCCs(goodWindowsI);
            tGood = t(startIndex(goodWindowsI));
            tdiff = seconds(diff(tGood));
            lstride = tdiff < nStrides*stride;
            sumlstride = sum(lstride);
        end
        summi = length(goodWindowsI);
    end

    %
    n = n + 1;
    allAmps = allAmps(goodWindowsI);
    dstack = dstack(:,goodWindowsI);
    shifts = shifts(:,goodWindowsI);
    meanCCs = meanCCs(goodWindowsI);
    medCCs = medCCs(goodWindowsI);
    plags = plags(:,goodWindowsI);
    vel = vel(goodWindowsI);
    baz = baz(goodWindowsI);
    waveformMain(:,n:n+summi-1) = dstack;
    if Gcolumns > 2
        bel = bel(goodWindowsI);
    end

    shiftsMain(:,n:n+summi-1) = shifts;
    ampMain(n:n+summi-1) = allAmps; %sensitivity*0.5*peak2peak(dstack)';
    meanCCsMain(n:n+summi-1) = meanCCs;
    medCCsMain(n:n+summi-1) = medCCs;
    tMain(n:n+summi-1) = t(startIndex(goodWindowsI));
    medratioMain(n:n+summi-1) = medratio(goodWindowsI);
    stdratioMain(n:n+summi-1) = stdratio(goodWindowsI);
    slownessMain(:,n:n+summi-1) = slownii(:,goodWindowsI);
    plagsMain(:,n:n+summi-1) = plags;
    bazMain(n:n+summi-1) = baz;
    velMain(n:n+summi-1) = vel;
    if Gcolumns > 2
        belMain(n:n+summi-1) = bel;
    end

    fprintf("-----------------------------------------------\n");
    fprintf("processed <strong>%d</strong> events on day: %s\n",summi,datestr(day_));
    fprintf("-----------------------------------------------\n");
    n = n+summi-1;
end

shiftsMain(:,n+1:end) = [];
ampMain(n+1:end) = [];
meanCCsMain(n+1:end) = [];
medCCsMain(n+1:end) = [];
medratioMain(n+1:end) = [];
stdratioMain(n+1:end) = [];
slownessMain(:,n+1:end) = [];
waveformMain(:,n+1:end) = [];
plagsMain(:,n+1:end) = [];
tMain(n+1:end) = [];
bazMain(n+1:end) = [];
velMain(n+1:end) = [];
if Gcolumns > 2
    belMain(n+1:end) = [];
end

%%
close all;
nUsed = sum(isfinite(shiftsMain))';
vI = velMain >= 230 & velMain <= 450 & ampMain >= ampThresh & ...
    medCCsMain >= thresh & nUsed >= 3;

figure(); 
plot(tMain(vI),velMain(vI),'.'); zoom on; grid on; title("Velocity");

figure(); 
%plot(tMain(vI),medCCsMain(vI),'.'); zoom on; grid on; hold on; 
plot(tMain(vI),meanCCsMain(vI),'.'); zoom on; grid on;
legend("Median","Mean");

figure();
bazMain(bazMain < 0) = bazMain(bazMain < 0) + 360;
plot(tMain(vI),bazMain(vI),'.'); zoom on; grid on; title("BAZ");
ax = gca;
ax.YLim = [0 360];
ax.YTick = (0:30:360)';

figure();
subplot(211);
plot(tMain(vI),(0:sum(vI)-1)','.'); zoom on; grid on; title("Cumulative Number");
subplot(212);
nHours = 1; 
rate = t2r(tMain(vI),hours(nHours));
plot(tMain(vI),rate,'.'); zoom on; grid on; title("Hourly Rate");

figure(); 
plot(tMain(vI),ampMain(vI),'.'); zoom on; grid on; title("Amplitudes");
hold on;
plot(tMain(~vI),ampMain(~vI),'.');
ax = gca; ax.YScale = 'log';

%%
% close all;
% lS = nComps-length(badComponents);
% n = size(tMain,1);
% for i = 1:nComps %lS
%     figure(1);
%     hold on; ll = plot(sort(shiftsMain(i,:))/Fs2,(0:n-1)'/n,'linewidth',2);
%     ll.Color(4) = 0.8; zoom on;
%     fprintf('Sensor %d: Averge Shift: %g\n',i,median(shiftsMain(i,:))/Fs2);
% end
% 
% %
% figure('units','normalized','outerposition',[0 0 1 1]);
% SS = scatter(tMain,ampMain*sensitivity,100,medCCsMain,'filled');
% zoom on; grid on; ax = gca; ax.YScale = 'log';
% SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5; SS.MarkerFaceAlpha = 0.5; colorbar; caxis([0.8 1])
% 
% %
% % shiftsMain = shiftsMain - shiftsMain(nComps-length(badComponents),:);
% % plagsMaster = getDD(shiftsMain);
% %
% % vel = NaN(size(medCCsMain));
% % baz = vel;
% % bel = vel;
% %
% % for i = 1:size(plagsMaster,2)
% %     %
% %     pm = -plagsMaster(:,i)/Fs2;
% %     %slow_ = lscov(G,pm);
% %     %slow_ = pinv(G)*pm;
% %     slow_ = pm\G;
% %     %slow_ = mean(G./pm);
% %     slownessMaster(:,i) = slow_';
% %     %
% %
% %     if size(G,2) > 2
% %         [azimuth,elevation,r] = cart2sph(slow_(1),slow_(2),slow_(3));
% %         vel(i,1) = r; %rssq(slow_);
% %         baz(i,1) = 90 - rad2deg(azimuth);
% %         bel(i,1) = rad2deg(elevation);
% %         if bel(i) > 90
% %             bel(i) = bel(i) - 90;
% %         end
% %         fprintf('%g %g %g\n',vel(i),baz(i),bel(i));
% %     else
% %         baz(i,1) = atan2d(slow_(1),slow_(2))+360;
% %         if baz(i) >360
% %             baz(i) = baz(i)-360;
% %         end
% %         vel(i,1) = rssq(slow_);
% %     end
% %
% %     if baz(i) > 180
% %         baz(i) = baz(i) - 360;
% %     end
% % end
% 
% strictI = medCCsMain >= 0.85 & vel >= 230 & vel <= 430 & ampMain*sensitivity >= 2*ampThresh; % & baz >= -65 & baz <= -5; % & bel <= 4 & baz >= -65 & baz <= -35;
% 
% figure('units','normalized','outerposition',[0 0 1/2 1]);
% 
% for i = 1:lS
%     ax(i) = subplot(lS,1,i);
%     plot(ax(i),tMain(strictI),shiftsMain(i,strictI)/Fs2,'.');
%     zoom on; grid on; ylabel('lag [sec.]'); title(['Node: ',num2str(i)]);
% end
% linkaxes(ax,'xy'); zoom on;
% 
% %
% figure('units','normalized','outerposition',[0 0 1 1]);
% plot(tMain(strictI),baz(strictI),'o'); zoom on; grid on; ylim([-180 180]);
% hold on;
% plot(tMain(~strictI),baz(~strictI),'.'); zoom on; grid on; ylim([-180 180]);
% 
% %
% figure('units','normalized','outerposition',[0 0 1 1]);
% ax9 = subplot(211);
% SS = scatter(tMain(strictI),ampMain(strictI)*sensitivity,100,baz(strictI),'filled');
% zoom on; grid on; ax = gca; ax.YScale = 'log';
% SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5; SS.MarkerFaceAlpha = 0.5; colorbar; caxis([-180 180]);
% ylabel('Presion [Pa.]');
% 
% ax9(2) = subplot(212);
% histogram(tMain(strictI),(dateshift(min(tMain(strictI)),'start','day'):hours(1):dateshift(max(tMain(strictI)),'end','day'))'); zoom on; grid on;
% ylabel('Numero cada Hora');
% cbar_dummy = colorbar;
% cbar_dummy.Visible = 'off';
% linkaxes(ax9,'x');
% 
% if size(Gorig,2) > 2
%     figure();
%     plot(baz(strictI),bel(strictI),'o'); zoom on;
%     hold on; plot(baz(~strictI),bel(~strictI),'.'); zoom on; grid on;
% end
% 
% figure();
% plot(tMain(strictI),(vel(strictI))','o'); zoom on;
% hold on;
% plot(tMain(~strictI),(vel(~strictI))','.'); zoom on;
% 
% tGood = tMain(strictI);
% aGood = ampMain(strictI)*sensitivity;
% 
% figure();
% clear sax;
% sax(1) = subplot(211);
% plot(tGood,(0:sum(strictI)-1)'/sum(strictI),'.'); zoom on; grid on;
% sax(2) = subplot(212);
% histogram(tGood,dateshift(min(tGood),'start','day'):days(1):dateshift(max(tGood),'end','day')+1); zoom on; grid on; linkaxes(sax,'x');
% 
% Ncut = 101;
% [~,startIndex,endIndex,~,~] = cutWindows(datenum(tGood),Ncut,Ncut-1,false);
% aCut = cutWindows(aGood,Ncut,Ncut-1,false);
% 
% figure();
% plot(tGood(endIndex),mean(aCut)'.*3600*Ncut./seconds(tGood(endIndex)-tGood(startIndex)),'.'); zoom on; grid on;
% close(2);
