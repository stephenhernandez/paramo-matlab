clear; close all;
strong_motion_dir = fullfile("~","masa","strong_motion");
cd(strong_motion_dir);

yearList = dir('2024*');
lYears = length(yearList);
eventListMain = [];
for i = 1:lYears
    thisYear = yearList(i).name;
    cd(thisYear);
    monthList = dir();
    monthList(ismember( {monthList.name}, {'.', '..'})) = [];
    lMonths = length(monthList);
    for j = 1:lMonths
        thisMonth = monthList(j).name;
        cd(thisMonth);
        eventList = dir();
        eventList(ismember( {eventList.name}, {'.', '..'})) = [];
        lEvents = length(eventList);
        for k = 1:lEvents
            thisEvent = eventList(k).name;
            cd(thisEvent);
            pngs = dir('*.png');
            %lia = ismember('AmplitudeVsDistance_ACC.png',pull(pngs,'name',""));
            lia = ismember('AmplitudeVsDistance_WA.png',pull(pngs,'name',""));
            if ~lia
                cd ..;
                continue;
            end
            eventid = split(thisEvent,'_');
            eventid = string(eventid(1));
            eventListMain = [eventListMain; eventid];
            %!cat summary_WA.txt
            cd ..
        end
        cd ..;
    end
    cd ..;
end

%%
lID = length(eventListMain);
E = populateSeisCompStructure(lID);
for i = 1:lID
    disp(i);
    event_ = eventListMain(i);
    E(i) = readSCBulletin(event_);
end

%%
stderr = pull(E,'rms');
nPhases = pull(E,'usedPhases');
nMag = pull(E,'nmag');
nSphases = pull(E,'nSphases');
magerr2 = pull(E,'magerr2');
timerr = pull(E,'timerr');
laterr = pull(E,'laterr');
lonerr = pull(E,'lonerr');
deptherr = pull(E,'deptherr');
magerr1 = pull(E,'magerr');
t = pull(E,'t');
mag = pull(E,'mag');
azgap = pull(E,'azgap');
depth = pull(E,'depth');
lon = pull(E,'lon');
lat = pull(E,'lat');
evType = pull(E,'evType');
totErr = sqrt(laterr.^2 + lonerr.^2 + deptherr.^2);

%%
allIDs = pull(E,'id');
minT = datetime(2000,01,01);
maxT = datetime(2024,12,31);
minNmlv = 10;
azGapThresh = 360;
minMag = -4.0; %4.5;
maxMag = 86.5;
minPhases = 20; %10;
timeErrThresh =  2; %1.5;
rmsThresh = timeErrThresh;
spatialErrThresh = 20; %10;
errorDepthRatio = 0.1; %0.15;
magerrThresh1 = 0.5; %0.4;
magerrThresh2 = 0.5; %0.4;
minDeptherr = -0.2;
minDepth = -6;
maxDepth = 250;
minSphases = 4; %2;
minTimeErr = 0;

%%
tI = t >= minT & t <= maxT & nMag >= minNmlv & ...
    nPhases >= minPhases & mag <= maxMag & mag >= minMag & azgap <= azGapThresh & ...
    abs(magerr1) <= magerrThresh1 & abs(magerr2) <= magerrThresh2 & ...
    timerr >= minTimeErr & timerr <= timeErrThresh & stderr <= rmsThresh & ...
    depth > minDepth & depth <= maxDepth & deptherr > minDeptherr & ...
    (totErr < spatialErrThresh | totErr./abs(depth) <= errorDepthRatio) & ...
    nSphases >= minSphases & ~(evType == "outsideofnetworkinterest");
sum(tI)

%%
IDfiltered = allIDs(tI);
lIDs = length(IDfiltered);
tF = t(tI);
magF = mag(tI);
depthF = depth(tI);
latF = lat(tI);
lonF = lon(tI);

VariableNames = {'NET_STATION_LOCID_CMPNM','MaxWADisp','MedianWADisp',...
    'TimeOfMaxDisp','DIST','AZ','BAZ','PeakVecDisp','EstimatedDuration'};
%VariableNames = {'NET_STATION_LOCID_CMPNM','MaxAcc','TimeOfMaxAcc','DIST',...
%    'AZ','BAZ','PeakVecACC','EstimatedDuration_'};

sensorTypeList = ["BH";"HN";"HH";"HN";"EN";"BL"];
minDur = 10; %2;
maxDur = 120; %180; %seconds
maxAmpAcc = 2e3; %12e3;
maxAmpSeis = 5e2;
minAmp = 0.1; %0.03;
MinStations = 10; %4;
maxDist = 8e2; %7e2;
MinComponents = 3;
magFact = 3;

%%
cd(strong_motion_dir);
rMain = [];
azMain = [];
bazMain = [];
peakVecMain = [];
durationMain = [];
idMain = [];
tMain = [];
magMain = [];
snclMain = [];
maxCmpAmpMain = [];
medAmpMain = [];
depthMain = [];
latMain = [];
lonMain = [];

for i = 1:lYears
    thisYear = yearList(i).name;
    cd(thisYear);
    monthList = dir();
    monthList(ismember( {monthList.name}, {'.', '..'})) = [];
    lMonths = length(monthList);
    for j = 1:lMonths
        thisMonth = monthList(j).name;
        cd(thisMonth);
        eventList = dir();
        eventList(ismember( {eventList.name}, {'.', '..'})) = [];
        lEvents = length(eventList);
        for k = 1:lEvents
            thisEvent = eventList(k).name;
            cd(thisEvent);
            pngs = dir('*.png');
            lia = ismember('AmplitudeVsDistance_WA.png',pull(pngs,'name',""));
            %lia = ismember('AmplitudeVsDistance_ACC.png',pull(pngs,'name',""));
            if ~lia
                cd ..;
                continue;
            end
            eventid = split(thisEvent,'_');
            t_ = datetime(string(eventid(end)));
            eventid = string(eventid(1));
            [lia,locb] = ismember(eventid,IDfiltered);
            if ~lia
                cd ..;
                continue;
            end

            %%
            disp([t_ tF(locb)]);
            T = readtable('summary_WA.txt');
            %T = readtable('summary_ACC.txt');
            T.Properties.VariableNames = VariableNames;
            NET_STATION_LOCID_CMPNM = string(T.NET_STATION_LOCID_CMPNM);
            NET_STATION_LOCID = NET_STATION_LOCID_CMPNM;
            sensorType = NET_STATION_LOCID_CMPNM;
            for ii = 1:length(NET_STATION_LOCID_CMPNM)
                NET_STATION_LOCID_ = NET_STATION_LOCID_CMPNM(ii);
                sensorType_ = split(NET_STATION_LOCID_,".");
                sensorType_ = char(sensorType_(end));
                sensorType(ii) = string(sensorType_(1:end-1));
                NET_STATION_LOCID_ = char(NET_STATION_LOCID_);
                NET_STATION_LOCID(ii) = string(NET_STATION_LOCID_(1:end-1));
            end
            cmpnmI = ismember(sensorType,sensorTypeList);
            duration_ = T.EstimatedDuration;
            amp_ = T.PeakVecDisp;
            amp2_ = T.MedianWADisp;
            dists_ = T.DIST;
            cmpnmI = (cmpnmI & duration_ >= minDur & duration_ <= maxDur & ...
                amp2_ >= minAmp & amp_ <= maxAmpAcc & dists_ <= maxDist) & ...
                ~((sensorType == "BH" | sensorType == "HH") & amp_ >= maxAmpSeis);
            NET_STATION_LOCID = NET_STATION_LOCID(cmpnmI);
            T = T(cmpnmI,:);

            %%
            uniqueMonitoringSensors = unique(NET_STATION_LOCID);
            lUniqueMonitoringSensors = length(uniqueMonitoringSensors);
            if lUniqueMonitoringSensors < MinStations
                cd ..;
                continue;
            end

            %%
            n = 0;
            for ii = 1:lUniqueMonitoringSensors
                thisSNCL = uniqueMonitoringSensors(ii);
                lia = ismember(NET_STATION_LOCID,thisSNCL);
                nSNCLs = sum(lia);
                if nSNCLs < MinComponents
                    fprintf('not enough channels for: %s\n',thisSNCL);
                    continue;
                end
                n = n+1;
            end

            if n < MinStations
                fprintf('not enough monitoring points for event: %s\n',eventid);
                cd ..;
                continue;
            end

            %% the bulk of the work
            for ii = 1:lUniqueMonitoringSensors
                thisSNCL = uniqueMonitoringSensors(ii);
                lia = ismember(NET_STATION_LOCID,thisSNCL);
                nSNCLs = sum(lia);
                if nSNCLs < MinComponents
                    fprintf('not enough channels for: %s\n',thisSNCL);
                    continue;
                end
                T_ = T(lia,:);
                rMain = [rMain; median(T_.DIST)];
                azMain = [azMain; median(T_.AZ)];
                bazMain = [bazMain; median(T_.BAZ)];
                peakVecMain = [peakVecMain; median(T_.PeakVecDisp)];
                durationMain = [durationMain; median(T_.EstimatedDuration)];
                idMain = [idMain; eventid];
                depthMain = [depthMain; depthF(locb)];
                latMain = [latMain; latF(locb)];
                lonMain = [lonMain; lonF(locb)];
                tMain = [tMain; tF(locb)];
                magMain = [magMain; magF(locb)];
                snclMain = [snclMain; thisSNCL];
                maxCmpAmpMain = [maxCmpAmpMain; max(T_.MaxWADisp)];
                medAmpMain = [medAmpMain; max(T_.MedianWADisp)];
            end
            cd ..
        end
        cd ..;
    end
    cd ..;
end

%%
MinGoodStationObservations = 10;
goodUniqueSNCLs = unique(snclMain);
uniqSNCLcount = groupcounts(snclMain);
goodUniqueSNCLs = goodUniqueSNCLs(uniqSNCLcount >= MinGoodStationObservations);

%%
%cd ~/products/EcuadorStrongMotion
cd(strong_motion_dir);
rMain = [];
azMain = [];
bazMain = [];
peakVecMain = [];
durationMain = [];
idMain = [];
tMain = [];
magMain = [];
snclMain = [];
maxCmpAmpMain = [];
medAmpMain = [];
depthMain = [];
latMain = [];
lonMain = [];

for i = 1:lYears
    thisYear = yearList(i).name;
    cd(thisYear);
    monthList = dir();
    monthList(ismember( {monthList.name}, {'.', '..'})) = [];
    lMonths = length(monthList);
    for j = 1:lMonths
        thisMonth = monthList(j).name;
        cd(thisMonth);
        eventList = dir();
        eventList(ismember( {eventList.name}, {'.', '..'})) = [];
        lEvents = length(eventList);
        for k = 1:lEvents
            thisEvent = eventList(k).name;
            cd(thisEvent);
            pngs = dir('*.png');
            %lia = ismember('AmplitudeVsDistance_ACC.png',pull(pngs,'name',""));
            lia = ismember('AmplitudeVsDistance_WA.png',pull(pngs,'name',""));
            if ~lia
                cd ..;
                continue;
            end
            eventid = split(thisEvent,'_');
            t_ = datetime(string(eventid(end)));
            eventid = string(eventid(1));
            [lia,locb] = ismember(eventid,IDfiltered);
            if ~lia
                cd ..;
                continue;
            end

            %%
            disp([t_ tF(locb)]);
            %T = readtable('summary_ACC.txt');
            T = readtable('summary_WA.txt');
            T.Properties.VariableNames = VariableNames;
            dists_ = T.DIST;
            [~,sortDistI] = sort(dists_);
            T = T(sortDistI,:);
            NET_STATION_LOCID_CMPNM = string(T.NET_STATION_LOCID_CMPNM);
            NET_STATION_LOCID = NET_STATION_LOCID_CMPNM;
            sensorType = NET_STATION_LOCID_CMPNM;
            for ii = 1:length(NET_STATION_LOCID_CMPNM)
                NET_STATION_LOCID_ = NET_STATION_LOCID_CMPNM(ii);
                sensorType_ = split(NET_STATION_LOCID_,".");
                sensorType_ = char(sensorType_(end));
                sensorType(ii) = string(sensorType_(1:end-1));
                NET_STATION_LOCID_ = char(NET_STATION_LOCID_);
                NET_STATION_LOCID(ii) = string(NET_STATION_LOCID_(1:end-1));
            end
            cmpnmI = ismember(sensorType,sensorTypeList) & ...
                ismember(NET_STATION_LOCID,goodUniqueSNCLs);
            duration_ = T.EstimatedDuration;
            amp_ = T.PeakVecDisp;
            amp2_ = T.MedianWADisp;
            dists_ = T.DIST;
            cmpnmI = (cmpnmI & duration_ >= minDur & duration_ <= maxDur & ...
                amp2_ >= minAmp & amp_ <= maxAmpAcc & dists_ <= maxDist) & ...
                ~((sensorType == "BH" | sensorType == "HH") & amp_ >= maxAmpSeis);
            NET_STATION_LOCID = NET_STATION_LOCID(cmpnmI);
            T = T(cmpnmI,:);

            %%
            uniqueMonitoringSensors = unique(NET_STATION_LOCID);
            lUniqueMonitoringSensors = length(uniqueMonitoringSensors);
            if lUniqueMonitoringSensors < MinStations
                cd ..;
                continue;
            end

            %%
            n = 0;
            for ii = 1:lUniqueMonitoringSensors
                thisSNCL = uniqueMonitoringSensors(ii);
                lia = ismember(NET_STATION_LOCID,thisSNCL);
                nSNCLs = sum(lia);
                if nSNCLs < MinComponents
                    fprintf('not enough channels for: %s\n',thisSNCL);
                    continue;
                end
                n = n+1;
            end

            if n < MinStations
                fprintf('not enough monitoring points for event: %s\n',eventid);
                cd ..;
                continue;
            end

            %% the bulk of the work
            for ii = 1:lUniqueMonitoringSensors
                thisSNCL = uniqueMonitoringSensors(ii);
                lia = ismember(NET_STATION_LOCID,thisSNCL);
                nSNCLs = sum(lia);
                if nSNCLs < MinComponents
                    fprintf('not enough channels for: %s\n',thisSNCL);
                    continue;
                end
                T_ = T(lia,:);
                rMain = [rMain; median(T_.DIST)];
                azMain = [azMain; median(T_.AZ)];
                bazMain = [bazMain; median(T_.BAZ)];
                peakVecMain = [peakVecMain; median(T_.PeakVecDisp)];
                durationMain = [durationMain; median(T_.EstimatedDuration)];
                idMain = [idMain; eventid];
                depthMain = [depthMain; depthF(locb)];
                latMain = [latMain; latF(locb)];
                lonMain = [lonMain; lonF(locb)];
                tMain = [tMain; tF(locb)];
                magMain = [magMain; magF(locb)];
                snclMain = [snclMain; thisSNCL];
                maxCmpAmpMain = [maxCmpAmpMain; max(T_.MaxWADisp)];
                medAmpMain = [medAmpMain; max(T_.MedianWADisp)];
            end
            cd ..
        end
        cd ..;
    end
    cd ..;
end

%%
close all;
figure(); SS = scatter(rMain,peakVecMain,magFact*exp(magMain),depthMain,'o','filled'); zoom on;
ax = gca;
ax.YScale = 'log'; ax.XScale = 'log';
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.2; SS.MarkerFaceAlpha = 0.5;
grid on; cbar = colorbar;
set(ax,"ColorScale",'log');
clim([1 200]);
colormap turbo;

%
uniqIDs = unique(idMain);
lUniqIDs = length(uniqIDs);
gCount = 1;
gFlag = true;
if gFlag
    gMax = 1e7;
    Gmain = zeros(gMax,length(goodUniqueSNCLs)+2);
    dMain = NaN(gMax,1);
end
cumG = 0;

for i = 1:lUniqIDs
    id_ = uniqIDs(i);
    uI = idMain == id_;
    nAmps = sum(uI);
    sncls_ = snclMain(uI);
    amp_ = medAmpMain(uI);
    %amp_ = peakVecMain(uI);
    r_ = rMain(uI);
    [r_,sortDistI] = sort(r_);
    amp_ = amp_(sortDistI);
    sncls_ = sncls_(sortDistI);

    [~,goodSnclsI] = ismember(sncls_,goodUniqueSNCLs);
    diffCount = nAmps*(nAmps - 1)*0.5;
    if gFlag
        Gtmp = full(Gvdcc(nAmps));
        Gtmp = Gtmp(1:end-1,:);
        difflog = -getDD(log10(r_));
        difflin = -getDD(r_);
        difflogamp = getDD(log10(amp_));
        Gmain(gCount:gCount+diffCount-1,1) = difflog;
        Gmain(gCount:gCount+diffCount-1,2) = difflin;
        dMain(gCount:gCount+diffCount-1) = difflogamp;
        for j = 1:nAmps
            Gmain(gCount:gCount+diffCount-1,goodSnclsI(j)+2) = Gtmp(:,j);
        end
    end
    gCount = gCount + diffCount;
    fprintf("%d/%d %d %d\n",i,lUniqIDs,nAmps,gCount-1);
end

if gFlag
    dMain = dMain(1:gCount-1);
    Gmain = Gmain(1:gCount-1,:);
end

Gmain = [Gmain; [0 0 ones(1,length(goodUniqueSNCLs))]];
dMain = [dMain; 0];
solution = Gmain\dMain;

rdum = logspace(log10(5),log10(700),401);
att1 = solution(1)*log10(rdum) + solution(2)*rdum;  %hernandez
att2 = 1.11*log10(rdum) + 0.00189*rdum + 0.591;     %uhrhammer

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

%%