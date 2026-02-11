clear; close all;
template_search_directory = fullfile("~","masa","template_search");
cd(template_search_directory);
custom_prefix = "FER2025"; %"SN2025";
if ~exist(custom_prefix,"dir")
    mkdir(custom_prefix);
end
cd(fullfile(template_search_directory,custom_prefix));

tStart = datetime(2023,01,01);
[evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = ...
    readCat1(tStart);

%E = readSCBulletin("igepn2025pzoo");
% origLat = -0.811195; %E.lon;
% origLon = -91.13496; %E.lat;
origLat = -0.372; %Feranndina
origLon = -91.54; %Fernandina
width = 0.35;
height = 0.25;
minMag = 0;
lI = eqlat >= origLat-height & eqlat <= origLat+height & ...
    eqlon >= origLon-width & eqlon <= origLon+width & eqmag >= minMag;

table(evDescription(lI),evStatus(lI),t(lI),eqlat(lI),eqlon(lI),eqdepth(lI),...
    eqmag(lI),ids(lI),stderr(lI),azgap(lI),phaseTot(lI),nMLv(lI),timerr(lI),...
    eqlaterr(lI),eqlonerr(lI),eqdeptherr(lI),eqmagerr(lI),magType(lI),meanMethod(lI),...
    locMethod(lI),earthModel(lI),creationTime(lI),agencyID(lI),evMode(lI),...
    scbullmag(lI),authorID(lI),evType(lI))

totEvents = sum(lI);
lI = find(lI);
lfc = 3;
hfc = 12;
newFs = 100;
noiseWin = 1;
secDur = 10;
diffFlag = true;

%kstnmList = ["VCH1";"SN12"];
%chanList = ["HHZ";"HHN";"HHE"];
kstnmList = ["FER1";"FER3"];
chanList = ["HHZ";"HHN";"HHE";"BHZ";"BHN";"BHE"];

n = 0;
E = populateSeisCompStructure(totEvents);
nStations = NaN(totEvents,1);
idFinal = repmat("",totEvents,1);
tFinal = NaT(totEvents,1);
magFinal = NaN(totEvents,1);

for i = 1:totEvents
    id_ = ids(lI(i));
    E_ = readSCBulletin(id_);
    Pphases = E_.Pphases;
    stnms = pull(Pphases,"stnm","");
    [pI,locb] = ismember(kstnmList,stnms);
    nStations_ = sum(pI);
    if nStations_ < length(kstnmList)
        continue;
    end

    %
    n = n+1;
    tTmp = pull(Pphases,"t");
    idFinal(n) = id_;
    nStations(n) = nStations_;
    tFinal(n) = min(tTmp(locb(pI)));
    magFinal(n) = E_.mag;
    E(n) = E_;
    sprintf("%d %d %s\n",n,nStations_,id_);
end
idFinal = idFinal(1:n);
nStations = nStations(1:n);
E = E(1:n);
tFinal = tFinal(1:n);
magFinal = magFinal(1:n);

close all;
magFact = 5;
figure();
SS = scatter(pull(E,'lon'),pull(E,'lat'),magFact.*exp(pull(E,'mag')),...
    pull(E,'depth'),'filled');
zoom on; grid on; axis equal;
hold on;
plot(origLon,origLat,'kp');

SS.MarkerEdgeColor = 'k';
SS.MarkerFaceAlpha = 0.5;
SS.MarkerEdgeAlpha = 0.5;

figure(); plot(pull(E,'t'),pull(E,'mag'),'.'); zoom on; grid on;
figure(); plot(pull(E,'t'),(0:n-1)','.'); zoom on; grid on;

%%
lS = length(kstnmList);
Pf = extractWaveforms(tFinal-seconds(noiseWin),seconds(secDur),kstnmList,...
    chanList,"EC","",true,true,2,true,[lfc,hfc,false,false,diffFlag,newFs]);

old_cmps = ["BHZ";"BHN";"BHE"];
rename_cmps = ["HHZ";"HHN";"HHE"];
for i = 1:length(old_cmps)
    old_ = old_cmps(i);
    new_ = rename_cmps(i);
    bh_kcmpnmI = strcmp(pull(Pf,"kcmpnm"),old_);
    if ~sum(bh_kcmpnmI)
        fprintf("no %s\n",old_);
        continue;
    end
    Pf_ = push(Pf(bh_kcmpnmI),"kcmpnm",repmat(new_,sum(bh_kcmpnmI,"all"),1),true);
    Pf(bh_kcmpnmI) = Pf_;
end
refs = pull(Pf,"ref");
refs1 = max(refs,[],2);
badI = logical(sum(isnat(refs1),2));
Pf(badI,:) = [];
for i = 1:size(Pf,1)
    Pf_ = Pf(i,:);
    refs_ = pull(Pf_,"ref");
    goodI = ~isnat(refs_);
    Pf_ = [Pf_(goodI) Pf_(~goodI)];
    Pf(i,:) = Pf_;
end

refs = pull(Pf,"ref");
refs1 = max(refs,[],2);
badI = logical(sum(isnat(refs1),2));

clear k;
for i = 1:size(refs,1)
    k(i) = find(~isnat(refs(i,:)),1,"last");
end
k = k';
nGood = max(k);
Pf = Pf(:,1:nGood);
refs = pull(Pf,'ref');

badI = badI | sum(isnat(refs),2) > 0;
Pf(badI,:) = [];
idFinal(badI) = [];
nStations(badI) = [];
E(badI) = [];
tFinal(badI) = [];
magFinal(badI) = [];
n = n - sum(badI);

%
Pf = Pf';
Pf = nanGapWaveforms(Pf,0);

for i = 1:n
    Pf_ = Pf(:,i);
    Pf_ = syncWaveforms(Pf_,false,true,true);
    Pf(:,i) = push(Pf_,"evid",repmat(idFinal(i),length(Pf_),1));
    Pf_ = Pf(:,i);
    Pf(:,i) = push(Pf_,"eqmag",repmat(magFinal(i),length(Pf_),1));
end

%%
tic;
CANUSEGPU = canUseGPU();
%dayStart = datetime(2018,09,01);
%dayEnd = datetime(2025,09,24);
dayStart = datetime(2025,12,23);
dayEnd = datetime(2025,12,23);

dayInc = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);
writeFlag = true;
plotFlag = false;
verboseFlag = true;
madThresh = 10;
fileVersionNumber = 2;

refreshData = true;
PfTmp = Pf(:,1);
Psncls = strcat(pull(PfTmp,"knetwk"),pull(PfTmp,"kstnm"),...
    pull(PfTmp,"khole"),pull(PfTmp,"kcmpnm"));

if refreshData
    for i = 1:lDays
        day_ = dayVec(i);
        S = loadWaveforms(day_,1,kstnmList,chanList,"EC","",false,false); %"~/masa/backups/");
        for j = 1:length(old_cmps)
            old_ = old_cmps(j);
            new_ = rename_cmps(j);
            bh_kcmpnmI = strcmp(pull(S,"kcmpnm"),old_);
            if ~sum(bh_kcmpnmI)
                fprintf("no %s\n",old_);
                continue;
            end
            S_ = push(S(bh_kcmpnmI),"kcmpnm",repmat(new_,sum(bh_kcmpnmI,"all"),1),true);
            S(bh_kcmpnmI) = S_;
        end

        Ssncls = strcat(pull(S,"knetwk"),pull(S,"kstnm"),pull(S,"khole"),pull(S,"kcmpnm"));
        lia = ismember(Ssncls,Psncls);
        dataExist = sum(lia);
        if ~dataExist
            continue;
        end

        fprintf("data exist (%d sncls) for day: %s\n",dataExist,day_);
        if diffFlag
            S = differentiateWaveforms(S);
        end
        
        Sf = syncWaveforms(S,false,true,true);
        Sf = filterWaveforms(detrendWaveforms(Sf),lfc,hfc);
        Sf = resampleWaveforms(Sf,newFs); %must come before merge
        Sf = mergeWaveforms(Sf); %must come after resample
        Sf = nanGapWaveforms(Sf,0);
        Sf = padWaveforms(Sf);

        Ssncls = strcat(pull(Sf,"knetwk"),pull(Sf,"kstnm"),...
            pull(Sf,"khole"),pull(Sf,"kcmpnm"));
        [lia,locb] = ismember(Ssncls,Psncls);

        templateSearch(Sf(lia),Pf(locb(lia),:),writeFlag,plotFlag,...
            verboseFlag,madThresh,custom_prefix,fileVersionNumber,CANUSEGPU);
    end
end
toc;

%%
MarkerEdgeAlpha = 0.5;
MarkerFaceAlpha = 0.5;
variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";...
    "templateNumber";"madMain";"nUsedMain";"polMain";"evidMain"};
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain,polMain,evidMain] = ...
    loadRepeaterCatalog(custom_prefix,variableNames);
[tMain,madMain,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain] = ...
    filterCatalog(tMain,madMain,6,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain);

%ccI = (ccMain >= 0.45 | madMain >= 27) & tMain >= datetime(2017,01,01) & ...
%    tMain < datetime(2026,01,01) & magMain >= -1 & nUsedMain >= 3;
% ccI = (ccMain >= 0.3 & madMain >= 15) & tMain >= datetime(2023,11,01) & ...
%     tMain < datetime(2026,01,01) & magMain >= -1 & ampMain >= 200 & nUsedMain >= 6;
ccI = (ccMain >= 0.25 | madMain >= 14) & tMain >= datetime(2023,11,01) & ...
    tMain < datetime(2027,01,01) & magMain >= -1 & ampMain >= 300 & nUsedMain >= 6;
ccIOrig = ccI;

close all;
axL = linkedPlot(tMain(ccI),[(0:sum(ccI)-1)' ccMain(ccI) ampMain(ccI) madMain(ccI) ...
    magMain(ccI) templateNumber(ccI) nUsedMain(ccI)],".","compact"); zoom on;
axL(3).YScale = "log";
axL(4).YScale = "log";

figure();
tiledlayout(2,1,'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile();
histogram(tMain(ccI),...
    (dateshift(min(tMain(ccI)),'start','hour'):hours(1):dateshift(max(tMain(ccI)),'end','hour'))');
ax(1).XTick = [];
yyaxis right;
plot(tMain(ccI),(0:sum(ccI)-1)','.');
ax(1).YAxis(1).Label.String = "Numero cada Hora";
ax(1).YAxis(2).Label.String = "Numero Acumulativo";

ax(2,1) = nexttile();
magFact = 5;
SS = scatter(tMain(ccI),magMain(ccI),magFact*exp(magMain(ccI)),ccMain(ccI),'o',"filled"); zoom on;
ax(2).YLabel.String = "magnitud";
linkaxes(ax,"x");
grid(ax,"on");
xlim tight
hold(ax(2),"on");
ax(2).Box = "on";
cbar = colorbar;
SS.MarkerEdgeColor = "k";
SS.MarkerEdgeAlpha = MarkerEdgeAlpha;
SS.MarkerFaceAlpha = MarkerFaceAlpha;
plot(tFinal,magFinal,'.'); zoom on;

%
Results = table((1:sum(ccI))',tMain(ccI),ccMain(ccI),ampMain(ccI),dMag(ccI),...
    magMain(ccI),templateNumber(ccI),madMain(ccI),polMain(ccI),evidMain(ccI),nUsedMain(ccI));
disp(Results)
sum(ccI)
sum(magMain(ccI) >= 3)

repeat_id_info = [unique(templateNumber(ccI)),groupcounts(templateNumber(ccI))];
[~,sI] = sort(groupcounts(templateNumber(ccI)),"descend");
repeat_id_info = repeat_id_info(sI,:);
table(repeat_id_info)

%%
ccI2 = ccIOrig & templateNumber == 47;
Results2 = table((1:sum(ccI2))',tMain(ccI2),ccMain(ccI2),ampMain(ccI2),dMag(ccI2),...
    magMain(ccI2),templateNumber(ccI2),madMain(ccI2),polMain(ccI2),evidMain(ccI2),nUsedMain(ccI2));
disp(Results2)
disp(sum(ccI))

%%
PfOrig = Pf;

%%
tic;
Pf = extractWaveforms(tMain(ccI),seconds(secDur),kstnmList,chanList,"EC","",...
    true,true,2,true,[lfc,hfc,false,false,diffFlag,newFs]); toc;

%%
old_cmps = ["BHZ";"BHN";"BHE"];
rename_cmps = ["HHZ";"HHN";"HHE"];
for i = 1:length(old_cmps)
    old_ = old_cmps(i);
    new_ = rename_cmps(i);
    bh_kcmpnmI = strcmp(pull(Pf,"kcmpnm"),old_);
    if ~sum(bh_kcmpnmI)
        fprintf("no %s\n",old_);
        continue;
    end
    Pf_ = push(Pf(bh_kcmpnmI),"kcmpnm",repmat(new_,sum(bh_kcmpnmI,"all"),1),true);
    Pf(bh_kcmpnmI) = Pf_;
end
toc;

refs = pull(Pf,"ref");
refs1 = max(refs,[],2);
badI = logical(sum(isnat(refs1),2));
Pf(badI,:) = [];
for i = 1:size(Pf,1)
    Pf_ = Pf(i,:);
    refs_ = pull(Pf_,"ref");
    goodI = ~isnat(refs_);
    Pf_ = [Pf_(goodI) Pf_(~goodI)];
    Pf(i,:) = Pf_;
end

%%
refs = pull(Pf,"ref");
refs1 = max(refs,[],2);
badI = logical(sum(isnat(refs1),2));

clear k;
for i = 1:size(refs,1)
    k(i) = find(~isnat(refs(i,:)),1,"last");
end
k = k';
nGood = max(k);
Pf = Pf(:,1:nGood);
refs = pull(Pf,'ref');

badI = badI | sum(isnat(refs),2) > 0;
Pf(badI,:) = [];

%%
Pf2 = Pf; clear Pf;
Pf2Orig = Pf2; Pf2Orig = Pf2Orig(:,1:6);

%%
tNew = tMain(ccI);
tOrig = pull(Pf2Orig(:,1),"ref");
[lia, locb] = ismember(tOrig,tNew);
tNew(locb(lia)) = [];
lNew = length(tNew);
if lNew > 0
    tic;
    Pf2 = extractWaveforms(tNew,seconds(secDur),kstnmList,chanList,"EC","",...
        true,true,2,true,[3,12,false,false,true,newFs]);
    old_cmps = ["BHZ";"BHN";"BHE"];
    rename_cmps = ["HHZ";"HHN";"HHE"];
    for i = 1:length(old_cmps)
        old_ = old_cmps(i);
        new_ = rename_cmps(i);
        bh_kcmpnmI = strcmp(pull(Pf2,"kcmpnm"),old_);
        if ~sum(bh_kcmpnmI)
            fprintf("no %s\n",old_);
            continue;
        end
        Pf2_ = push(Pf2(bh_kcmpnmI),"kcmpnm",repmat(new_,sum(bh_kcmpnmI,"all"),1),true);
        Pf2(bh_kcmpnmI) = Pf2_;
    end
    toc;
    Pf2 = Pf2(:,1:6);
    Pf2 = [Pf2Orig; Pf2];
end

%%
Pf2 = Pf2(:,1:6);
nEvents = size(Pf2,1);
P_ = Pf2(1,:);
d_ = pull(normalizeWaveforms(P_(:)));
d_ = normalizeWaveforms(d_(:));

d = zeros(length(d_),nEvents);
for i = 1:nEvents
    P_ = detrendWaveforms(Pf2(i,:));
    d_ = pull(normalizeWaveforms(P_(:)));
    badI = ~isfinite(d_);
    d_(badI) = 0;
    d_ = normalizeWaveforms(detrend(d_(:)));
    ld_ = length(d_);
    d(1:ld_,i) = d_;
end

%%
MAXLAG = round(0.5*secDur*newFs);
d = normalizeWaveforms(detrend(d));
CANUSEGPU = canUseGPU();
N = size(d,2);
d1 = normalizeWaveforms(detrend(d(1:3003,1:N)));
d2 = normalizeWaveforms(detrend(d(3004:end,1:N)));
if CANUSEGPU
    dGPUsingle1 = gpuArray(single(d1));
    dGPUsingle2 = gpuArray(single(d2));

    tic;
    [maxccp_,plags_,pols_] = doccFreqCircShiftPolarities(dGPUsingle1(:,1:N),...
        true,[],MAXLAG);
    toc;
    pause(5)
    tic;
    [maxccp_2,plags_2,pols_2] = doccFreqCircShiftPolarities(dGPUsingle2(:,1:N),...
        true,[],MAXLAG);
    toc;
end

[maxccp_,maxI] = max([maxccp_ maxccp_2],[],2,"linear");
plags_ = [plags_ plags_2];
pols_ = [pols_ pols_2];
plags_ = plags_(maxI);  % this is probably not usable, do not use
pols_ = pols_(maxI);    % this is probably not usable, do not use

maxccp = squareform(maxccp_);
thresh = 0.3;
figure('units','normalized','outerposition',[0 0 1 1]);
h = imagesc(maxccp);
title('Similarity Matrix');
ylabel('Event Number');
xlabel('Event Number');
set(h, 'alphadata', maxccp >= thresh);
axis square;
colorbar;
clim([thresh 1]);
zoom on;

LOGFLAG = false;
MAGFLAG = true;
method = "complete";
[family,l_uniq_indices,nSingletons,multipletNumber] = generateFamilies(maxccp_,0.7,method);
[multipletI,lf] = plotClusterTimeSeries(family,tMain(ccI),magMain(ccI),MAGFLAG,LOGFLAG);

%%
NCHAN = 6;
PWSFLAG = true;
% final parameters commented out to save
% [newEvents,newFamilies,lPerMultiplet,...
%     singletonI,lSingletons,maxccp,plags] = pruneAndMergeEvents(d,...
%    maxccp_,[0.9; 0.9; 0.8; 0.7; 0.6; 0.5],[2; 2; 2; 2; 2; 3],"complete",true);

threshVec = [0.9; 0.9; 0.8; 0.8; 0.7; 0.7; 0.7; 0.7; 0.6; 0.6; 0.6];
[newEvents,newFamilies,lPerMultiplet,singletonI,lSingletons,maxccp] = ...
    pruneAndMergeEvents(d,maxccp_,threshVec,...
    [2*ones(8,1); 3; 3; 3],"complete",PWSFLAG,CANUSEGPU,MAXLAG,NCHAN);

[multipletI,lf] = plotClusterTimeSeries(newFamilies,...
    tMain(ccI),magMain(ccI),MAGFLAG,LOGFLAG);

maxccp = squareform(maxccp);
thresh = 0.7;
figure('units','normalized','outerposition',[0 0 1 1]);
h = imagesc(maxccp);
title('Similarity Matrix');
ylabel('Event Number');
xlabel('Event Number');
set(h, 'alphadata', maxccp >= thresh);
axis square;
colorbar;
clim([thresh 1]);
zoom on;