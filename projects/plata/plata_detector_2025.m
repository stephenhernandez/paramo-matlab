clear; close all;
template_search_directory = fullfile("~","masa","template_search");
custom_prefix = "ISPT2025";
cd(fullfile(template_search_directory,custom_prefix));

[evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = ...
    readCat1(datetime(2025,01,01));

%E = readSCBulletin("igepn2025pzoo");
origLat = -1.45; %E.lon;
origLon = -81.05; %E.lat;
width = 0.3;
height = 0.3;
lI = eqlat >= origLat-height & eqlat <= origLat+height & ...
    eqlon >= origLon-width & eqlon <= origLon+width;

table((1:sum(lI))',evDescription(lI),evStatus(lI),t(lI),eqlat(lI),eqlon(lI),eqdepth(lI),...
    eqmag(lI),ids(lI),stderr(lI),azgap(lI),phaseTot(lI),nMLv(lI),timerr(lI),...
    eqlaterr(lI),eqlonerr(lI),eqdeptherr(lI),eqmagerr(lI),magType(lI),meanMethod(lI),...
    locMethod(lI),earthModel(lI),creationTime(lI),agencyID(lI),evMode(lI),...
    scbullmag(lI),authorID(lI),evType(lI))

totEvents = sum(lI);
lI = find(lI);
lfc = 2;
hfc = 10;
newFs = 50; %hfc*8;
noiseWin = 0;
secDur = 35;
diffFlag = true;

kstnmList = ["JIPI";"SALI";"ISPT";"PPLP"];
chanList = ["HHZ";"HHN";"HHE";"BLZ";"BLN";"BLE"];
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
    if ~nStations_
        continue;
    end

    %
    n = n+1;
    tTmp = pull(Pphases,"t");
    idFinal(n) = id_;
    nStations(n) = nStations_;
    tFinal(n) = E_.t; %min(tTmp(locb(pI)));
    magFinal(n) = E_.mag;
    E(n) = E_;
    sprintf("%d %d %s\n",n,nStations_,id_);
end
idFinal = idFinal(1:n);
nStations = nStations(1:n);
E = E(1:n);
tFinal = tFinal(1:n);
magFinal = magFinal(1:n);

%
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
    chanList,"EC","",true,true,2,true,[lfc,hfc,false,false,false,newFs]);

delSNCL = "EC.PPLP.HH";
for i = 1:n
    Pf_ = Pf(i,:);
    cmp = char(pull(Pf_(:),'kcmpnm'));
    cmp = string(cmp(:,1:2));
    allSNCLs_ = compose("%s.%s.%s",pull(Pf_(:),"knetwk"),pull(Pf_(:),"kstnm"),cmp);
    ref_ = pull(Pf_(:),'ref');
    goodI = ~isnat(ref_) & (allSNCLs_ ~= delSNCL);
    Pf_(1:sum(goodI)) = Pf_(goodI);
    Pf_(sum(goodI)+1:end) = populateWaveforms(1);
    Pf(i,:) = Pf_;
end

refs = pull(Pf,'ref');
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
if diffFlag
    Pf = differentiateWaveforms(Pf);
end
Pf = nanGapWaveforms(Pf,0);

for i = 1:n
    Pf_ = Pf(:,i);
    Pf_ = syncWaveforms(Pf_,false,true,true);
    Pf(:,i) = push(Pf_,"evid",repmat(idFinal(i),length(Pf_),1));
    Pf_ = Pf(:,i);
    Pf(:,i) = push(Pf_,"eqmag",repmat(magFinal(i),length(Pf_),1));
end

%%
dayStart = datetime(2025,01,01);
dayEnd = datetime(2025,11,20);

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
Psncls = strcat(pull(PfTmp,'knetwk'),pull(PfTmp,'kstnm'),pull(PfTmp,'khole'),pull(PfTmp,'kcmpnm'));

if refreshData
    parfor i = 1:lDays
        day_ = dayVec(i);
        S = loadWaveforms(day_,1,kstnmList,chanList);

        Ssncls = strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm'));

        [lia,locb] = ismember(Ssncls,Psncls);
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
        Sf = resampleWaveforms(Sf,newFs);
        Sf = nanGapWaveforms(Sf,0);
        Sf = padWaveforms(Sf);

        [~,~,~,~,~,~,~,~,~,ccnorm,tLong] = ...
            templateSearch(Sf(lia),Pf(locb(lia),:),writeFlag,plotFlag,...
            verboseFlag,madThresh,custom_prefix,fileVersionNumber);
    end
end

%%
[evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = ...
    readCat1(datetime(2025,01,01));

%E = readSCBulletin("igepn2025pzoo");
origLat = -1.45; %E.lon;
origLon = -81.05; %E.lat;
width = 0.3;
height = 0.3;
lI = eqlat >= origLat-height & eqlat <= origLat+height & ...
    eqlon >= origLon-width & eqlon <= origLon+width;

table((1:sum(lI))',evDescription(lI),evStatus(lI),t(lI),eqlat(lI),eqlon(lI),eqdepth(lI),...
    eqmag(lI),ids(lI),stderr(lI),azgap(lI),phaseTot(lI),nMLv(lI),timerr(lI),...
    eqlaterr(lI),eqlonerr(lI),eqdeptherr(lI),eqmagerr(lI),magType(lI),meanMethod(lI),...
    locMethod(lI),earthModel(lI),creationTime(lI),agencyID(lI),evMode(lI),...
    scbullmag(lI),authorID(lI),evType(lI))

variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";...
    "templateNumber";"madMain";"nUsedMain";"polMain";"evidMain"};
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain,polMain,evidMain] = ...
    loadRepeaterCatalog(custom_prefix,variableNames);
[tMain,madMain,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain] = ...
    filterCatalog(tMain,madMain,15,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain);

%
ccI = (ccMain >= 0.3 & madMain >= 20) & tMain >= datetime(2025,01,01) & ...
    tMain < datetime(2026,01,01) & magMain >= 0 & nUsedMain >= 6;
ccIOrig = ccI;

close all;
axL = linkedPlot(tMain(ccI),[(0:sum(ccI)-1)' ccMain(ccI) ampMain(ccI) madMain(ccI) ...
    magMain(ccI) templateNumber(ccI) nUsedMain(ccI)]); zoom on;
axL(3).YScale = 'log';
axL(4).YScale = 'log';

%%
figure();
plot(tMain(ccI),nUsedMain(ccI),'.','linewidth',3); zoom on; grid on;
ylabel('n_used')

%
figure();
histogram(tMain(ccI),...
    (dateshift(min(tMain(ccI)),'start','hour'):days(1):dateshift(max(tMain(ccI)),'end','hour'))');
zoom on; title('numero de detecciones cada dia'); grid on;

%fig = helicorder(clipWaveforms(Tf(1),2e3)); fig.Visible = 'on';

figure();
tiledlayout(2,1,'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile();
histogram(tMain(ccI),...
    (dateshift(min(tMain(ccI)),'start','hour'):days(1):dateshift(max(tMain(ccI)),'end','hour'))');
ax(1).XTick = [];
yyaxis right;
plot(tMain(ccI),(0:sum(ccI)-1)','.');
ax(1).YAxis(1).Label.String = "Numero cada Dia";
ax(1).YAxis(2).Label.String = "Numero Acumulativo";

ax(2,1) = nexttile();
plot(tMain(ccI),magMain(ccI),'o','linewidth',1); zoom on;
ax(2).YLabel.String = "magnitud";
linkaxes(ax,"x");
grid(ax,"on");
xlim tight
hold(ax(2),"on");
plot(t(lI),eqmag(lI),'.'); zoom on;

%
Results = table((1:sum(ccI))',tMain(ccI),ccMain(ccI),ampMain(ccI),dMag(ccI),...
    magMain(ccI),templateNumber(ccI),madMain(ccI),polMain(ccI),evidMain(ccI),nUsedMain(ccI));
disp(Results)
sum(ccI)
sum(magMain(ccI) >= 3)

%
repeat_id_info = [unique(templateNumber(ccI)),groupcounts(templateNumber(ccI))];
[~,sI] = sort(groupcounts(templateNumber(ccI)),"descend");
repeat_id_info = repeat_id_info(sI,:);
table(repeat_id_info)

%%
ccI2 = ccIOrig & templateNumber == 274;
Results2 = table((1:sum(ccI2))',tMain(ccI2),ccMain(ccI2),ampMain(ccI2),dMag(ccI2),...
    magMain(ccI2),templateNumber(ccI2),madMain(ccI2),polMain(ccI2),evidMain(ccI2),nUsedMain(ccI2));
disp(Results2)

% Pf2 = extractWaveforms(tMain(ccI2),seconds(secDur),kstnmList,...
%     chanList,"EC","",true,true,2,true,[lfc,hfc,false,false,false,newFs]);