clear; close all;
template_search_directory = fullfile("~","masa","template_search");
cd(template_search_directory);
custom_prefix = "napo";
if ~exist(custom_prefix,"dir")
    mkdir(custom_prefix);
end
cd(fullfile(template_search_directory,custom_prefix));

[evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = readCat1(datetime(2025,01,01));

E = readSCBulletin("igepn2025ceos");
width = 0.15;
height = 0.15;
lI = eqlat >= E.lat-height & eqlat <= E.lat+height & eqlon >= E.lon-width & eqlon <= E.lon+width;

table(evDescription(lI),evStatus(lI),t(lI),eqlat(lI),eqlon(lI),eqdepth(lI),...
    eqmag(lI),ids(lI),stderr(lI),azgap(lI),phaseTot(lI),nMLv(lI),timerr(lI),...
    eqlaterr(lI),eqlonerr(lI),eqdeptherr(lI),eqmagerr(lI),magType(lI),meanMethod(lI),...
    locMethod(lI),earthModel(lI),creationTime(lI),agencyID(lI),evMode(lI),...
    scbullmag(lI),authorID(lI),evType(lI))

%%
totEvents = sum(lI);
lI = find(lI);
lfc = 1;
hfc = 16;
newFs = 100;
noiseWin = 2;
secDur = 18;
diffFlag = true;

kstnmList = "PIAT"; %["PIAT";"PIS1";"VCES"];
chanList = ["HHZ";"HHN";"HHE"];
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
    tFinal(n) = min(tTmp(locb(pI)));
    magFinal(n) = E_.mag;
    E(n) = E_;
    %disp([n nStations_ id_]);
end

idFinal = idFinal(1:n);
nStations = nStations(1:n);
E = E(1:n);
tFinal = tFinal(1:n);
magFinal = magFinal(1:n);
table(tFinal,idFinal,nStations,magFinal)

close all;
magFact = 2;
figure(); SS = scatter(pull(E,'lon'),pull(E,'lat'),magFact.*exp(pull(E,'mag')),...
    pull(E,'depth'),'filled'); zoom on; grid on; axis equal;
SS.MarkerEdgeColor = 'k'; SS.MarkerFaceAlpha = 0.5; SS.MarkerEdgeAlpha = 0.5;
figure(); plot(pull(E,'t'),pull(E,'mag'),'.'); zoom on; grid on;
figure(); plot(pull(E,'t'),(0:n-1)','.'); zoom on; grid on;

%%
lS = length(kstnmList);
%lSNCLs = lS * 3; %there are always 3 channels
Pf = extractWaveforms(tFinal-seconds(noiseWin),seconds(secDur),kstnmList,...
    chanList,"EC","",true,true,2,true,[lfc,hfc,false,false,false,newFs]);
refs = pull(Pf,'ref');
refs1 = max(refs,[],2);
badI = logical(sum(isnat(refs),2));

idBad = idFinal(badI);
Ebad = E(badI);

Pf(badI) = [];
idFinal(badI) = [];
nStations(badI) = [];
E(badI) = [];
tFinal(badI) = [];
magFinal(badI) = [];
Pf = Pf';
if diffFlag
    Pf = differentiateWaveforms(Pf);
end
Pf = nanGapWaveforms(Pf);

for i = 1:n
    Pf_ = Pf(:,i);
    Pf(:,i) = push(Pf_,"evid",repmat(idFinal(i),length(Pf_),1));
    Pf_ = Pf(:,i);
    Pf(:,i) = push(Pf_,"eqmag",repmat(magFinal(i),length(Pf_),1));
end

%%
dayStart = datetime(2025,01,01);
dayEnd = datetime(2025,11,05);
dayInc = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

writeFlag = true;
plotFlag = false;
verboseFlag = true;
madThresh = 10;
fileVersionNumber = 2;

refreshData = true;
if refreshData
    parfor i = 1:lDays
        day_ = dayVec(i);
        S = loadWaveforms(day_,1,kstnmList,chanList);
        PfTmp = Pf(:,1);
        Psncls = strcat(pull(PfTmp,'knetwk'),pull(PfTmp,'kstnm'),pull(PfTmp,'khole'),pull(PfTmp,'kcmpnm'));
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
variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";...
    "templateNumber";"madMain";"nUsedMain";"polMain";"evidMain"};
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain,polMain,evidMain] = ...
    loadRepeaterCatalog(custom_prefix,variableNames);
[tMain,madMain,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain] = ...
    filterCatalog(tMain,madMain,20,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain);

%%
%ccI = (ccMain >= 0.16 & madMain >= 16) & magMain >= -1 & nUsedMain >= 3;
ccI = (ccMain >= 0.5 & madMain >= 30) & magMain >= 0 & nUsedMain >= 3;

close all;
axL = linkedPlot(tMain(ccI),[(0:sum(ccI)-1)' ccMain(ccI) ampMain(ccI) madMain(ccI) ...
    magMain(ccI) templateNumber(ccI) nUsedMain(ccI)],".","compact"); zoom on;
axL(3).YScale = 'log';
axL(4).YScale = 'log';

%
% figure();
% plot(tMain(ccI),nUsedMain(ccI),'.','linewidth',3); zoom on; grid on;
% ylabel('n_used')
% %
% figure();
% histogram(tMain(ccI),...
%     (dateshift(min(tMain(ccI)),'start','hour'):hours(1):dateshift(max(tMain(ccI)),'end','hour'))');
% zoom on; title('numero de detecciones cada hora'); grid on;

%
%fig = helicorder(clipWaveforms(Tf(1),2e3)); fig.Visible = 'on';

figure();
tiledlayout(2,1,'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile();
histogram(tMain(ccI),...
    (dateshift(min(tMain(ccI)),'start','day'):days(1):dateshift(max(tMain(ccI)),'end','day'))');
ax(1).XTick = [];
yyaxis right;
plot(tMain(ccI),(0:sum(ccI)-1)','.');
ax(1).YAxis(1).Label.String = "Numero Diario";
ax(1).YAxis(2).Label.String = "Numero Acumulativo";

ax(2,1) = nexttile();
plot(tMain(ccI),magMain(ccI),'o','linewidth',1); zoom on;
ax(2).YLabel.String = "magnitud";
linkaxes(ax,"x");
grid(ax,"on");
xlim tight
hold(ax(2),"on")

%%
Results = table((1:sum(ccI))',tMain(ccI),ccMain(ccI),ampMain(ccI),dMag(ccI),...
    magMain(ccI),templateNumber(ccI),madMain(ccI),polMain(ccI),evidMain(ccI),nUsedMain(ccI));
disp(Results)
sum(ccI)
sum(magMain(ccI) >= 3.5)

%%
ccIOrig = ccI;
ccI = ccI & tMain <= datetime(2025,01,31,23,00,00);
Results2 = table((1:sum(ccI))',tMain(ccI),ccMain(ccI),ampMain(ccI),dMag(ccI),magMain(ccI),templateNumber(ccI),madMain(ccI));
disp(Results2)
ccI = ccIOrig;

%%
Pf2 = extractWaveforms(tMain(ccI),seconds(secDur),"PIAT",["HHZ";"HHN";"HHE"],"EC","",true,true,2,true,[lfc,hfc,false,false,false,newFs]);

%%
snrBig = NaN(size(Pf2)); locsBig = snrBig;
for i = 1:3
    P_ = Pf2(:,i);
    for j = 1:size(P_,1)
        [locs,snr,staOlta,sosSTA] = stalta(P_(j),0.8*noiseWin,0.8*noiseWin,1.5,true,false);
        if isempty(locs)
            continue;
        end
        snr_ = snr./abs(locs-(noiseWin*newFs));
        [~,maxI] = max(snr_);
        snrBig(j,i) = snr(maxI);
        locsBig(j,i) = locs(maxI);
    end
end

ccI = ccIOrig & isfinite(locsBig(:,1)) & locsBig(:,1) < 300;
figure(); semilogy(tMain,max(snrBig,[],2),'.'); zoom on; grid on;
figure(); semilogy(tMain,locsBig(:,1),'.'); zoom on; grid on;