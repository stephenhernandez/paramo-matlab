function [exitStatus,ccnorm_,t_] = processSNCLGroup(snclGroup,tStart)
if nargin < 2
    tStart = NaT;
end
exitStatus = true;
ccnorm_ = [];
t_ = [];
snclList = snclGroup.snclList;
principal = snclList(1,1);
nrows = size(snclList,1);

%
iFlag = false;
if nrows > 1
    iFlag = strcmp(snclList(2,1),principal) & strcmp(snclList(2,2),"BDF");
end

%%
saveDir = snclGroup.saveDir;
refEllipse = snclGroup.refEllipse;
pauseTime = snclGroup.pauseTime;
lfc = snclGroup.lfc;
hfc = snclGroup.hfc;
newFs = snclGroup.newFs;
faceAlpha = snclGroup.faceAlpha;
templateFileName = snclGroup.templateFileName;
temporaryCatalogFile = snclGroup.temporaryCatalogFile;
longTermCatalogFile = snclGroup.longTermCatalogFile;
shortSnippet = snclGroup.shortSnippet;
threshold = snclGroup.threshold;
maxTemplates = snclGroup.maxTemplates;
recordLength = snclGroup.recordLength;
maxN = snclGroup.maxN;
mpd = snclGroup.mpd;
linearccnorm= snclGroup.linearccnorm;
plotFlag= snclGroup.plotFlag;
verboseFlag = snclGroup.verboseFlag;
diffFlag = snclGroup.diffFlag;

%%
debuggingMode = snclGroup.debuggingMode;
if debuggingMode
    nHours = 6;
    nDays = 4;
    verboseFlag = true;
else
    nHours = snclGroup.nHours;
    nDays = snclGroup.nDays;
end

if isnat(tStart)
    nowTime = datetime('now')+hours(5); %current date/time in UTC
    loopStart = nowTime - days(nDays);
else
    loopStart = tStart;
end

%%
lk = size(snclList,1);

%% read data
Sorig = populateWaveforms(lk);
nn = 0;
for ii = 1:lk
    SNCL_ = snclList(ii,:);
    kstnm_ = SNCL_(1);
    kcmpnm_ = SNCL_(2);
    knetwk_ = SNCL_(3);
    khole_ = SNCL_(4);

    %
    S_ = extractWaveforms(loopStart,seconds(days(nDays)),kstnm_,kcmpnm_,knetwk_,khole_,true,verboseFlag,nDays+1);
    if isnat(S_.ref) %oops, bad found! skip...
        continue;
    end
    nn = nn + 1;

    Sorig(nn,1) = S_;
    if strcmp(saveDir,'sangay')
        if strcmp(S_.kstnm,"SAGA") && strcmp(S_.kcmpnm,"HHZ")
            Shf = filterWaveforms(S_,10);
            SAGA_orig = S_;
        end
    end
end

%%
clear S_ SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
Sorig = Sorig(1:nn,1);
lS = length(Sorig); %<-- also equal to nn
if verboseFlag
    for i = 1:lS
        fprintf(1,'%s: %s - %s\n',Sorig(i).kstnm,datestr(Sorig(i).ref,31),datestr(Sorig(i).ref+Sorig(i).e,31));
    end
    disp(' ');
end

%%
if ~lS
    disp('no waveform data!!');
    exitStatus = false;
    return;
end

%%
try
    [ccnorm_,t_] = updateCustomCatalog(Sorig,templateFileName,temporaryCatalogFile,longTermCatalogFile,...
        threshold,maxTemplates,recordLength,maxN,diffFlag,mpd,...
        linearccnorm,plotFlag,verboseFlag);
    pause(pauseTime);
    if ~isnat(tStart)
        fprintf('short circuit\n')
        try
            updateAllSangayFeatures(true);
        catch
            warning('couldnt update sangay feature set for some reason');
            return;
        end
        return;
    end
catch ME
    fprintf(2,'couldnt update custom catalog\n');
    %warning(ME.message);
    %exitStatus = false;
    rethrow(ME)
    return;
end

%%
[tPotential,~,NCC] = filterUniqueEvents(longTermCatalogFile,mpd);

%% run STA/LTA
tw = 100;
fillValue = NaN;
if strcmp(Sorig(1).kstnm,principal)
    if strcmp(saveDir,'sangay')
        %
        tic;
        SAGA_templateFileName = '~/research/now/sangay/saga_svd_basis_functions_4';
        SAGA_temporaryCatalogFile = '~/research/now/sangay/sagaTmpUpdate';
        SAGA_longTermCatalogFile = '~/research/now/sangay/sangaySubspaceDetectorSAGA_v2';

        if verboseFlag
            disp('processing SAGA sensor, using subspace detector');
        end

        try
            updateCustomCatalog(Sorig(1),SAGA_templateFileName,SAGA_temporaryCatalogFile,...
                SAGA_longTermCatalogFile,0.08,5,150,maxN,diffFlag,20,...
                linearccnorm,plotFlag,verboseFlag);
        catch
            warning('couldnt update SAGA.HHZ');
        end
        toc;

        if verboseFlag
            disp('done SAGA sensor');
        end

        %
        tic;
        Sorig(1) = resampleWaveforms(...
            intWaveforms(...
            filterWaveforms(...
            nanGapWaveforms(...
            taperWaveforms(...
            detrendWaveforms(...
            differentiateWaveforms(Sorig(1))),tw),fillValue),0.25,2)),newFs);
        toc;
        principalSeismic = filterUniqueEvents(SAGA_longTermCatalogFile,20);

        %
        locs1 = principalSeismic >= Sorig(1).ref;
        if sum(locs1)
            principalSeismic = principalSeismic(locs1);
        end
        pause(pauseTime);
        toc;
    else
        Sorig(1) = resampleWaveforms(...
            intWaveforms(...
            filterWaveforms(...
            nanGapWaveforms(...
            taperWaveforms(...
            detrendWaveforms(...
            differentiateWaveforms(Sorig(1))),tw),fillValue),lfc,hfc)),newFs);
        if verboseFlag
            disp('applying sta/lta detector on principal SEISMIC sensor');
        end
        locs1 = stalta(nanGapWaveforms(Sorig(1),0),20,20,2,true,0,true,1/5);

        %
        if ~isempty(locs1)
            tSeismic = getTimeVec(Sorig(1));
            principalSeismic = tSeismic(locs1);
            clear tSeismic;
        end
    end

    %
    if iFlag %if principal has requested infrasound...
        if lS > 1 && strcmpi(Sorig(2).kstnm,principal)
            % infrasound data found...
            Sorig(2) = resampleWaveforms(...
                intWaveforms(...
                filterWaveforms(...
                nanGapWaveforms(...
                taperWaveforms(...
                detrendWaveforms(...
                differentiateWaveforms(Sorig(2))),tw),fillValue),lfc,hfc)),newFs);

            if strcmp(saveDir,'sangay')
                if verboseFlag
                    disp('assuming SAGA infrasound is same as SAGA seismic');
                end
                tPrincipalPossible = principalSeismic;
            else
                if verboseFlag
                    disp('applying sta/lta detector on principal INFRASOUND sensor');
                end

                %
                locs2 = stalta(nanGapWaveforms(Sorig(2),0),10,20,5,true,0,true,1/5);
                if ~isempty(locs2)
                    %% infrasound detections found
                    tInfrasound = getTimeVec(Sorig(2));
                    principalInfrasound = tInfrasound(locs2);
                    clear tInfrasound;

                    %
                    if ~isempty(locs1)
                        %% merge time data and find possible matches
                        tAll = [principalSeismic; principalInfrasound];
                        iAll = [zeros(size(principalSeismic)); ones(size(principalInfrasound))];
                        [tAll,sI] = sort(tAll);
                        iAll = iAll(sI);
                        difft = seconds(diff(tAll));
                        tPrincipalPossible = difft <= 25 & diff(iAll) == 1;
                        tPrincipalPossible = tAll(tPrincipalPossible);
                    end
                end
            end

            %
            if lS > 2
                Sorig(3:lS) = resampleWaveforms(...
                    intWaveforms(...
                    filterWaveforms(...
                    nanGapWaveforms(...
                    taperWaveforms(...
                    detrendWaveforms(...
                    differentiateWaveforms(Sorig(3:lS))),tw),fillValue),lfc,hfc)),newFs);
            end

        elseif lS > 1
            Sorig(2:lS) = resampleWaveforms(...
                intWaveforms(...
                filterWaveforms(...
                nanGapWaveforms(...
                taperWaveforms(...
                detrendWaveforms(...
                differentiateWaveforms(Sorig(2:lS))),tw),fillValue),lfc,hfc)),newFs);
        end
    elseif lS > 1
        % infrasound not requested
        Sorig(2:lS) = resampleWaveforms(...
            intWaveforms(...
            filterWaveforms(...
            nanGapWaveforms(...
            taperWaveforms(...
            detrendWaveforms(...
            differentiateWaveforms(Sorig(2:lS))),tw),fillValue),lfc,hfc)),newFs);
    end
else
    % principal not found, filtering all with same parameters
    Sorig = resampleWaveforms(...
        intWaveforms(...
        filterWaveforms(...
        nanGapWaveforms(...
        taperWaveforms(...
        detrendWaveforms(...
        differentiateWaveforms(Sorig)),tw),fillValue),lfc,hfc)),newFs);
end

%%
kstnms = pull(Sorig,'kstnm');
[stla,stlo] = metaDataFromStationList(kstnms);
if strcmp(saveDir,'sangay')
    d_ = distance(stla,stlo,-2.00535,-78.341294,refEllipse)*1e-3;
elseif strcmp(saveDir,'reventador')
    d_ = distance(stla,stlo,-0.080850,-77.657995,refEllipse)*1e-3;
elseif strcmp(saveDir,'cotopaxi')
    d_ = distance(stla,stlo,-0.683727,-78.436542,refEllipse)*1e-3;
elseif strcmp(saveDir,'fernandina')
    d_ = distance(stla,stlo,-0.371528,-91.539583,refEllipse)*1e-3;
end

%%
[~,sI] = sort(d_);
Sorig = Sorig(sI);

%%
refTimes = pull(Sorig,'ref');
endTimes = refTimes + pull(Sorig,'e');

%%
maxEnd = max(endTimes);
cutStart = maxEnd - seconds(3600*nHours);

%%
Scut = populateWaveforms(lS);
for ii = 1:lS
    Scut(ii) = cutWaveforms(Sorig(ii),cutStart,0,endTimes(ii) - cutStart);
end

%% plot shorter record section
close all;
clipFlag = true;
if clipFlag
    if verboseFlag
        disp('clipping and plotting waveforms...');
    end

    %%
    SclippedPlus = clipWaveforms(Scut,4,true);
    for jj = 1:(iFlag+1)
        SclippedPlus(jj) = Scut(jj);
    end

    %%
    for ii = 1:lS
        sc_ = SclippedPlus(ii);
        if ~isnat(sc_.ref)
            Stmp = populateWaveforms();
            Stmp.kstnm = 'Gamma';
            newDT = median(seconds(diff(t_)));
            Stmp = dealHeader(Stmp,ccnorm_,1./newDT,t_(1) - seconds(recordLength));
            Stmp = cutWaveforms(Stmp,cutStart,0,t_(end) - cutStart,true,0);
            SclippedPlus = [SclippedPlus; Stmp];
            break;
        end
    end

    %%
    visibilityMode = 'off';
    if debuggingMode
        visibilityMode = 'on';
    end

    [~,ax] = plotWaveforms(SclippedPlus,...
        [],[],[],[],[],[],[],...
        visibilityMode);

    hold(ax(end),'on'); %<-- prep for plotting

    if isempty(threshold)
        threshold = median(ccnorm_,"omitnan") + 8*mad(ccnorm_,1);
    end

    plot(ax(end),[SclippedPlus(end).ref SclippedPlus(end).ref+SclippedPlus(end).e],...
        [threshold threshold],'--','Color',[0.5 0.5 0.5],'linewidth',1);

    tI = tPotential >= SclippedPlus(end).ref;
    plot(ax(end),tPotential(tI),NCC(tI),'p','linewidth',1);
    if sum(tI)
        titleStr = ['Ultima Actualizacion: ',datestr(now,31),' (Tiempo Local); Numero de Posibles Eventos Detectados: ',num2str(sum(tI))];
    else
        titleStr = ['Ultima Actualizacion: ',datestr(now,31),' (Tiempo Local)'];
    end

    %%
    jj = 0;
    for ii = 1:lS
        sc_ = SclippedPlus(ii);
        if ~isnat(sc_.ref)
            jj = jj + 1;
            title(ax(jj),titleStr); %,'FontSize',16);
            break;
        end
    end

else
    [~,ax] = plotWaveforms(Scut,...
        [],[],[],[],[],[],[],...
        'off');
end
xlim([cutStart maxEnd]);

%% highlight Principal Detections
if exist('tPrincipalPossible','var') % there are events to highlight
    if strcmp(saveDir,'sangay')
        snippetDuration = 150; % seconds
    else
        snippetDuration = shortSnippet; % seconds
    end

    %%
    snippetDuration = round(snippetDuration*newFs); % samples
    findI = tPrincipalPossible >= cutStart;
    tP_ = tPrincipalPossible(findI);

    %%
    if sum(findI)
        for jj = 1:(iFlag+1)
            dCut = SclippedPlus(jj).d;
            if ~isempty(dCut)
                tLastHour = getTimeVec(Scut(jj));
                lMax = Scut(jj).npts;
                lMax = lMax - snippetDuration + 1;
                for kk = 1:sum(findI)
                    findI_ = find(tLastHour >= tP_(kk),1);
                    if findI_ <= lMax
                        cutI = (findI_:findI_+snippetDuration-1)';
                        hold(ax(jj),'on');
                        plot(ax(jj),tLastHour(cutI),dCut(cutI),'k','Linewidth',1);
                    end
                end
            end
        end
    end
    %clear tPrincipalPossible;
else
    disp('No data from designated principal sensor(s), nothing to highlight');
end

%% count events from regional sensors
if exist('tPotential','var')
    %%
    snippetDuration = (recordLength + 30);
    snippetDuration = round(snippetDuration*newFs);

    %%
    findI = tPotential >= cutStart & tPotential <= cutStart + seconds(nHours*3600 - 600);
    tPot = tPotential(findI);

    %% display some statistics
    nFound = sum(findI);
    if nFound > 1
        medIET = median(seconds(diff(tPot)),"omitnan");
        fprintf('Total number of detected events using proxy sensors within last <strong>%d</strong> hour(s): <strong>%d</strong>\n',...
            nHours,nFound);
        fprintf('Median inter-event waiting time for this time period: <strong>%f</strong> [sec.]\n',...
            medIET);
        fprintf('Median projected daily rate for this time period: <strong>%f</strong> [#/day]\n',...
            86400./median(medIET,"omitnan"));
    elseif nFound
        disp(['lone event found at: ',datestr(tPot(1))])
    else
        disp('no events detected');
        if verboseFlag
            for kk = 1:length(Scut)
                disp(' ');
                disp(Scut(kk).ref);
                disp(Scut(kk).e);
                disp(Scut(kk).ref + Scut(kk).e);
                disp(' ')
            end
        end
    end

    %%
    hStart = 1;
    if strcmp(Scut(1).kstnm,principal)
        hStart = hStart + 1;
        if iFlag
            if strcmp(Scut(2).kstnm,principal)
                hStart = hStart + 1;
            end
        end
    end

    %%
    if nFound
        Nax = hStart;
        for jj = hStart:lS
            dCut = Scut(jj).d;
            if isempty(dCut)
                continue;
            end

            %%
            tLastHour = getTimeVec(Scut(jj));
            lMax = Scut(jj).npts;
            lMax = lMax - snippetDuration + 1;
            for kk = 1:sum(findI)
                %%
                findI_ = find(tLastHour >= tPot(kk),1);
                if findI_ <= lMax
                    %%
                    cutI = (findI_:findI_+snippetDuration-1)';
                    hold(ax(Nax),'on');

                    %%
                    ap = area(ax(Nax),[min(tLastHour(cutI)) max(tLastHour(cutI))],...
                        [max(ax(Nax).YLim) max(ax(Nax).YLim)],min(ax(Nax).YLim),...
                        'FaceColor',ax(Nax).ColorOrder(2,:));
                    ap.FaceAlpha = faceAlpha;
                    ap.LineStyle = 'none';
                end
            end
            Nax = Nax + 1;

        end
    end
end

%% optional debugging break
% if debuggingMode
%     exitStatus = false;
%     return;
% end

%%
if strcmp(saveDir,'sangay')
    fname1 = fullfile('~','public_html','helis','lastHour.jpg');
else
    fname1 = fullfile('~','public_html','helis',saveDir,'lastHour.jpg');
end
print('-djpeg',fname1);
pause(pauseTime);

if ~debuggingMode
    close all;
end
clear Scut;

%% now get individual sensor helicorders
cutStart = maxEnd - seconds(86400); %one day of helicorder data
for ii = 1:lS
    tic;
    S_ = Sorig(ii);

    if isnat(S_.ref)
        continue;
    end
    kstnm_ = S_.kstnm;
    kcmpnm_ = S_.kcmpnm;
    if strcmp(saveDir,'cotopaxi')
        tmpChan = char(kcmpnm_);
        tmpChan = tmpChan(end);
        if tmpChan ~= 'Z'
            continue;
        end
    end

    %%
    S_ = cutWaveforms(S_,cutStart,0,endTimes(ii) - cutStart,verboseFlag);

    %% generate helicorder
    if isnat(S_.ref)
        continue;
    end

    try
        fig = helicorder(S_);
    catch
        continue;
    end
    fig.Visible = 'off';
    if debuggingMode
        fig.Visible = 'on';
    end

    %%
    if strcmp(saveDir,'sangay')
        if strcmpi(kstnm_,"SAG1")
            khole_ = lower(S_.khole);
            fname = fullfile('~','public_html','helis',strcat(char(lower(kstnm_)),'_',char(lower(kcmpnm_)),'_',char(khole_),'_filtered.jpg'));
        else
            fname = fullfile('~','public_html','helis',strcat(char(lower(kstnm_)),'_',char(lower(kcmpnm_)),'_filtered.jpg'));
        end
    else
        fname = fullfile('~','public_html','helis',saveDir,strcat(char(upper(kstnm_)),'_',char(upper(kcmpnm_)),'_filtered.jpg'));
    end

    print(fig,fname,'-djpeg');
    pause(pauseTime);
    if ~debuggingMode
        close all;
    end

    toc;
end

%% sangay specific updates
if strcmp(saveDir,'sangay')
    %%
    if exist('Shf','var')
        kstnm_ = Shf.kstnm;
        kcmpnm_ = Shf.kcmpnm;
        Shf = cutWaveforms(Shf,cutStart,0,Shf.ref + Shf.e - cutStart,verboseFlag);
        fig = helicorder(Shf);
        fig.Visible = 'off';
        fname = fullfile('~','public_html','helis',strcat(char(lower(kstnm_)),'_',char(lower(kcmpnm_)),'_hf_filtered.jpg'));
        print(fig,fname,'-djpeg');
        pause(pauseTime);
        if ~debuggingMode
            close all;
        end
    end

    %%
    try
        updateAllSangayFeatures(true);
    catch
        warning('couldnt update sangay feature set for some reason')
    end
end