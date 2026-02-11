% make puyo helicorders
clear; close all; clc;

%%
debugginMode = true; %default should be off;
if ~debugginMode
    cd ~/matlab/
    PWD = pwd;
    PWD = strcat(PWD,'/working/');
    
    setenv('TZ','America/Guayaquil');
    datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
    javaaddpath([PWD,'/IRIS-WS-2.0.15.jar']);
    javaaddpath([PWD,'/matTaup/lib/matTaup.jar']);
    format long g;
    addpath(genpath(PWD),'-end');
    rng('shuffle');
    
    set(groot,'defaultAxesFontSize',18);
    set(groot,'defaultColorbarFontSize',12);
    set(groot,'defaultLineLineWidth',0.1);
    set(groot,'defaultScatterLineWidth',0.1);
    set(groot,'defaultStemLineWidth',0.1);
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultColorbarTickLabelInterpreter','latex');
    set(groot,'defaultTextarrowshapeInterpreter','latex');
    set(groot,'defaultTextboxshapeInterpreter','latex');
    set(groot,'defaultLineMarkerSize',8);
    set(groot,'defaultFigurePaperPositionMode','auto');
    set(groot,'defaultFigureRenderer','Painters');
    set(groot,'defaultFigureColor',[1 1 1]);
else
    nHours = 20;
    nDays = 0;
end

%%
refEllipse = referenceEllipsoid('wgs84');
pauseTime = 3;
lfc = 0.6;
hfc = 1.2;
newFs = 16;
faceAlpha = 0.2;
nHours = 3;
nDays = 1;

%%
% debugginMode = true; %default should be off;
% if debugginMode
%     nHours = 20;
%     nDays = 0;
% end

%%
% sangaySNCLs = ["SAGA" "HHZ" "EC" "";...
%     "SAGA" "BDF" "EC" "01";...
%     "PUYO" "HHZ" "EC" "";...
%     "TAIS" "HHZ" "EC" "";...
%     "PKYU" "HHZ" "EC" "";...
%     "TAMH" "HHZ" "EC" "";...
%     "BPAT" "BHZ" "EC" "";...
%     "BMAS" "BHZ" "EC" "";...
%     "BRTU" "HHZ" "EC" "";...
%     "BULB" "BHZ" "EC" "";...
%     "BRUN" "BHZ" "EC" "";...
%     "BBIL" "BHZ" "EC" "";...
%     "PORT" "HHZ" "EC" "";...
%     "JSCH" "HHZ" "EC" "";...
%     "CHSH" "HHZ" "EC" "";...
%     "PIS1" "HHZ" "EC" "";...
%     "COHC" "HHZ" "EC" "";...
%     "PIAT" "HHZ" "EC" "";...
%     "BOSC" "HHZ" "EC" "";...
%     "BREF" "BHZ" "EC" ""];

sangaySNCLs = ["SAGA" "HHZ" "EC" "";...
    "SAGA" "BDF" "EC" "01";...
    "PUYO" "HHZ" "EC" "";...
    "TAIS" "HHZ" "EC" "";...
    "PKYU" "HHZ" "EC" "";...
    "TAMH" "HHZ" "EC" "";...
    "BPAT" "BHZ" "EC" "";...
    "BMAS" "BHZ" "EC" "";...
    "BRTU" "HHZ" "EC" "";...
    "BULB" "BHZ" "EC" "";...
    "BRUN" "BHZ" "EC" "";...
    "PORT" "HHZ" "EC" "";...
    "JSCH" "HHZ" "EC" "";...
    "COHC" "HHZ" "EC" "";...
    "PIAT" "HHZ" "EC" ""];

%%
lk = size(sangaySNCLs,1);
cd ~/public_html/helis/

%%
while true
    %%
    nowTime = dn2dt(now)+hours(5);
    load('~/research/now/sangay/SangayRegionalAnalysis_v2','tabs');
    lastT = max(tabs)-nDays;
    clear tabs;
    
    %% gen new
    loopStart = dateshift(lastT,'start','day') - nDays;
    loopEnd = dateshift(nowTime,'start','day');
    templateFileName = '~/research/now/sangay/nineTemplatesEighteenSensorsSangay';
    [indiv_events,tabs,pks,templateIndex,maxAmpRMS,pksOrig,saveStds,Neff] = ...
        singleStationRepeaterSearch(4,templateFileName,9,1,150,loopStart,loopEnd,1e5,false,false);
    save('~/research/now/sangay/sangay_eighteenSensors_update','indiv_events','tabs','pks','templateIndex','maxAmpRMS','pksOrig','saveStds','Neff');
    clear indiv_events tabs pks templateIndex maxAmpRMS pksOrig saveStds Neff
    
    %% sync datasets
    load ~/research/now/sangay/SangayRegionalAnalysis_v2.mat
    
    tNew = load('~/research/now/sangay/sangay_eighteenSensors_update','tabs'); tNew = tNew.tabs;
    mar = load('~/research/now/sangay/sangay_eighteenSensors_update','maxAmpRMS'); mar = mar.maxAmpRMS;
    tINew = load('~/research/now/sangay/sangay_eighteenSensors_update','templateIndex'); tINew = tINew.templateIndex;
    stds = load('~/research/now/sangay/sangay_eighteenSensors_update','saveStds'); stds = stds.saveStds;
    pksONew = load('~/research/now/sangay/sangay_eighteenSensors_update','pksOrig'); pksONew = pksONew.pksOrig;
    pksNew = load('~/research/now/sangay/sangay_eighteenSensors_update','pks'); pksNew = pksNew.pks;
    neffNew = load('~/research/now/sangay/sangay_eighteenSensors_update','Neff'); neffNew = neffNew.Neff;
    
    tGood = tNew > lastT;
    tKeep = tabs <= lastT;
    
    maxAmpRMS = [maxAmpRMS(tKeep,:); mar(tGood,:)];
    tabs = [tabs(tKeep); tNew(tGood)];
    templateIndex = [templateIndex(tKeep); tINew(tGood)];
    pksOrig = [pksOrig(tKeep); pksONew(tGood)];
    pks = [pks(tKeep); pksNew(tGood)];
    Neff = [Neff(tKeep); neffNew(tGood)];
    
    %%
    clear tNew mar tINew stds pksONew pksNew neffNew lastT;
    save('~/research/now/sangay/SangayRegionalAnalysis_v2','indiv_events','tabs','pks','templateIndex','maxAmpRMS','pksOrig','saveStds','Neff');
    clear indiv_events tabs pks templateIndex maxAmpRMS pksOrig saveStds Neff
    
    %% read data
    Sorig = populateWaveforms(lk);
    day_ = dateshift(dn2dt(floor(datenum(nowTime))),'start','day');
    [yyyy,mm,dd] = datevec(day_);
    doy = day(day_,'doy');
    
    if doy < 10
        doyStr = ['00',num2str(doy)];
    elseif doy < 100
        doyStr = ['0',num2str(doy)];
    else
        doyStr = num2str(doy);
    end
    
    %%
    n = 0;
    for ii = 1:lk
        SNCL_ = sangaySNCLs(ii,:);
        kstnm_ = SNCL_(1);
        kcmpnm_ = SNCL_(2);
        knetwk_ = SNCL_(3);
        khole_ = SNCL_(4);
        
        f = fullfile('~','rawdata',num2str(yyyy),char(knetwk_),char(kstnm_),strcat(char(kcmpnm_),'.D'));
        miniseedfile = dir(strcat(f,'/',char(knetwk_),'.',char(kstnm_),'.',char(khole_),'.',char(kcmpnm_),'.*.',doyStr));
        if ~isempty(miniseedfile)
            S_ = loadWaveforms(day_-1,2,kstnm_,kcmpnm_,knetwk_,khole_);
            if ~isnat(S_.ref)
                n = n + 1;
                Sorig(n,1) = S_;
            end
        end
    end
    Sorig = Sorig(1:n,1);
    lS = length(Sorig);
    
    %%
    if lS
        %% plot last hour
        if strcmp(Sorig(1).kstnm,"SAGA")
            Sorig(1) = filterWaveforms(Sorig(1),0.6,4.8);
            locs1 = stalta(Sorig(1),20,20,2,true,0,true,1/5);
            if ~isempty(locs1)
                tSeismic = getTimeVec(Sorig(1));
                sagaSeismic = tSeismic(locs1);
                clear tSeismic;
            end
            
            if lS > 1 && strcmp(Sorig(2).kstnm,"SAGA")
                Sorig(2) = filterWaveforms(Sorig(2),0.3,1.2);
                locs2 = stalta(Sorig(2),10,20,5,true,0,true,1/5);
                if ~isempty(locs2)
                    tInfrasound = getTimeVec(Sorig(2));
                    sagaInfrasound = tInfrasound(locs2);
                    clear tInfrasound;
                    
                    if ~isempty(locs1)
                        tAll = [sagaSeismic; sagaInfrasound];
                        iAll = [zeros(size(sagaSeismic)); ones(size(sagaInfrasound))];
                        [tAll,sI] = sort(tAll);
                        iAll = iAll(sI);
                        difft = seconds(diff(tAll));
                        tSagaPossible = difft <= 25 & diff(iAll) == 1;
                        tSagaPossible = tAll(tSagaPossible);
                    end
                end
            end
            
            if lS > 2
                Sorig(3:lS) = filterWaveforms(Sorig(3:lS),lfc,hfc);
            end
        else
            Sorig = filterWaveforms(Sorig,lfc,hfc);
        end
        
        %%
        Sorig = resampleWaveforms(Sorig,newFs);
        
        %%
        kstnms = pull(Sorig,'kstnm');
        [stla,stlo] = metaDataFromStationList(kstnms);
        d_ = distance(stla,stlo,-2.0053,-78.34078,refEllipse)*1e-3;
        [~,sI] = sort(d_);
        
        d_ = d_(sI);
        Sorig = Sorig(sI);
        
        refTimes = pull(Sorig,'ref');
        endTimes = refTimes + pull(Sorig,'e');
        
        minStart = min(refTimes);
        maxEnd = max(endTimes);
        
        cutStart = maxEnd - seconds(3600*nHours);
        Scut = populateWaveforms(lS);
        for ii = 1:lS
            Scut(ii) = cutWaveforms(Sorig(ii),cutStart,0,endTimes(ii) - cutStart);
        end
        
        %%
        close all;
        if debugginMode
            [~,ax] = plotWaveforms(Scut,...
                [],[],[],[],[],[],[],...
                'on');
        else
            [~,ax] = plotWaveforms(Scut,...
                [],[],[],[],[],[],[],...
                'off');
        end
        xlim([cutStart maxEnd]);
        
        %% highlight SAGA
        if exist('tSagaPossible','var') %SAGAHHZ and SAGABDF exist, and there are events to highlight
            snippetDuration = 60; % seconds
            snippetDuration = round(snippetDuration*newFs); % samples
            findI = tSagaPossible >= cutStart;
            tSaga_ = tSagaPossible(findI);
            
            if sum(findI)
                for jj = 1:2
                    dCut = Scut(jj).d;
                    tLastHour = getTimeVec(Scut(jj));
                    lMax = Scut(jj).npts;
                    lMax = lMax - snippetDuration + 1;
                    for kk = 1:sum(findI)
                        findI_ = find(tLastHour >= tSaga_(kk),1);
                        if findI_ <= lMax
                            cutI = (findI_:findI_+snippetDuration-1)';
                            hold(ax(jj),'on');
                            plot(ax(jj),tLastHour(cutI),dCut(cutI),'k','Linewidth',1);
                        end
                    end
                end
            end
            clear tSagaPossible;
        end
        
        %% highlight Regional (using BB regional sensors)
        tPotential = sangayRegionalUniqueEvents_v2();
        if exist('tPotential','var') %SAGAHHZ and SAGABDF exist, and there are events to highlight
            snippetDuration = 150; % seconds
            snippetDuration = round(snippetDuration*newFs); % samples
            findI = tPotential >= cutStart & tPotential <= cutStart + seconds(nHours*3600 - 600);
            tPot = tPotential(findI);
            
            hStart = 1;
            if strcmp(Scut(1).kstnm,"SAGA")
                hStart = hStart + 1;
                if strcmp(Scut(2).kstnm,"SAGA")
                    hStart = hStart + 1;
                end
            end
            
            if sum(findI)
                Nax = hStart;
                for jj = hStart:lS
                    dCut = Scut(jj).d;
                    if ~isempty(dCut)
                        tLastHour = getTimeVec(Scut(jj));
                        lMax = Scut(jj).npts;
                        lMax = lMax - snippetDuration + 1;
                        for kk = 1:sum(findI)
                            findI_ = find(tLastHour >= tPot(kk),1);
                            if findI_ <= lMax
                                cutI = (findI_:findI_+snippetDuration-1)';
                                hold(ax(Nax),'on');
                                %plot(ax(Nax),tLastHour(cutI),dCut(cutI),'Linewidth',2,'Color',ax(Nax).ColorOrder(2,:));
                                
                                %%
                                ap = area(ax(Nax),[min(tLastHour(cutI)) max(tLastHour(cutI))],...
                                    [max(ax(Nax).YLim) max(ax(Nax).YLim)],min(ax(Nax).YLim),...
                                    'FaceColor',ax(Nax).ColorOrder(2,:));
                                ap.FaceAlpha = faceAlpha;
                                ap.LineStyle = 'none';
                            end
                        end
                        Nax = Nax + 1;
                    else
                        continue;
                    end
                end
            end
            clear tPotential;
        end
        
        %% optional debugging break
        if debugginMode
            break;
        end
        
        %%
        fname1 = fullfile('~','public_html','helis','lastHour.jpg');
        print('-djpeg',fname1);
        pause(pauseTime);
        close all;
        clear Scut;
        
        %% now get individual sensor helicorders
        cutStart = maxEnd - seconds(86400);
        for i = 1:lS
            S_ = Sorig(i);
            if ~isnat(S_.ref)
                kstnm_ = S_.kstnm;
                kcmpnm_ = S_.kcmpnm;
                
                %%
                S_ = cutWaveforms(S_,cutStart,0,endTimes(i) - cutStart);
                
                %% generate helicorder
                if ~isnat(S_.ref)
                    [fig,ax] = helicorder(S_);
                    fig.Visible = 'off';
                    
                    %%
                    fname = fullfile('~','public_html','helis',strcat(lower(kstnm_),"_",lower(kcmpnm_),"_filtered.jpg"));
                    print('-djpeg',fname);
                    pause(pauseTime);
                    close all;
                end
            end
        end
    end
end
