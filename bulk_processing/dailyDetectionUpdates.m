function dailyDetectionUpdates
cd ~/
clear; close all;

%%
saveFlag = false;
updateCatalogs();
updateSangay(saveFlag);
[Sall,E,goodI] = loadEventWaveforms(fnames,varargin);

end

%%
function updateCatalogs
updateGlobalCatalog();
%%
[t,eqlat,eqlon,eqdepth,eqmag,id] = readCat1();
end

function plotEcuadorSeismicity
end

function 

end

function generateAnimations
end

function updateSangay(saveFlag)
%% Sangay (using PUYO)
load('~/research/now/sangay/puyo_match_filter_analysis');
thresh = 0.6;
diffFlag = true;
lastDay = dn2dt(floor(datenum(max(tabs))));
tStart = lastDay;
tEnd = dn2dt(floor(now));
[~,tabsNew,pksNew,templateIndexNew,maxAmpRMSNew,pksOrigNew] = ...
    singleStationRepeaterSearch(thresh,'~/research/now/sangay/puyo_template.mat',1,1,60,tStart,tEnd,1e4,diffFlag);

% clobber redundant
tI = tabs >= lastDay;
tabs(tI) = [];
maxAmpRMS(tI) = [];
pks(tI) = [];
pksOrig(tI) = [];
templateIndex(tI) = [];

% append (fill with) new
tabs = [tabs; tabsNew];
maxAmpRMS = [maxAmpRMS; maxAmpRMSNew];
pks = [pks; pksNew];
pksOrig = [pksOrig; pksOrigNew];
templateIndex = [templateIndex; templateIndexNew];
clearvars -except tabs maxAmpRMS pks pksOrig saveStds templateIndex
if saveFlag
save('~/research/now/sangay/puyo_match_filter_analysis')
end
end

function updateFernandina
tStart = datetime(2012,01,01); %datetime(2012,11,08);
tEnd = dn2dt(floor(now));
days = tStart:tEnd;
tMaster = NaT(1e6,1);
maxAmpRMS = NaN(1e6,1);
ntmp = 1;

for i = 1:length(days)
    day_ = days(i);
    disp(' ');
    disp(day_);
    disp(' ');
    
    %%
    mph = 2;
    sta = 5;
    lta = 5;
    lfc = 5;
    hfc = 10;
    pickRange = 4;
    newFs = 100;
    diffFlag = 0;
    hFlag = 0;
    staltaPlotFlag = false;
    tw = 0.002;
    minSeparation = 2*sta;
    envFlag = 0;
    envHFC = 10;
    npoles = 8;
    
    %%
    S = loadWaveforms(day_,1,"FER2","HHZ");
    S(2,1) = loadWaveforms(day_,1,"FER1","BHZ");
    if isnat(S(1).ref) || isnat(S(2).ref)
        disp(' ');
        fprintf('skipping: %s\n',datestr(day_));
        disp(' ');
    else
        if diffFlag
            S = taperWaveforms(differentiateWaveforms(S),tw);
        end
        S = filterWaveforms(S,lfc,hfc,npoles);
        S = resampleWaveforms(S,newFs);
        
        %%
        locsFER2 = stalta(S(1),sta,lta,mph,hFlag,staltaPlotFlag,envFlag,envHFC);
        tFER2 = getTimeVec(S(1));
        tFER2 = tFER2(locsFER2);
        
        locsFER1 = stalta(S(2),sta,lta,mph,hFlag,staltaPlotFlag,envFlag,envHFC);
        tFER1 = getTimeVec(S(2));
        tFER1 = tFER1(locsFER1);
        
        %%
        FER1 = zeros(length(tFER1),1);
        FER2 = ones(length(tFER2),1);
        
        %%
        indices = [FER1; FER2];
        t = [tFER1; tFER2];
        
        %%
        [t,sI] = sort(t);
        indices = indices(sI);
        
        %%
        minT = min(t);
        t = seconds(t - minT);
        
        %%
        difft = diff(t);
        diffIndices = abs(diff(indices));
        
        %%
        I = (difft < pickRange) & (diffIndices > 0);
        sumI = sum(I);
        
        %%
        if sumI
            disp(['Number of Events Found: ',num2str(sumI)]);
            
            t = t(1:end-1);
            indices = indices(1:end-1);
            t = t(I);
            indices = indices(I);
            
            %%
            if sumI > 1
                I = find(diff(t) < minSeparation);
                t(I) = [];
                indices(I) = [];
            end
            
            sumI = sumI - length(find(I));
            newT = minT + seconds(t);
            Scut = cutWaveforms(S(2),newT,0,2*sta,0);
            d = double(pull(Scut));
            maxAmp = rms(detrend(diff(d)))';
            
            tMaster(ntmp:ntmp+sumI-1) = newT;
            maxAmpRMS(ntmp:ntmp+sumI-1) = maxAmp;
            ntmp = ntmp + sumI;
            
        end
    end
end
ntmp = ntmp - 1;
tMaster = tMaster(1:ntmp);
maxAmpRMS = maxAmpRMS(1:ntmp);
end