clear; close all; clc;
%cd ~/data/fer/

%%
tStart = datetime(2020,01,01); %datetime(2012,11,08);
tEnd = datetime(2020,03,01); 
%tEnd = %dateshift(dn2dt(now) + hours(5),'start','day');
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
    
    %%
    lfc = 5;
    hfc = 10;
    
    %%
    pickRange = 4;
    newFs = 100;
    diffFlag = 0;
    hFlag = true;
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

%%
close all;
figure(); semilogy(tMaster,maxAmpRMS,'.'); zoom on;
figure(); semilogy(tMaster(1:end-1),3600./medfiltSH(seconds(diff(tMaster)),40,true),'.'); zoom on;
grid on
yyaxis right;
plot(tMaster,1:length(tMaster),'.'); zoom on;

%%
% clc;
S = loadWaveforms(datetime(2020,01,14),01,"FER1","BHZ");
S = resampleWaveforms(S,newFs);
S = filterWaveforms(S,lfc,hfc,npoles);
tt = getTimeVec(S);
dd = S.d;
tI = tMaster >= datetime(2020,01,14) & tMaster <= datetime(2020,01,15);
mm = movmax(dd,(2*sta*newFs)+1);
t_ = tMaster(tI);
[~,otherIndex] = min(abs(tt - t_(1)));
otherIndex = otherIndex + cumsum([0; round(100*seconds(diff(t_)))]);

figure(); plot(tt,dd); hold on;

% Scut = cutWaveforms(S(1),tMaster(1:min([ntmp 1000])),-seconds(2),18,0);
% d = double(pull(Scut));

S = loadWaveforms(datetime(2020,01,14),01,"FER2","HHZ");
S = resampleWaveforms(S,newFs);
S = filterWaveforms(S,lfc,hfc,npoles);
tt = getTimeVec(S);
dd = S.d;
plot(tt,dd);
plot(t_,mm(otherIndex).*ones(sum(tI),1),'kp','linewidth',1,'markersize',14); zoom on;

% Scut = cutWaveforms(S(1),tMaster(1:min([ntmp 1000])),-seconds(2),18,0);
% d = double(pull(Scut));
