%july29_2020
clear; close all; clc;


%%
refEllipse = referenceEllipsoid('wgs84');
center_lon = -77.6581;
center_lat = -0.0810;

%%
%listOfSensors = ["REVN","LAV4","REVS","CASC"]; %if you uncomment, make
%sure to uncomment the rsam curve below as well

listOfSensors = ["REVN","REVS","CASC"];
lS = length(listOfSensors);
[stla,stlo,stel] = metaDataFromStationList(listOfSensors);
d_ = distance(stla,stlo,center_lat,center_lon,refEllipse)*1e-3; %distance in km

%% load RSAM for future identification of gaps
cd ~/products/rsam/;
load EC.REVN..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat
X = S;
% load EC.LAV4..SHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat
% X = [X;S];
load EC.REVS.01.BDF_4Sec1Hz_60DUR_Pa_RMS_SLIDEMED_MedPreFilt3Points.mat
X = [X;S];
load EC.CASC..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat
X = [X;S];

%%
S = X;
clear X;

%%
cd ~/research/now/reventador/
%log10(1) + 1.110*log10(d_) + 0.00189*d_ + 0.591

%%
revn_mag_correction = 1.23547175560585;
lav4_mag_correction = 1.18706049708942;
revs_mag_correction = 1.37775559719201;

%%
REVN_ResponseInfo(1).tstart = datetime(2012,01,257);
REVN_ResponseInfo(1).sensitivity = 3.141950e+08;
REVN_ResponseInfo(2).tstart = datetime(2015,01,079);
REVN_ResponseInfo(2).sensitivity = 3.141950e+08;
REVN_ResponseInfo(3).tstart = datetime(2019,01,256);
REVN_ResponseInfo(3).sensitivity = 1.256780e+09;

%%
LAV4_ResponseInfo(1).tstart = datetime(2008,01,226);
LAV4_ResponseInfo(1).sensitivity = 6.000000e+07;
LAV4_ResponseInfo(2).tstart = datetime(2014,01,079);
LAV4_ResponseInfo(2).sensitivity = 1.174310e+07;
LAV4_ResponseInfo(3).tstart = datetime(2018,01,004);
LAV4_ResponseInfo(3).sensitivity = 7.339440e+08;

%%
REVS_ResponseInfo(1).tstart = datetime(2012,01,257);
REVS_ResponseInfo(1).sensitivity = 3.141950e+08;
REVS_ResponseInfo(2).tstart = datetime(2015,01,105);
REVS_ResponseInfo(2).sensitivity = 3.141950e+08;

%%
CASC_ResponseInfo(1).tstart = datetime(2012,07,12);
CASC_ResponseInfo(1).sensitivity = 3.141950e+08;

%%
tStart = [datetime(2014,04,01); ...
    datetime(2016,11,10);...
    datetime(2017,09,23);...
    datetime(2018,04,02);...
    datetime(2018,07,11);...
    datetime(2018,09,19);...
    datetime(2018,11,24)];

tEnd = [datetime(2014,07,23);
    datetime(2017,09,22);...
    datetime(2018,04,01);...
    datetime(2018,07,10);...
    datetime(2018,09,18);...
    datetime(2018,11,23);...
    datetime(2019,03,27)];

%%
lTimes = length(tStart);

%%
durs = days(tEnd - tStart);

%%

nCount = NaN(lTimes,lS);
completeness = nCount;
energy = nCount;
energy2 = energy;
energy3 = energy;
maxAmp = cell(size(nCount));
t = maxAmp;

%%
f = 1/3;

%%
for i = 1:lS
    disp(i);
    S_ = S(i);
    kstnm_ = listOfSensors(i);
    
    %% loop through periods with distinct sensor sensitivities
    if strcmp(kstnm_,"REVN") %REVN
        C = load('REVN_HHZ_2');
        t_ = C.t;
        p2p_ = C.p2p;
        clear C;
        
        %%
        respStart = pull(REVN_ResponseInfo,'tstart');
        respSensitivity = pull(REVN_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = t_ >= respStart(ii) & t_ <= respStart(ii+1);
                else
                    tII = t_ >= respStart(ii);
                end
                if sum(tII)
                    p2p_(tII) = p2p_(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_ = p2p_/respSensitivity;
        end
    elseif strcmp(kstnm_,"LAV4") %LAV4
        C = load('LAV4_SHZ_2');
        t_ = C.t;
        p2p_ = C.p2p;
        clear C;
        
        %%
        respStart = pull(LAV4_ResponseInfo,'tstart');
        respSensitivity = pull(LAV4_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = t_ >= respStart(ii) & t_ <= respStart(ii+1);
                else
                    tII = t_ >= respStart(ii);
                end
                if sum(tII)
                    p2p_(tII) = p2p_(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_ = p2p_/respSensitivity;
        end
    elseif strcmp(kstnm_,"REVS") %REVS
        C = load('REVS_HHZ_2');
        t_ = C.t;
        p2p_ = C.p2p;
        clear C;
        
        %%
        respStart = pull(REVS_ResponseInfo,'tstart');
        respSensitivity = pull(REVS_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = t_ >= respStart(ii) & t_ <= respStart(ii+1);
                else
                    tII = t_ >= respStart(ii);
                end
                if sum(tII)
                    p2p_(tII) = p2p_(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_ = p2p_/respSensitivity;
        end
    elseif strcmp(kstnm_,"CASC") %CASC
        C = load('CASC_HHZ_2');
        t_ = C.t;
        p2p_ = C.p2p;
        clear C;
        
        %%
        respStart = pull(CASC_ResponseInfo,'tstart');
        respSensitivity = pull(CASC_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = t_ >= respStart(ii) & t_ <= respStart(ii+1);
                else
                    tII = t_ >= respStart(ii);
                end
                if sum(tII)
                    p2p_(tII) = p2p_(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_ = p2p_/respSensitivity;
        end
    end
    
    %% loopp through predifined ash sampling periods
    for j = 1:lTimes
        disp(j);
        tStart_ = tStart(j);
        tEnd_ = tEnd(j);
        dur_ = durs(j); %total gap in days
        catStart = t_(1);
        %
        disp('counting number of events');
        goodI = t_ >= tStart_ & t_ <= tEnd_;
        sumGood = sum(goodI);
        if sumGood
            nCount(j,i) = sum(goodI);
            %maxAmp{j,i} = p2p_(goodI);
            t{j,i} = t_(goodI);
            
            %%
            B = pi*f/(50*2);
            magEst = log10(p2p_);
            magEst = magEst + log10(d_(i));       %body wave model
            magEst = magEst + B*d_(i)*log10(exp(1)) + 7;
            
            %magEst = log10(p2p_) + 1.110*log10(d_(i)) + 0.00189*d_(i) + 7;
            
            %%
            magEst_ = magEst(goodI & magEst >= nanmedian(magEst) & magEst <= nanmedian(magEst)+1);
            maxAmp{j,i} = magEst(goodI);
            
            %%
            %energy2(j,i) = sum(p2p_(goodI));
            energy2(j,i) = log10(sum(10.^(1.5*magEst_)));
            energy3(j,i) = mean(magEst_);
            disp(nCount(j,i))
            
            %
            disp('cutting window');
            Scut = cutWaveforms(S_,tStart_,0,86400*dur_);
            tRSAM = getTimeVec(Scut);
            
            %
            disp('energy');
            Sd = Scut.d;
            
            if catStart > tStart_
                energy(j,i) = nansum(Sd(tRSAM >= catStart));
            else
                energy(j,i) = nansum(Sd);
            end
            
            %
            gapFlag = Scut.gapFlag;
            if gapFlag
                gaps = Scut.gapInfo;
                gapStartIndex = gaps(:,1);
                
                tGAP = tRSAM(gapStartIndex);
                gapsAfterCatStart = tGAP >= tStart_;
                sumGaps = sum(gaps(gapsAfterCatStart,2))/1440; %gap in days
                if catStart > tStart_
                    sumGaps = sumGaps + days(catStart - tStart_);
                end
                completeness(j,i) = 1-(sumGaps/dur_);
            else
                completeness(j,i) = 1;
            end
            disp(completeness(j,i));
        end
    end
end

%%
%scaledEnergyDensity = (((energy2./completeness)./durs)./(nCount./completeness)); %energy per day per event
scaledEnergyDensity = energy2./durs./nCount; %energy per day per event

% for i = 1:3
% scaledEnergyDensity(:,i) = scaledEnergyDensity(:,i)/min(scaledEnergyDensity(:,i));
% end

close all;
figure();
plot(scaledEnergyDensity,'o-'); zoom on;
legend(listOfSensors);

%%
titleLabels = listOfSensors; %["REVN","LAV4","REVS"];
maxN = NaN(size(maxAmp));

%%
for i = 1:lS
    disp(i);
    for j = 1:7
        maxN(j,i) = length(maxAmp{j,i});
    end
end
maxNorig = maxN;
maxN = max(maxN);

%%
fig = figure();
for i = 1:lS
    disp(i);
    figure(fig);
    ax(i) = subplot(lS,1,i); hold on;
    nowData = NaN(maxN(i),7);
    for j = 1:length(tStart)
        data_ = maxAmp{j,i};
        nowData(1:length(data_),j) = data_;
    end
    boxplot(nowData,num2str(maxNorig(:,i)),'Notch','on');
    title(titleLabels(i));
    ylabel('pseudo-mag');
end
zoom on;

%%
[r,c] = size(t);
figure('units','normalized','outerposition',[0 0 1 1]);
n = 0;
for i = 1:r
    tStart_ = tStart(i);
    tEnd_ = tEnd(i);
    for j = 1:c
        %%
        n = n+1;
        tt_ = t{i,j};
        data_ = maxAmp{i,j};
        
        %%
        if ~isempty(tt_)
            ax_(n) = subplot(r,c,n);
            plot(tt_,data_,'o'); zoom on; grid on;
            xlim([tStart_ tEnd_]);
        end
    end
    
    %%
    if n > c
        linkaxes(ax_(n-c+1:n),'xy');
    end
end

%%
if n > c
    linkaxes(ax_(n-c+1:n),'xy');
end