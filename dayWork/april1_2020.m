% april1_2020
clear; close all; clc;

%%
refEllipse = referenceEllipsoid('wgs84');
center_lon = -77.6581;
center_lat = -0.0810;

%%
[stla,stlo,stel] = metaDataFromStationList(["REVN","LAV4","REVS"]);
d_ = distance(stla,stlo,center_lat,center_lon,refEllipse)*1e-3;

%%
cd ~/products/rsam/;
load EC.REVN..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat
X = S;
load EC.LAV4..SHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat
X(2,1) = S;
%load EC.REVS..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat
load EC.REVS.01.BDF_4Sec1Hz_60DUR_Pa_RMS_SLIDEMED_MedPreFilt3Points.mat
X(3,1) = S;
S = X;
clear X;

%%
%log10(1) + 1.110*log10(d_) + 0.00189*d_ + 0.591
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
lTimes = length(tStart);

%%
durs = days(tEnd - tStart);

%%
lS = length(S);
nCount = NaN(lTimes,lS);
completeness = nCount;
energy = nCount;
energy2 = energy;
energy3 = energy;
maxAmp = cell(size(nCount));
t = maxAmp;

for i = 1:lS
    disp(i);
    S_ = S(i);
    t_ = getTimeVec(S_);
    
    [locs_,snr,~,sosSTA_] = stalta(S_,2*60,5*60,1.75,false,0,false);
    if i == 3
        locs_ = locs_(snr <= 100);
    end
    %locs{i} = locs_;
    %snr{i} = snr_;
    sosSTA_ = sosSTA_(locs_);
    t_ = t_(locs_);
    
    if i == 1 %REVN
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
                    sosSTA_(tII) = sosSTA_(tII) / respSensitivity(ii);
                end
            end
        else
            sosSTA_ = 1e6*sosSTA_/respSensitivity;
        end
    elseif i == 2 %LAV4
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
                    sosSTA_(tII) = sosSTA_(tII) / respSensitivity(ii);
                end
            end
        else
            sosSTA_ = 1e6*sosSTA_/respSensitivity;
        end
    else
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
                    sosSTA_(tII) = sosSTA_(tII) / respSensitivity(ii);
                end
            end
        else
            sosSTA_ = 1e6*sosSTA_/respSensitivity;
        end
    end
    
    
    for j = 1:lTimes
        disp(j);
        tStart_ = tStart(j);
        tEnd_ = tEnd(j);
        dur_ = durs(j); %total gap in days
        
        %
        disp('counting number of events');
        goodI = t_ >= tStart_ & t_ <= tEnd_;
        nCount(j,i) = sum(goodI);
        maxAmp{j,i} = sosSTA_(goodI);
        t{j,i} = t_(goodI);
        energy2(j,i) = sum(sosSTA_(goodI));
        energy3(j,i) = nanmedian(sosSTA_(goodI));
        disp(nCount(j,i))
        
        %
        disp('cutting window');
        Scut = cutWaveforms(S_,tStart_,0,86400*dur_);
        
        %
        disp('energy');
        energy(j,i) = nansum(Scut.d);
        
        %
        gaps = Scut.gapInfo;
        if ~isempty(gaps)
            sumGaps = sum(gaps(:,2))/1440; %gap in days
            completeness(j,i) = 1-(sumGaps/dur_);
        else
            completeness(j,i) = 1;
        end
        disp(completeness(j,i));
        
    end
end

%%
scaledEnergyDensity = 1e8*(((energy2./completeness)./durs)./(nCount./completeness));
% for i = 1:3
% scaledEnergyDensity(:,i) = scaledEnergyDensity(:,i)/min(scaledEnergyDensity(:,i));
% end

close all; figure(); plot(scaledEnergyDensity,'o-'); zoom on;
ax = gca;
ax.YScale = 'log';
legend('REVN','LAV4','REVS');

%%
titleLabels = ["REVN","LAV4","REVS"];
maxN = NaN(size(maxAmp)); 
%maxNorig = maxN;
for i = 1:3
disp(i);
for j = 1:7
maxN(j,i) = length(maxAmp{j,i});
end
end
maxNorig = maxN; maxN = max(maxN);

for i = 1:3
disp(i);
figure(1);
ax(i) = subplot(3,1,i); hold on;
nowData = NaN(maxN(i),7);
for j = 1:7
data_ = maxAmp{j,i};
nowData(1:length(data_),j) = data_;
end
boxplot(log10(1e6*nowData),num2str(maxNorig(:,i)),'Notch','on');
title(titleLabels(i));
ylabel('pseudo-mag');
end

%%
tt_ = t{2,3};
data_ = maxAmp{2,3};
figure(); 
plot(tt_,log10(1e6*data_),'o'); zoom on; grid on;





