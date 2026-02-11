clear;
close all;

%% preallocate
maxN = 10000;
detection_times = NaT(maxN,1);
fund_freq = NaN(maxN,1);

%%
harm_freq_1 = fund_freq;
harm_freq_2 = fund_freq;
harm_freq_3 = fund_freq;

%%
HSI_1 = fund_freq;
HSI_2 = fund_freq;
HSI_3 = fund_freq;

%%
% tStart = datetime(2009,01,01);
% tEnd = dn2dt(ceil(now));

%tStart = datetime(2015,08,13);
%tEnd = datetime(2015,08,13);

%tStart = datetime(2021,01,05);
%tEnd = datetime(2021,01,05);

%tStart = datetime(2012,12,17);
%tEnd = datetime(2012,12,17);

%tStart = datetime(2009,10,31); %06h44, 20h38, 23h01-03 %<-- dud, due to gaps, not real
%tEnd = datetime(2009,10,31);

% tStart = datetime(2008,11,16); %lots on this day
% tEnd = datetime(2008,11,16);

%tStart = datetime(2010,05,28); %16h29-16h30 %day of tungu eruption?
%tEnd = datetime(2010,05,28);

%tStart = datetime(2015,05,28);
%tEnd = datetime(2015,05,28);

%tStart = datetime(2015,07,26);
%tEnd = datetime(2015,07,26);

% tStart = datetime(2021,05,12);
% tEnd = datetime(2021,05,15);

%tStart = datetime(2021,01,01);
%tEnd = datetime(2021,03,31);

%tStart = datetime(2006,01,01);
%tEnd = datetime(2021,04,01);

%tEnd = datetime(2012,01,03);

% tStart = datetime(2021,05,18);
% tEnd = datetime(2021,05,22);

% tStart = datetime(2017,06,01);
% tEnd = datetime(2017,07,01);
% 
% tStart = datetime(2017,06,27);
% tEnd = datetime(2017,06,27);

%tStart = datetime(2020,12,16);
%tEnd = datetime(2020,12,16);

% tStart = datetime(2020,11,23);
% tEnd = datetime(2020,11,23);

tStart = datetime(2025,04,14);
tEnd = tStart;

%%
days = (tStart:tEnd)';
ldays = length(days);
minHSI = 30;
plotFlag = 0;
lfc = 0.25; %for Reventador, might be better to go down to 0.4;
%lfc = 0.5; %for Reventador, might be better to go down to 0.4;
hfc = 64*lfc;
newFs = 50;
if newFs < hfc*2
    return
end
welchFlag = ~true;
maxFundFreq = 2.5;

%%
n = 1;
for i = 1:ldays
    %S = loadWaveforms(days(i),1,"BREF","BHZ");
    %S = loadWaveforms(days(i),1,"REVN","HHZ");
    %S = loadWaveforms(days(i),1,"REVS","HHZ");
    %S = loadWaveforms(days(i),1,"CONE","SHZ");
    S = loadWaveforms(days(i),1,"LAV4","SHZ");
    
    %
    if isnat(S.ref)
        disp(['no data, skipping: ',datestr(days(i))]);
        continue;
    end
    
    %
    if isempty(S)
        disp(['is empty, skipping: ',datestr(days(i))]);
        continue;
    end
    
    disp(['processing: ',datestr(days(i))]);
    %S_ = transferWaveforms(S,0.5,20,4,50,'disp',true,true);
    %S_ = transferWaveforms(S,lfc,hfc,4,50,'vel');
    
    S = resampleWaveforms(S,newFs);
    S = filterWaveforms(S,lfc,hfc);
    
    %%
    try
        [detection_times_,...
            fund_freq_,...
            harm_freq_1_,...
            harm_freq_2_,...
            harm_freq_3_,...
            HSI_1_,...
            HSI_2_,...
            HSI_3_] = ...
            tremometer_control(S,plotFlag,minHSI,lfc,welchFlag);
        
        %%
        ldetections = length(detection_times_);
        if ldetections
            detection_times(n:n+ldetections-1) = detection_times_;
            fund_freq(n:n+ldetections-1) = fund_freq_;
            harm_freq_1(n:n+ldetections-1) = harm_freq_1_;
            harm_freq_2(n:n+ldetections-1) = harm_freq_2_;
            harm_freq_3(n:n+ldetections-1) = harm_freq_3_;
            
            %%
            HSI_1(n:n+ldetections-1) = HSI_1_;
            HSI_2(n:n+ldetections-1) = HSI_2_;
            HSI_3(n:n+ldetections-1) = HSI_3_;
            n = n + ldetections;
        end
    catch
        warning('day %s bad, skipping...',datestr(days(i)))
    end
end

%%
if n > 1
    n = n - 1;
    detection_times = detection_times(1:n);
    fund_freq = fund_freq(1:n);
    harm_freq_1 = harm_freq_1(1:n);
    harm_freq_2 = harm_freq_2(1:n);
    harm_freq_3 = harm_freq_3(1:n);
    HSI_1 = HSI_1(1:n);
    HSI_2 = HSI_2(1:n);
    HSI_3 = HSI_3(1:n);
    
    %%
    minHSIs = min([HSI_1 HSI_2 HSI_3],[],2);
    minHSIs = median([HSI_1 HSI_2 HSI_3],2,'omitnan');
    
    %%
    figure();
    fI = fund_freq <= maxFundFreq;
    ax = subplot(311);
    plot(detection_times(fI),1:sum(fI),'.'); zoom on; grid on;
    ax(2) = subplot(312);
    plot(detection_times(fI),fund_freq(fI),'.'); zoom on; grid on;
    ax(3) = subplot(313);
    semilogy(detection_times(fI),minHSIs(fI),'.'); zoom on; grid on;
    linkaxes(ax,'x');
end

%%
fprintf('done\n');