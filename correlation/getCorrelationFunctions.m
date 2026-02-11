function [dayData,nIndependentChannels,t,lags,maxWindows,lDays,causalLength,...
    tEnd,checkLS,interstation_distance,sameSensorFlag,refTrace,newFs,H1_,H2_,...
    Dorig,lfc,hfc,Nsmooth,nfft] = getCorrelationFunctions(varargin)
%
% [dayData,nIndependentChannels,t,lags,maxWindows,lDays,causalLength,...
%    tEnd,checkLS,interstation_distance,sameSensorFlag] = getCorrelationFunctions(tStart,tEnd,...
%   SNCL1,SNCL2,dW,lfc,hfc,newFs,...
%   secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory);
%

%
% Written by Stephen Hernandez, Instituto Geofisico, Quito, Ecuador
% v1 by Stephen Hernandez: 24 Feb 2021
% 20JUN2025: added decon version for single station cross-component case

%%
nVarargin = length(varargin);

% functionDefaults = {...
%     datetime(2023,01,131),...       % tStart
%     datetime(2023,01,131),...       % tEnd
%     ["CTX1","DPZ","EC",""],...      % SNCL1
%     ["CTX2","DPZ","EC",""],...      % SNCL2
%     false,...                       % allFlag
%     0,...                           % dW
%     1,...                           % lfc
%     2,...                           % hfc
%     25,...                          % newFs
%     3600,...                        % secDur
%     10,...                          % maxLag
%     1,...                           % Nsmooth
%     false,...                       % transferFlag
%     true,...                        % verboseFlag
%     "~/masa/backups"};  % pathToWaveformServerOrDirectory

functionDefaults = {...
    datetime(2016,04,01),...        % tStart
    datetime(2016,04,30),...        % tEnd
    ["VCH1","HHZ","EC",""],...      % SNCL1
    ["VCH1","HHZ","EC",""],...      % SNCL2
    true,...                        % allFlag
    0,...                           % dW
    2,...                           % lfc
    4,...                           % hfc
    40,...                          % newFs
    1800,...                        % secDur
    2^04,...                        % maxLag
    1,...                           % Nsmooth
    false,...                       % transferFlag
    false,...                       % verboseFlag
    "~/masa/backups"};              % pathToWaveformServerOrDirectory

% % nice results
% functionDefaults = {...
%     datetime(2018,09,01),...        % tStart
%     datetime(2018,10,01),...        % tEnd
%     ["SN14","HHZ","9D",""],...      % SNCL1
%     ["SN11","HHZ","9D",""],...      % SNCL2
%     true,...                        % allFlag
%     0,...                           % dW
%     0.5,...                         % lfc
%     1.5,...                         % hfc
%     32,...                          % newFs
%     2^11,...                        % secDur
%     2^05,...                        % maxLag
%     1,...                           % Nsmooth
%     true,...                        % transferFlag
%     true,...                        % verboseFlag
%     '/Users/stephen/data/galapagos/iguana/BROADBAND'};  % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2014,01,001),...       % tStart
%     datetime(2015,01,001),...       % tEnd
%     ["BREF","BHZ","EC",""],...      % SNCL1
%     ["BTAM","BHZ","EC",""],...      % SNCL2
%     false,...                       % allFlag
%     0,...                           % dW
%     1,...                           % lfc
%     2,...                           % hfc
%     2^4,...                         % newFs
%     2^10,...                        % secDur
%     2^05,...                        % maxLag
%     84,...                          % Nsmooth
%     false,...                       % transferFlag
%     true,...                        % verboseFlag
%     '~/rawdata'};                   % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2019,01,001),...       % tStart
%     datetime(2021,01,001),...       % tEnd
%     ["FER1","BHZ","EC",""],...      % SNCL1
%     ["VCH1","HHZ","EC",""],...      % SNCL2
%     false,...                       % allFlag
%     0,...                           % dW
%     1/32,...                        % lfc
%     1,...                           % hfc
%     2^3,...                         % newFs
%     2^12,...                        % secDur
%     2^07,...                        % maxLag
%     21,...                          % Nsmooth
%     false,...                       % transferFlag
%     true,...                        % verboseFlag
%     '~/rawdata'};                   % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2023,01,270),...        % tStart
%     datetime(2023,01,300),...       % tEnd
%     ["SN12","HHZ","EC",""],...      % SNCL1
%     ["SN07","HHZ","EC",""],...      % SNCL2
%     false,...                       % allFlag
%     0,...                           % dW
%     0.125,...                       % lfc
%     4,...                           % hfc
%     2^5,...                         % newFs
%     2^7,...                         % secDur
%     2^06,...                        % maxLag
%     675,...                         % Nsmooth
%     false,...                       % transferFlag
%     true,...                        % verboseFlag
%     '~/data/ForAndy/'};  % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2023,01,01),...        % tStart
%     datetime(2023,01,311),...        % tEnd
%     ["SN12","HHN","EC",""],...      % SNCL1
%     ["SN12","HHE","EC",""],...      % SNCL2
%     ~true,...                       % allFlag
%     0,...                           % dW
%     1,...                           % lfc
%     2,...                           % hfc
%     2^4,...                         % newFs
%     2^11,...                        % secDur
%     2^05,...                        % maxLag
%     [],...                          % Nsmooth
%     false,...                       % transferFlag
%     true,...                        % verboseFlag
%     '~/data/ForAndy/'};             % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2019,10,01),...        % tStart
%     datetime(2021,12,31),...        % tEnd
%     ["PTLC","HHZ","CM","00"],...    % SNCL1
%     ["PPLP","HHZ","EC",""],...      % SNCL2
%     false,...                           % allFlag
%     0,...                         % dW
%     1/128,...                       % lfc
%     1/32,...                       % hfc
%     1,...                         % newFs
%     2^12,...                        % secDur
%     2^11,...                        % maxLag
%     1,...                          % Nsmooth
%     true,...                        % transferFlag
%     true,...                       % verboseFlag
%     '~/rawdata/'};                  % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2018,04,20),...        % tStart
%     datetime(2018,04,30),...        % tEnd
%     ["SN02","HHZ","9D",""],...      % SNCL1
%     ["SN06","HHZ","9D",""],...      % SNCL2
%     0,...                           % allFlag
%     0,...                           % dW
%     1/8,...                        % lfc
%     2,...                           % hfc
%     16,...                          % newFs
%     2^10,...                        % secDur
%     2^07,...                        % maxLag
%     1,...                           % Nsmooth
%     false,...                       % transferFlag
%     true,...                        % verboseFlag
%     '~/rawdata/'};                  % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2020,03,01),...        % tStart
%     datetime(2020,05,10),...        % tEnd
%     ["SAGA","HHZ","EC",""],...      % SNCL1
%     ["TAIS","HHZ","EC",""],...      % SNCL2
%     false,...                       % allFlag
%     0,...                           % dW
%     1/512,...                       % lfc
%     1/2,...                         % hfc
%     16,...                          % newFs
%     2^10,...                        % secDur
%     2^07,...                        % maxLag
%     1,...                           % Nsmooth
%     true,...                        % transferFlag
%     true,...                        % verboseFlag
%     '~/rawdata/'};                  % pathToWaveformServerOrDirectory

% functionDefaults = {...
%     datetime(2019,10,01),...        % tStart
%     datetime(2021,01,01),...        % tEnd
%     ["PTLC","HHZ","CM","00"],...    % SNCL1
%     ["PPLP","HHZ","EC",""],...      % SNCL2
%     1,...                           % allFlag
%     0.1,...                         % dW
%     1/128,...                       % lfc
%     1/064,...                       % hfc
%     1/4,...                         % newFs
%     2^12,...                        % secDur
%     2^10,...                        % maxLag
%     21,...                          % Nsmooth
%     true,...                        % transferFlag
%     false,...                       % verboseFlag
%     '~/rawdata/'};                  % pathToWaveformServerOrDirectory

% SAGA and ZUMB give nice results in 8 - 16 sec. band.
% SAGA and PUYO give nice results in 8 - 16 sec. band.

% functionDefaults = {...
%     datetime(2019,10,01),...        % tStart
%     datetime(2021,01,01),...        % tEnd
%     ["PTLC","HHZ","CM","00"],...    % SNCL1
%     ["PPLP","HHZ","EC",""],...      % SNCL2
%     1,...                           % allFlag
%     0.1,...                         % dW
%     1/128,...                       % lfc
%     1/064,...                       % hfc
%     0.5,...                         % newFs
%     2^12,...                        % secDur
%     2^10,...                        % maxLag
%     21,...                          % Nsmooth
%     true,...                        % transferFlag
%     false,...                       % verboseFlag
%     '~/rawdata/'};                  % pathToWaveformServerOrDirectory


optsToUse = functionDefaults;
if nVarargin
    optsToUse(1:nVarargin) = varargin;
end

[tStart,tEnd,SNCL1,SNCL2,allFlag,dW,lfc,hfc,newFs,secDur,maxLag,Nsmooth,...
    transferFlag,verboseFlag,pathToWaveformServerOrDirectory] = deal(optsToUse{:});
Dorig = [];

%%
if isfloat(tEnd)
    days = tStart:tStart+tEnd-1;
else
    days = tStart:tEnd;
end
lDays = length(days);

%
stnm1 = char(SNCL1(1));
chan1 = SNCL1(2);
net1 = SNCL1(3);
locid1 = SNCL1(4);

stnm2 = char(SNCL2(1));
chan2 = SNCL2(2);
net2 = SNCL2(3);
locid2 = SNCL2(4);

charchan1 = char(chan1);
charchan1 = charchan1(1:2);

%%
az12 = 0;
rotateFlag = false;                 % default true
sameSensorFlag = sum(SNCL1 == SNCL2) == 4;
if sameSensorFlag
    %autocorrelation?
    if allFlag
        if strcmp(charchan1,"BD")
            transferFlag = false;   % force false
        end
    end
    interstation_distance = 0;
else % two different sensors
    %cross-correlation?
    refEllipse = referenceEllipsoid("wgs84");
    [stlat1,stlon1] = metaDataFromStationList(SNCL1(1));
    [stlat2,stlon2] = metaDataFromStationList(SNCL2(1));
    [interstation_distance,az12] = distance(stlat1,stlon1,stlat2,stlon2,refEllipse);
    if allFlag
        if verboseFlag
            fprintf('sensors are different, all channels requested, must rotate?\n');
        end
        transferFlag = true;        % force true
        rotateFlag = true;          % force true
    end
end

%%
rawDataDir = pathToWaveformServerOrDirectory;
rawDataDir1 = rawDataDir;
rawDataDir2 = rawDataDir;

%%
cornersfin = isfinite([lfc hfc]);
if ~any(cornersfin)
    fprintf(2,'no valid corners input, doing nothing\n')
    return;
end

%%
totN = secDur*newFs;
nfft = totN;
mxl = nfft/2;
lags = (-mxl:mxl-1)';
lags = lags/newFs;
mI = abs(lags) <= maxLag;
k = find(mI,1,'last');
mI(k) = false;
sum_mI = sum(mI);
lags = lags(mI);
nOverlap = 0; %.5;
stride = (1-nOverlap)*secDur;
maxWindows = floor((86400-secDur)/stride)+1; %floor(86400/secDur);

%% filter info
npoles = 4;
detrendFlag = true;
zeroPhaseFlag = true;
Hbu = freqOperator(nfft,lfc,hfc,newFs,npoles); % column vector

%% notes to the user
fprintf('Filter desired: <strong>%5.3f - %5.3f</strong> (%05.3f sec.,%05.3f sec.)\n',lfc,hfc,1/hfc,1/lfc);
fprintf('Total number of days to process: <strong>%d</strong>\n',lDays);
fprintf('Time range: <strong>%s</strong> - <strong>%s</strong>\n',tStart,tEnd);
fprintf('Total number of points in truncated lag time vector (+/-): %d\n',sum_mI);
fprintf('Maximum number of windows per day: <strong>%d</strong>\n',maxWindows);

%% preallocation
if allFlag
    if sameSensorFlag
        % 12: 3 for ZZ,NN,EE, 6 for caus/acaus each of: ZE, ZN, EN, + 3 for fold/sym stacks of cross terms
        nIndependentChannels = 12;
        checkLS = 3;
    else
        % 27: 9 for 3*(ZZ,RR,TT), 9 for 3*(ZR,ZT,TZ), and 9 for 3*(ZT,RZ,RT) [3 because caus/acaus/sym]
        nIndependentChannels = 27;
        checkLS = 6;
    end
elseif sameSensorFlag
    % autocorrelation _without_ allFlag (1 stream read, so no cross terms)
    % only 1 channel since caus==acaus
    nIndependentChannels = 1;
    checkLS = 1;
else
    % 2 different sensors, but no all flag (2 streams read)
    % 3 for caus/acaus/symmetric
    nIndependentChannels = 3;
    checkLS = 2;
end

%
causalLength = maxLag*newFs;
causalChannelLength = nIndependentChannels*causalLength;
numberElements = causalChannelLength*lDays*maxWindows;
totSubWindows = causalChannelLength*maxWindows;
totalIndividualTimeSlices = nIndependentChannels*maxWindows*lDays;

fprintf('\n');
if numberElements >= 1e9 && numberElements < 1e12
    neStr = sprintf('%5.2fG',numberElements/1e9);
elseif numberElements >= 1e6
    neStr = sprintf('%5.2fM',numberElements/1e6);
elseif numberElements >= 1e3
    neStr = sprintf('%5.2fK',numberElements/1e3);
else
    neStr = sprintf('%f',numberElements);
end

fprintf('Total number of elements: <strong>%s</strong>\n',neStr);
fprintf('number of "independent" channels: <strong>%d</strong>\n',nIndependentChannels);
fprintf('\n');

dayData2 = NaN(totSubWindows,lDays);
if checkLS == 2
    Dorig = NaN(sum_mI*maxWindows,lDays);
end

si = (1:nIndependentChannels*maxWindows:totalIndividualTimeSlices)';
ei = flipud((totalIndividualTimeSlices:-nIndependentChannels*maxWindows:nIndependentChannels*maxWindows)');
t = NaT(maxWindows,lDays);

%% figure out what channels we need to eventually read...
chanList1 = chan1;
locList1 = locid1;
chanList2 = chan2;
locList2 = locid2;

if ~sameSensorFlag  % different sensors
    locList2 = locid2;
    if allFlag
        %rotateFlag true
        chanList1 = [string([charchan1,'Z']); string([charchan1,'N']); string([charchan1,'E'])];
        charchan2 = char(chan2);
        charchan2 = charchan2(1:2);
        chanList2 = [string([charchan2,'Z']); string([charchan2,'N']); string([charchan2,'E'])];
    else
        chanList2 = chan2;
    end
else                % same sensor
    if allFlag
        if strcmp(charchan1,"BD")
            locList1 = ["01"; "02"; "03"];
        else
            %chanList1 = [string([charchan1,'Z']); string([charchan1,'1']); string([charchan1,'2'])];
            chanList1 = [string([charchan1,'Z']); string([charchan1,'N']); string([charchan1,'E'])];
        end
    end
end

%% get response info independent of whether deconvolution asked for (very little overhead)
if transferFlag
    sncls = [SNCL1; SNCL2];
    responseStructure = singleSNCLFreqResponse(sncls,tStart,tEnd,nfft,newFs,'disp');
    clear sncls;

    %% a bit redundant but this is how matlab wants broadcast variables to be handled within parfor loops
    Tend1 = responseStructure(1).Tend;
    ltend1 = responseStructure(1).N;
    H1_ = responseStructure(1).H;

    Tend2 = responseStructure(end).Tend;
    ltend2 = responseStructure(end).N;
    H2_ = responseStructure(end).H;
else
    H1_ = [];
    H2_ = H1_;
    Tend1 = NaT;
    Tend2 = Tend1;
    ltend1 = 0;
    ltend2 = ltend1;
end


tdwFlag = true; %good for phase resolution
deconFlag = ~true;
if deconFlag
    tdwFlag = false;
end

parfor i = 1:lDays
    dayStart = days(i);
    fprintf('%s\n',dayStart);

    % load data (1 or 3 sncls for this first step)
    H1 = Hbu;
    H2 = Hbu;

    S = loadWaveforms(dayStart,1,string(stnm1),chanList1,net1,locList1,false,verboseFlag,rawDataDir1);
    if ~sameSensorFlag %different sensors, read the next sensor...
        S = [S; ...
            loadWaveforms(dayStart,1,string(stnm2),chanList2,net2,locList2,false,...
            verboseFlag,rawDataDir2)];
    end
    % at this point i possibly have 1, 2, 3, or 6 sncls read!

    % i have to sync to make cross term sections as efficient as possible...
    % i risk possibly throwing data away when doing the autocorr section...
    % however, in theory, a single sensor should not have incongruous gaps between sncls anyways, so i should be OK

    S = detrendWaveforms(S);
    S = syncWaveforms(S,false,true,true);
    S = padWaveforms(S);
    S = resampleWaveforms(S,newFs);         % doesnt remove mean, keeps NaNs
    S = nanGapWaveforms(S,0);

    %%
    sTime = pull(S,'ref');
    goodTimes = ~isnat(sTime);
    if ~all(goodTimes)
        if verboseFlag
            fprintf("data do not sufficiently overlap, skipping\n");
        end
        continue;
    end

    %%
    lS = length(S);
    if lS ~= checkLS
        fprintf(1,"mismatch between expected (%d) and actual (%d) number of streams to read.\n",checkLS,lS);
        fprintf("double check whats going on on day: %s (if you so desire...)\n",dayStart);
        fprintf('\n');
        continue;
    end

    %%
    npts = pull(S,"npts");
    minnpts = min(npts);
    if ~all(npts == minnpts)
        if verboseFlag
            fprintf("not all sncls the same length. fixing.\n");
        end

        fI = find(npts ~= minnpts);
        for nn = 1:length(fI)
            d = S(fI(nn)).d;
            S(fI(nn)).npts = minnpts;
            S(fI(nn)).e = seconds((minnpts-1)/newFs);
            S(fI(nn)).d = d(1:minnpts);
        end
    end

    %%
    if rotateFlag
        % since most of these sensors are equatorial, im just going to
        % assume azimuth is linearly related to backazimuth
        if verboseFlag
            fprintf("rotating waveforms from 1 to 2\n");
        end
        % S = rotateWaveforms(S,theta,eci,nci,renameCmp)
        S(1:3) = rotateWaveforms(S(1:3),az12,3,2,true); % ZRT
        S(4:6) = rotateWaveforms(S(4:6),az12,3,2,true); % ZRT
    end

    %%
    if transferFlag
        H1__ = H1_;
        H2__ = H2_;
        if ltend1 > 1
            ti = find(dayStart <= Tend1,1,'first');
            H1 = H1.*H1__(:,ti);
        else
            H1 = H1.*H1__;
        end

        %%
        if ~sameSensorFlag
            if ltend2 > 1
                ti = find(dayStart < Tend2,1,'first');
                H2 = H2.*H2__(:,ti);
            else
                H2 = H2.*H2__;
            end
        end
    end

    if zeroPhaseFlag
        H1 = H1.*conj(H1);
        H2 = H2.*conj(H2);
    end

    %%
    d = pull(S);

    % here we cut data from first sensor or channel first...
    [dcut,~,endIndex,badFlag,nwindows] = cutWindows(d(:,1),totN,nOverlap,detrendFlag);
    if badFlag
        if verboseFlag
            fprintf("requested cut window is too long.\n");
        end
        continue;
    end

    %% preallocate
    tToday = getTimeVec(S(1));
    D = complex(zeros(nfft,maxWindows*lS));
    tdown = tToday(endIndex);
    tdown2 = NaT(maxWindows,1);
    tdown2(1:nwindows) = tdown;

    %%
    if verboseFlag
        fprintf('number of windows for this day: %d\n',nwindows);
    end
    sis = maxWindows*(0:lS-1)+1; %start of each block (for each sensor)

    %% here we apply all desired filters
    Dtmp = fft(dcut,nfft,1); %just the first trace of D...
    Dtmp = H1.*Dtmp;
    D(:,sis(1):sis(1)+nwindows-1) = Dtmp;

    % new preallocation of correlated streams (not necessary when lS == 1)
    corrStreams = NaN(causalLength,nIndependentChannels*maxWindows);

    % cut the rest of the channels and populate the rest of D...
    if lS > 1
        for nn = 2:lS
            si_ = sis(nn);
            dcut = cutWindows(d(:,nn),totN,nOverlap,detrendFlag);
            Dtmp = fft(dcut,nfft,1);

            if nn == 2 || nn == 3
                if sameSensorFlag
                    Dtmp = H1.*Dtmp;
                    D(:,si_:si_+nwindows-1) = Dtmp;
                else
                    Dtmp = H2.*Dtmp;
                    D(:,si_:si_+nwindows-1) = Dtmp;
                end
            else %nn==4,5,6
                Dtmp = H2.*Dtmp;
                D(:,si_:si_+nwindows-1) = Dtmp;
            end
        end
    end

    if deconFlag && lS == 2
        %just whiten the first station
        %leave the second station alone
        D1 = D(:,1:maxWindows); %<-- first station
        D1 = ifft(D1,[],1,"symmetric");
        D1 = sign(detrend(D1)); %time domain
        D1 = fft(D1,nfft,1); %back to frequency domain
        D1 = fdWhiten(D1,lfc,hfc,dW,newFs,false,verboseFlag); %whiten 1st station n freq. dom.
        D(:,1:maxWindows) = conj(D1); %take conj. transpose of 1st station
    elseif lS == 2 || lS == 6
        if tdwFlag
            if verboseFlag
                fprintf("applying one-bit time-domain normalization\n");
            end
            D = ifft(D,[],1,'symmetric');
            D = sign(detrend(D)); %time domain
            D = fft(D,nfft,1); %back to frequency domain
        end
    end

    %% cross correlation (with possible spectral whitening)
    % NOTE: at this point, D is cut, and is in freq dom and has all channels (assuming more than one channel was requested...)
    switch lS
        case 1 %(single station autocorr)
            % autocorr (no spectral whitening)
            %D = conj(D).*conj(D);
            D = D.*conj(D);
            D = ifft(D,[],1,"symmetric");
            %D = sign(D);
            %D = tdNorm(D,2/lfc,1,newFs);

            % distribute to master matrix
            D = D(1:causalLength,:); % special case where D is auto-correlated in-situ
            dayData2(:,i) = D(:);
        case 2 %(no autocorr)
            if ~deconFlag
                D = fdWhiten(D,lfc,hfc,dW,newFs,false,verboseFlag); %whiten both stations
                D2 = conj(D(:,maxWindows+1:end));     % second sensor (conjugate of copy of second half of D matrix)
                D = D(:,1:maxWindows);                % first sensor (clobbers D)
            else
                D2 = D(:,maxWindows+1:end); % second sensor
                D = D(:,1:maxWindows);      % first sensor
            end

            % first iteration
            D3 = D.*D2; % this is a cross-correlation
            fprintf("sizeD3: %d, nfft: %d\n",size(D3,1),nfft);
            D3 = ifft(D3,[],1,'symmetric'); % convert back to time domain
            D3 = fftshift(D3,1);
            D3 = D3(2:end,:);
            %D4 = D3(mI,:);
            %Dorig(:,i) = D4(:);

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3);
            corrStreams(:,1:maxWindows) = caus_(1:causalLength,:);
            corrStreams(:,maxWindows+1:2*maxWindows) = acaus_(1:causalLength,:);
            corrStreams(:,2*maxWindows+1:3*maxWindows) = sym_(1:causalLength,:);

            % distribute to master matrix
            dayData2(:,i) = corrStreams(:);
        case 3 %(single station auto- and cross-corr)
            autocorr = D.*conj(D);
            medD = median(autocorr,"all","omitnan"); %<-- scalar
            D = D./medD;
            autocorr = autocorr./medD;
            Dwhitened = D; % make copy of D, will need later
            % do autcorrs
            D = ifft(autocorr,[],1,"symmetric"); %time-domain
            corrStreams(:,1:3*maxWindows) = normalizeWaveforms(D(1:causalLength,:));        %ZZ,NN,EE (autocorrs)

            % make copy of D and whiten for cross-component correlations (but original D stays unwhitened...)
            if tdwFlag & ~deconFlag
                Dwhitened = ifft(Dwhitened,[],1,"symmetric"); %convert all to time-domain
                Dwhitened = sign((Dwhitened)); %time domain whitening of _copy_
                Dwhitened = fft(Dwhitened,nfft); % pass to frequency domain (whitened in timedom, but not in freqdom)
                Dwhitened = fdWhiten(Dwhitened,lfc,hfc,dW,newFs,false,verboseFlag); % returned in frequency domain
            end

            % do cross component corrs
            if deconFlag
                %Davg = medfiltSH(autocorr',11,true)';              %average across time
                Davg = medfiltSH(autocorr,ceil(nfft/32),true);      %average across frequency
                DavgW = Davg;
                DavgW(:,1:maxWindows) = DavgW(:,maxWindows+1:2*maxWindows); %N^2, N^2, E^2
                Davg(:,maxWindows+1:2*maxWindows) = Davg(:,2*maxWindows+1:end); % move third (E) to second (N)
                Davg(:,2*maxWindows+1:end) = Davg(:,1:maxWindows); %move first (Z) to third (E), ends up at: Z^2, E^2, Z^2
                D = Dwhitened; %return D to its original state....
                D = circshift(D,-maxWindows,2);
                Dunfilt = D;
                D(:,maxWindows+1:end) = conj(D(:,maxWindows+1:end)); %E and Z are conj...
                Dw = Hbu.*Dwhitened./DavgW;
                Dwhitened(:,1:maxWindows) = conj(Dwhitened(:,1:maxWindows)); %Z is conj...
                D = Hbu.*D./Davg;
                
            else
                D = conj(Dwhitened); %<--depending on options, could be td- and fd-whitened, or just one of the two, or none
                D = circshift(D,-maxWindows,2);
            end

            %D = Dwhitened.*D;
            %D = Dwhitened.*D - conj(Dwhitened.*conj(D));
            D = Dwhitened.*D; %<--cross correlations (deconvolution)
            if deconFlag
                D = D - Dw.*Dunfilt; %<-- cross convolutions (normalized)
            else
                D(:,1:maxWindows) = conj(D(:,1:maxWindows)); %Z is conj...
            end

            D = ifft(D,[],1,"symmetric"); % now in time domain
            D = fftshift(D,1);
            D = normalizeWaveforms(D,true);
            D = D(2:end,:);

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D(:,1:maxWindows));                  %N.Z*
            corrStreams(:,3*maxWindows+1:6*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D(:,maxWindows+1:2*maxWindows));     %N.E*
            corrStreams(:,6*maxWindows+1:9*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D(:,2*maxWindows+1:3*maxWindows));   %E.Z*
            corrStreams(:,9*maxWindows+1:12*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            % distribute to master matrix
            dayData2(:,i) = corrStreams(:);
        case 6 %(no autocorr)
            D = fdWhiten(D,lfc,hfc,dW,newFs,false,true); %verboseFlag);
            D2 = conj(D(:,3*maxWindows+1:end));     % conjugate form of 3 components of second sensor
            D = D(:,1:3*maxWindows);                % 3 components of first sensor (clobber D again)

            % first iteration, Z1Z2,R1R2,T1T2
            D3 = D.*D2;
            D3 = ifft(D3,[],1,'symmetric'); % now in time domain
            D3 = fftshift(D3,1);
            D3 = D3(2:end,:);

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,1:maxWindows));                  %Z1Z2 (1:3)
            corrStreams(:,1:3*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,maxWindows+1:2*maxWindows));     %R1R2 (4:6)
            corrStreams(:,3*maxWindows+1:6*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,2*maxWindows+1:3*maxWindows));   %T1T2 (7:9)
            corrStreams(:,6*maxWindows+1:9*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];


            % second iteration, Z1R2,R1T2,T1Z2 (10:18)
            D2 = circshift(D2,-maxWindows,2);
            D3 = D.*D2;
            D3 = ifft(D3,[],1,'symmetric'); % now in time domain
            D3 = fftshift(D3,1);
            D3 = D3(2:end,:);
            %fprintf('sizeD3: %d, nfft: %d\n',size(D3,1),nfft);

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,1:maxWindows));                  %Z1R2 (10:12)
            corrStreams(:,9*maxWindows+1:12*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,maxWindows+1:2*maxWindows));     %R1T2 (13:15)
            corrStreams(:,12*maxWindows+1:15*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,2*maxWindows+1:3*maxWindows));   %T1Z2 (16:18)
            corrStreams(:,15*maxWindows+1:18*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            % third iteration, Z1T2,R1Z2,T1R2 (19:27)
            D2 = circshift(D2,-maxWindows,2);
            D3 = D.*D2;
            D3 = ifft(D3,[],1,'symmetric'); % now in time domain
            D3 = fftshift(D3,1);
            D3 = D3(2:end,:);

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,1:maxWindows));                  %Z1T2
            corrStreams(:,18*maxWindows+1:21*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,maxWindows+1:2*maxWindows));     %R1Z2
            corrStreams(:,21*maxWindows+1:24*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            [caus_,acaus_,sym_] = decomposeNoiseCorrMatrix(D3(:,2*maxWindows+1:3*maxWindows));   %T1R2
            corrStreams(:,24*maxWindows+1:27*maxWindows) = ...
                [caus_(1:causalLength,:) acaus_(1:causalLength,:) sym_(1:causalLength,:)];

            % distribute to master matrix
            dayData2(:,i) = corrStreams(:);
    end
    t(:,i) = tdown2;
end

%%
fprintf('here, i am populating the dayData matrix (via a reshape command)\n');
dayData = NaN(causalLength,totalIndividualTimeSlices);
for i = 1:lDays
    dayData_ = dayData2(:,i);
    dayData_ = reshape(dayData_,causalLength,nIndependentChannels*maxWindows);
    dayData(:,si(i):ei(i)) = dayData_;
end
t = t(:);
clear dayData_ H1__ H2__ responseStructure %free up some memory!!

%%
if exist('interstation_distance','var')
    interstation_distance = interstation_distance*1e-3;
    fprintf('interstation distance: %f km.\n',interstation_distance);
end

% % % continueFlag = false;
% % % if ~continueFlag
% % %     return;
% % % end

%%
if exist('Nsmooth','var')
    if isempty(Nsmooth)
        Nsmooth = max([4 2^(nextpow2(maxWindows)-2)]);
    end
else
    Nsmooth = max([4 2^(nextpow2(maxWindows)-2)]);
end

plotFlag = true;
%plotFlag = true; stackMethod = 'pws';
%Nsmooth = max([8 2^(nextpow2(maxWindows)-3)]); %-3 --> approximately 4.5 hour-long windows, -2 --> 9.1 hours, etc


smoothFlag = true;
% here loop through each "independent" channels first then over number of days to apply smoothing
if Nsmooth > 1 & smoothFlag
    widthInHours = (1-nOverlap)*(secDur*Nsmooth)/3600;
    fprintf('maxWindows: %u; Nsmooth: %u; Width in Approximate Hours: <strong>~%05.2f</strong>\n',...
        maxWindows,Nsmooth,widthInHours);
    box = ones(Nsmooth,1)/Nsmooth;
    dayData2 = NaN(causalLength,maxWindows*lDays);
    for i = 1:nIndependentChannels
        if verboseFlag
            fprintf('working on smoothing channel %d, Nsmooth: %d\n',i,Nsmooth);
        end

        startBlockI = maxWindows*(i-1)+1;
        endBlockI = maxWindows*i;
        for j = 1:lDays
            si = startBlockI + (j-1)*maxWindows*nIndependentChannels;
            ei = endBlockI + (j-1)*maxWindows*nIndependentChannels;

            si2 = maxWindows*(j-1)+1;
            ei2 = maxWindows*j;
            dayData2(:,si2:ei2) = dayData(:,si:ei);
        end

        %
        nanI = ~isfinite(dayData2);
        dayData2(nanI) = 0;

        dayData2 = fftfilt(box,dayData2')';
        wienerFilter = true;
        if wienerFilter
            wienerSmoother = [ceil(2*newFs/lfc),Nsmooth];
            dayData2 = wiener2(dayData2,wienerSmoother);
            disp(wienerSmoother);
        end
        %dayData2 = fftfilt(box,dayData2')';

        %
        for j = 1:lDays
            si = startBlockI + (j-1)*maxWindows*nIndependentChannels;
            ei = endBlockI + (j-1)*maxWindows*nIndependentChannels;

            si2 = maxWindows*(j-1)+1;
            ei2 = maxWindows*j;
            dayData(:,si:ei) = dayData2(:,si2:ei2);
        end
    end
else
    fprintf("not applying smoothing...\n");
end

%%
if ~plotFlag
    refTrace = [];
    return
end

%%
close all;
stackMethod = 'pws';
refTrace = NaN(causalLength,nIndependentChannels);
dayData2 = NaN(causalLength,maxWindows*lDays);
referenceStartTime = tEnd + 1;
if checkLS == 6
    titles = ["z1z2 caus";"z1z2 acaus";"z1z2 sym";...
        "r1r2 caus";"r1r2 acaus";"r1r2 sym";...
        "t1t2 caus";"t1t2 acaus";"t1t2 sym";...
        "z1r2 caus";"z1r2 acaus";"z1r2 sym";...
        "r1t2 caus";"r1t2 acaus";"r1t2 sym";...
        "t1z2 caus";"t1z2 acaus";"t1z2 sym";...
        "z1t2 caus";"z1t2 acaus";"z1t2 sym";...
        "r1z2 caus";"r1z2 acaus";"r1z2 sym";...
        "t1r2 caus";"t1r2 acaus";"t1r2 sym"];
end

axY1 = gobjects(nIndependentChannels,1);
for i = 1:nIndependentChannels
    fprintf('processing channel %d\n',i);
    startBlockI = maxWindows*(i-1)+1;
    endBlockI = maxWindows*i;
    for j = 1:lDays
        si = startBlockI + (j-1)*maxWindows*nIndependentChannels;
        ei = endBlockI + (j-1)*maxWindows*nIndependentChannels;
        si2 = maxWindows*(j-1)+1;
        ei2 = maxWindows*j;
        %disp([j si ei si2 ei2])
        dayData_ = dayData(:,si:ei);
        zeroI = dayData_ == 0;
        dayData_(zeroI) = NaN;
        dayData2(:,si2:ei2) = dayData_;
    end

    %
    if strcmp(stackMethod,'pws')
        nu = 2;
        refTrace_ = pws(dayData2(:,t <= referenceStartTime),false,false,nu);
    elseif strcmp(stackMethod,'med')
        refTrace_ = median(dayData2(:,t <= referenceStartTime),2,'omitnan');
    elseif strcmp(stackMethod,'mean')
        refTrace_ = mean(dayData2(:,t <= referenceStartTime),2,'omitnan');
    end

    % if ismember(i,[1 2]) %rayleigh wave components
    %     %ismember(i,[2 5 11 23]) %rayleigh wave components
    %     %if ismember(i,[1 2 4 5 10 11 22 23]) %rayleigh wave components
    %     dayData2 = dayData2./rssq(dayData2);
    %     dayData3 = cat(3,dayData3,dayData2);
    % end
    refTrace(:,i) = refTrace_;

    figure('units','normalized','outerposition',[0 0 2/5 1]);
    ax__(1) = subplot(5,6,[1 2 3 4 5 6]);
    mp = find(lags >= 0,1);
    causalLength1 = causalLength - 1;
    plot(lags(mp:mp+causalLength1-1),refTrace_(1:causalLength1),'linewidth',3);
    cc = colorbar;
    cc.Visible = 'off';
    grid on;
    zoom on;
    axis tight;

    if checkLS == 6
        title(titles(i));
    end

    ax__(2) = subplot(5,6,sort([7 8 13 14 19 20 25 26 9 10 11 12 15 16 17 18 21 22 23 24 27 28 29 30]));
    imagesc(ax__(2),lags(mp:mp+causalLength1-1),(1:length(t))',sign(dayData2'));
    ylim(ax__(2),[0 size(dayData2,2)]);

    colorbar;
    zoom on;
    linkaxes(ax__,'x');
    axY1(i) = ax__(1);
end
linkaxes(axY1,'y');

for i = 1:nIndependentChannels
    figure(i);
    ylim([0 size(dayData2,2)]);
end