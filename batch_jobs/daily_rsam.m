% processed up to 2011.188 on macmini. need 189 onwards
clear; close all; clc;
dayStart = datetime(2009,01,339);
dayEnd = datetime(2023,01,320); %dateshift(datetime('now')+hours(5),'start','day');
dayInc = 1;
%dayVec = (dayStart:dayInc:dayEnd)';
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

homeDir = fullfile('~','products','all_rsam');
cd(homeDir);

chanList = ["HHZ";"SHZ";"BHZ";"BLZ";"ENZ";"EHZ";"HNZ"];
lfcs = [1/50; 0.2; 0.25; 0.3; 0.5; 0.6;  1; 2;  5;   10];
hfcs = [1/20; 0.8;    2; 0.5;   4; 1.2; 16; 5; 15; -inf];
minLfc = min(lfcs);
rsamDur = 60;

%%
meanFlags = [true; true; false; false];
rmsFlags = [true; false; true; false];
llfc = length(lfcs);
lFlags = length(meanFlags);
nAmpCombos = llfc*lFlags;
npoles = 4;
HORFLAG = ~true;
resampleFs = 50;
MAX_RSAM_OBS = ceil(86400/rsamDur);
EF = 4;

%%
for i = 1:lDays
    day_ = dayVec(i);
    yyyy = year(day_);
    doy = day(day_,'doy');

    timerStart1 = tic;
    S = seedData2(day_,NaT,NaT,HORFLAG,true,chanList,["CM";"EC"]);
    telapsed = toc(timerStart1);
    disp(telapsed)
    lS = length(S);
    if ~lS
        fprintf('no data, skipping\n');
        continue;
    end

    refs = pull(S,'ref');
    badI = isnat(refs);
    S(badI) = [];
    lSNCLs = length(S);
    if ~lSNCLs
        fprintf('too many files, possibly multiplexed, skipping\n');
        continue;
    end

    S = differentiateWaveforms(S); %toc; fprintf('finish diff\n');
    S = nanGapWaveforms(S,0); %toc; fprintf('finish fill gaps\n');
    S = resampleWaveforms(S,resampleFs); %toc; fprintf('finish resample\n');
    S = taperWaveforms(S,resampleFs/minLfc); %toc; fprintf('finish taper\n');
    S = intWaveforms(S); %toc; fprintf('finish int\n');
    S = padWaveforms(S); %toc; fprintf('finish pad\n');

    knetwks = pull(S,'knetwk');
    kstnms = pull(S,'kstnm');
    kholes = pull(S,'khole');
    kcmpnms = pull(S,'kcmpnm');
    refs = pull(S,'ref');
    dOrig = pull(S);

    %keep metadata in S
    for j = 1:lSNCLs
        S(j).d = [];
    end
    fprintf('finish emptying S\n');

    npts = size(dOrig,1);
    dOrig = fft(dOrig);
    fprintf('finished FFT\n');

    Hbu = NaN(npts,llfc);
    for k = 1:llfc
        lfc = lfcs(k);
        hfc = hfcs(k);
        Hbu_ = freqOperator(npts,lfc,hfc,resampleFs,npoles);
        Hbu(:,k) = Hbu_;
    end
    fprintf('finished all filter coeffs\n');

    for j = 1:lSNCLs
        kstnm_ = kstnms(j);
        kcmpnm_ = kcmpnms(j);
        knetwk_ = knetwks(j);
        khole_ = kholes(j);
        ref = refs(j);
        dCurrent = dOrig(:,j);
        Scopy = S(j);
        mainAmpMatrix = zeros(MAX_RSAM_OBS,nAmpCombos);
        n = 0;
        tic;
        for k = 1:llfc
            Hbu_ = Hbu(:,k); %freqOperator(npts,lfc,hfc,resampleFs,npoles);
            df = dCurrent.*Hbu_;
            df = ifft(df,'symmetric');

            for l = 1:lFlags
                n = n+1;
                meanFlag_ = meanFlags(l);
                rmsFlag_ = rmsFlags(l);

                [ampVec,winlen] = amplitudeVector(df,resampleFs,rsamDur,meanFlag_,rmsFlag_);
                newFs = 1/rsamDur;
                newRef = dateshift(ref,'end','minute');
                iStart = t2i(newRef,ref,rsamDur);
                ampVec = ampVec(iStart:winlen:end);
                lAmpVec = length(ampVec);
                mainAmpMatrix(1:lAmpVec,n) = ampVec;
            end
        end

        %%
        mainAmpMatrix = mainAmpMatrix(:);
        S2 = dealHeader(Scopy,mainAmpMatrix,1,newRef);

        %%
        newDir = fullfile(homeDir,num2str(yyyy),knetwk_,kstnm_,sprintf('%s.D',kcmpnm_));
        if ~exist(newDir,'dir')
            mkdir(newDir);
        end
        cd(newDir);

        %%
        fname = sprintf("%s.%s.%s.%s",knetwk_,kstnm_,khole_,kcmpnm_);
        fname = char(fname);
        mkmseed(fname,S2.d,datenum(S2.ref),1./S2.delta,EF); toc;
    end
    telapsed2 = toc(timerStart1);
    fprintf('elapsed time for %s: %f\n',day_,telapsed2);
    fprintf('\n');
end
