clear; %close all; clc;

Fs = 20;
secDur = 60;
%tStart = datetime(2022,01,01);
tStart = datetime(2006,01,180);
tEnd = datetime(2022,05,26);
dayVec = (tStart:tEnd)';
ldays = length(dayVec);

winlen = secDur.*Fs;
kstnm = ["BRUN";"BULB";"BPAT";"BMAS";"BBIL";"BREF";"BNAS";"BTAM";"BMOR";"BVC2"];
kcmpnms = "BHZ";
lkstnms = length(kstnm);

%
nOverlap = 0;
maxWindows = 86400/secDur;
tMain = NaT(maxWindows,ldays);
ampMain = NaN(maxWindows,lkstnms,ldays);
winPerDay = NaN(ldays,1);

%%
lfc = 2;
hfc = 3;
detrendFlag = true;
R = singleSNCLFreqResponse(["BULB","BHZ","EC",""],datetime(2010,01,01),datetime(2010,01,02),0,Fs,'disp');
zeroes = R.zeros;
poles = R.poles;
constant = R.constant;
npoles = 4;
direction = true;

parfor i = 1:ldays
    day_ = dayVec(i);
    S = loadWaveforms(day_,1,kstnm,kcmpnms,"EC","");

    refs = pull(S,'ref');
    badI = isnat(refs);
    S(badI) = [];
    lS = length(S);

    if ~lS
        fprintf('no data on day: %s\n',datestr(day_));
        continue;
    end

    S = differentiateWaveforms(S);
    S = resampleWaveforms(S,Fs);
    S = interpolateWaveforms(S);
    S = transfer(S,zeroes,poles,constant,lfc,hfc,npoles,direction);
    S = scaleWaveforms(S,1e9);  % convert to nanometers
    S = nanGapWaveforms(S,0);   % in theory removes non-finites
    S = padWaveforms(S);
    S = syncWaveforms(S);

    minB = min(pull(S,'ref'));
    maxE = max(pull(S,'ref')+pull(S,'e'));

    tStart = dateshift(minB,'end','minute');
    tEnd = dateshift(maxE,'start','minute');

    dur_ = seconds(tEnd-tStart);
    if dur_ < secDur
        fprintf('snippet too short, skipping: %s\n',datestr(day_));
        continue;
    end
    S = cutWaveforms(S,tStart,0,dur_);
    refs = pull(S,'ref');
    badI = isnat(refs);
    S(badI) = [];
    lS = length(S);

    if ~lS
        fprintf('no data on day: %s\n',datestr(day_));
        continue;
    end

    dcut = [];
    t = getTimeVec(S);
    amp_ = NaN(maxWindows,lkstnms);
    t_ = NaT(maxWindows,1);

    %%
    lia = ismember(kstnm,pull(S,'kstnm'));
    locbs = find(lia)';
    n = 0;
    for j = locbs
        n = n+1;
        d_ = double(pull(S(n)));
        [dcut_,~,ei,~,nw_] = cutWindows(d_,winlen,nOverlap,detrendFlag);
        dcut_ = rms(dcut_)';
        amp_(1:nw_,j) = dcut_;
        t_(1:nw_) = t(ei);
    end

    %%
    tMain(:,i) = t_;
    ampMain(:,:,i) = amp_;
    fprintf('done with: <strong>%s</strong>\n',datestr(day_));
end

%%
cd ~/research/now/hvsr/
save('JICA_correctedRSAM_2Hz3Hz','ampMain','tMain');
