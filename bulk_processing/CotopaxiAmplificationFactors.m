clear;
close all;
% 2-4 didnt really pan out, not a lot of coherent phases.
% but i still kept the average ratios
lfc = 1/2;
hfc = 1;
Fs = hfc*8;

dayStart = datetime(2021,01,01);
dayEnd = datetime(2021,01,31);
dayVec = (dayStart:dayEnd)';
lDays = length(dayVec);

vecLength = 32767; %65534; %491505;
%ccmean = NaN(491505,lDays);
%ccmean = NaN(245745,lDays);
ccmean = NaN(vecLength,lDays);
ccpws = ccmean;

%kstnms = ["BREF";"BTAM";"BVC2";"BMOR";"VCES";"BNAS"];
kstnms = ["BREF";"BTER"]; %"BTAM";"BVC2";"BMOR";"VCES";"BNAS"];

nn = 0;
%kk = 1;
for kk = 1:lDays
    tic;
    day_ = dayVec(kk);
    disp(day_);
    S = loadWaveforms(day_,1,kstnms,["BHZ";"HHZ"],"EC");
    refs = pull(S,'ref');
    goodI = ~isnat(refs);
    if sum(goodI) < length(kstnms)
        continue;
    end

    Sorig = S;
    S = resampleWaveforms(detrendWaveforms(S),Fs);
    S = syncWaveforms(S,false,true,false);
    %S = filterWaveforms(S,2,4);
    S = scaleWaveforms(transferWaveforms(S,lfc,hfc,4,Fs,"vel",1,false),1e9);
    S = nanGapWaveforms(S,0);
    S = detrendWaveforms(S);
    totDur = seconds(min(pull(S,'e')));
    if totDur < 3600
        continue;
    end

    %%
    secDur = 3600/2;
    winlen = secDur*Fs;
    nOverlap = 0;
    refs = pull(S,'ref');
    badI = isnat(refs);
    S(badI) = [];
    lS = length(S);
    toc;

    %
    nWindows = 50; %1440/2;
    t = getTimeVec(S);
    lT = length(t) - winlen + 1;
    nI = sort(randsample(lT,nWindows));

    %
    dcut = [];
    denv = [];

    for i = 1:lS
        d_ = double(pull(S(i)));
        denv_ = abs(hilbert(d_));
        denv_ = zpkFilter(denv_,-inf,1/80,1,1,1);
        dcut_ = NaN(winlen,nWindows); %d_(nI); %cutWindows(d_,winlen,nOverlap,true);
        denv2 = dcut_;
        for j= 1:nWindows
            denv__ = denv_(nI(j):nI(j)+winlen-1);
            d__ = detrend(d_(nI(j):nI(j)+winlen-1));
            d__ = d__./denv__;      % time-domain normalization
            d__ = fdWhiten(d__,lfc,hfc,0,Fs,true,false);
            d__ = detrend(d__);
            dcut_(:,j) = d__./rssq(d__); %rssq(denv__); %normalize
            denv2(:,j) = denv__;
        end
        dcut = [dcut; dcut_];
        denv = [denv; denv2];
    end
    toc;

    %%
    dcut = dcut./rssq(dcut);
    %cc = NaN(245745,nWindows);
    cc = NaN(vecLength,nWindows);
    %cc = NaN(491505,nWindows);
    for i = 1:nWindows
        %disp(i);
        dcut_ = dcut(:,i);
        dcut_ = reshape(dcut_,[winlen,lS]);
        %dcutNum = [repmat(dcut_(:,1),1,lS-1) repmat(dcut_(:,2),1,lS-2) repmat(dcut_(:,3),1,lS-3) ...
        %    repmat(dcut_(:,4),1,lS-4) dcut_(:,5)];
        %dcutDenom = [dcut_(:,2:end) dcut_(:,3:end) dcut_(:,4:end) dcut_(:,5:end) dcut_(:,6)];
        cc_ = doCrossCorrFreqDom(dcut_(:,1),dcut_(:,2));
        %cc_ = doCrossCorrFreqDom(dcutDenom,dcutNum);
        cc(:,i) = cc_(:);
    end
    nn = nn+1;
    cc(~isfinite(cc)) = 0;
    cc__ = mean(cc,2,"omitnan");
    ccmean(:,nn) = cc__;
    cc__ = pws(cc);
    ccpws(:,nn) = cc__;
    toc;

    disp(" ");

end
