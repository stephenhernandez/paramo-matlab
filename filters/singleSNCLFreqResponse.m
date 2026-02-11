function [responseStructure,f] = singleSNCLFreqResponse(sncls,t1,t2,nfft,newFs,units,lfc,hfc,npoles)
%
% i have to confirm, but i believe this code assumes that the response is
% the same amongst all channels for any given sensor. there may be subtle
% differences in the response values, but for the purposes of my work, we
% can assume they're the same. In this logic, I also assume that the
% station location doesnt change either.
%

%
% clear; close all;
% sncls = ["PTLC","HHZ","CM","00"];
% sncls = [sncls; ["BTER","HHZ","EC",""]];
% sncls = [sncls; ["SAGA","HHZ","EC",""]];
% sncls = [sncls; ["BREF","BHZ","EC",""]];
%
% t1 = datetime(2019,01,01); t2 = dn2dt(floor(now)); %,2^12,64);
% nfft = 2^12;
% newFs = 64;
%

if nargin < 6
    units = 'vel';
end

if nargin < 7
    lfc = -inf;
end

if nargin < 8
    hfc = -inf;
end

if nargin < 9
    npoles = 4;
end

%%
try
    load('~/igdata/ecuadorSensorDataTable10','ZPT','tStart','constant',...
        'stla','stlo','stel','allSNCLs','tEnd','sensitivity','a0');
catch
    try
        load('~/Desktop/ecuadorSensorDataTable10','ZPT','tStart','constant',...
            'stla','stlo','stel','allSNCLs','tEnd','sensitivity','a0');
    catch
        return;
    end
end

kcmpnms = sncls(:,2);
kholes = sncls(:,4);
kstnms = sncls(:,1);
knetwks = sncls(:,3);

%%
SNCLs = strcat(knetwks,kstnms,kholes,kcmpnms); %this sorts my vector... could be bad
lSNCLs = length(SNCLs);

%%
if isscalar(nfft)
    nfft = repmat(nfft,lSNCLs,1);
else
    if length(nfft) ~= lSNCLs
        fprintf('nfft and sncls are not the same size\n');
        return;
    end
end

%%
if isscalar(newFs)
    newFs = repmat(newFs,lSNCLs,1);
else
    if length(newFs) ~= lSNCLs
        fprintf('newFs and sncls are not the same size\n');
        return;
    end
end

%%
instTypes = ZPT.instType;
gains = ZPT.instGain;
responseStructure = populateResponseStructure(lSNCLs);
f = 0;

%%
for i = 1:lSNCLs
    SEISMOMETER = true;
    responseStructure(i).knetwk = knetwks(i);
    responseStructure(i).kstnm = kstnms(i);
    responseStructure(i).khole = kholes(i);
    responseStructure(i).kcmpnm = kcmpnms(i);

    chanCode = char(kcmpnms(i));
    chanCode = string(chanCode(1,1:2));

    thisSNCL = SNCLs(i);
    lia = strcmp(allSNCLs,thisSNCL);
    if ~sum(lia)
        fprintf('no response found, returning nothing\n');
        return;
    end

    lia = find(lia);
    tStart_ = tStart(lia);
    tEnd_ = tEnd(lia);

    gains_ = gains(lia);
    instTypes_ = instTypes(lia);
    ZPT_ = ZPT(lia,:);
    constant_ = constant(lia);
    stla_ = stla(lia);
    stlo_ = stlo(lia);
    stel_ = stel(lia);
    sensitivity_ = sensitivity(lia);
    a0_ = a0(lia);

    [tStart_,sortI] = sort(tStart_);
    tEnd_ = tEnd_(sortI);
    ZPT_ = ZPT_(sortI,:);
    constant_ = constant_(sortI);
    stla_ = stla_(sortI);
    stlo_ = stlo_(sortI);
    stel_ = stel_(sortI);
    sensitivity_ = sensitivity_(sortI);
    a0_ = a0_(sortI);
    gains_ = gains_(sortI);
    instTypes_ = instTypes_(sortI);

    %%
    %disp([t1 t2])
    tIs = find(t1 >= tStart_,1,'last');
    tIe = find(t2 <= tEnd_,1,'first');

    %%
    if isempty(tIs)
        tIs = 1;
    end

    if tIe < tIs % added on 24 sep 2021 to apply to BMAS BDF
        tIe = tIs;
    end

    goodI = (tIs:tIe)';
    lgood = length(goodI);

    if ~lgood
        continue;
    end
    Tstart = tStart_(goodI);
    Tend = tEnd_(goodI);

    Tstart(1) = t1;
    Tend(end) = t2;

    Stla = stla_(1);
    Stlo = stlo_(1);
    Stel = stel_(1);

    %%
    n_ = nfft(i);
    fs_ = newFs(i);
    H = NaN(n_,lgood);
    for j = 1:lgood
        lia_ = goodI(j);
        Normalization = a0_(lia_);
        constTmp = constant_(lia_);
        Sensitivity = sensitivity_(lia_);
        Gain = gains_(lia_);
        if Normalization == 57150620 || Normalization == 571506000
            Normalization = 571506200; %Normalization*10;
            constTmp = Sensitivity*Normalization;
        end

        zeroTmp = ZPT_.zeros{lia_};
        polesTmp = ZPT_.poles{lia_};

        responseStructure(i).instType  = instTypes_(lia_);
        responseStructure(i).poles = polesTmp;
        responseStructure(i).constant = constTmp;
        responseStructure(i).gain = Gain;
        responseStructure(i).sensitivity = Sensitivity;
        responseStructure(i).Fs = fs_;
        responseStructure(i).A0 = Normalization;

        if (strcmp(chanCode,"HN") || strcmp(chanCode,"EN")) ...
                && ~contains(instTypes_(lia_),"RT")
            disp(instTypes_(lia_));
            %strcmp(instTypes_(lia_),"RT")
            SEISMOMETER = false;
            zeroTmp = [];
            if (constTmp < 1 || Sensitivity < 1) && Gain %CMG
                Sensitivity= 312500;
                if Normalization < 0
                    Normalization = -Normalization;
                end
                responseStructure(i).A0 = Normalization;
                constTmp = Gain*Sensitivity*Normalization;
                responseStructure(i).constant = constTmp;
                responseStructure(i).sensitivity = Sensitivity;
            end

            constTmp = 1/constTmp;
            if strcmp(units,'disp')         %dispFlag
                zeroTmp = [0 0 zeroTmp];    %#ok<AGROW>
            elseif strcmp(units,'vel')      %velFlag
                zeroTmp = [0 zeroTmp];      %#ok<AGROW>
            elseif strcmp(units,'acc')      %accFlag
                %zeroTmp = zeroTmp;
            end
        elseif (strcmp(chanCode,"HN") || strcmp(chanCode,"EN")) && ...
                (contains(instTypes_(lia_),"RT") || contains(instTypes_(lia_),"5TC"))
            SEISMOMETER = false;
            constTmp = 1/constTmp;
            if strcmp(units,'disp')         %dispFlag
                zeroTmp = zeroTmp(1:end);
            elseif strcmp(units,'vel')      %velFlag
                zeroTmp = zeroTmp(2:end);
            elseif strcmp(units,'acc')      %accFlag
                zeroTmp = zeroTmp(3:end);
            end
        end

        if SEISMOMETER
            if ~isfinite(constTmp) || constTmp < 1
                constTmp = 1;
            end
            constTmp = 1/constTmp;


            if strcmp(units,'disp')         %dispFlag
                zeroTmp = zeroTmp(1:end);
            elseif strcmp(units,'vel')      %velFlag
                zeroTmp = zeroTmp(2:end);
            elseif strcmp(units,'acc')      %accFlag
                zeroTmp = zeroTmp(3:end);
            end
        end

        %%
        responseStructure(i).zeros = zeroTmp;
        if n_
            [Hdecon,f] = cmplxResp(n_,polesTmp,zeroTmp,constTmp,fs_);
            if any(isfinite([lfc hfc]))
                Hbu = freqOperator(n_,lfc,hfc,fs_,npoles);
                Hdecon = Hdecon.*Hbu;
            end
            H(:,j) = Hdecon;
        end

    end

    %%
    responseStructure(i).N = lgood;
    responseStructure(i).Tstart = Tstart;
    responseStructure(i).Tend = Tend;
    responseStructure(i).Stel = Stel;
    responseStructure(i).Stla = Stla;
    responseStructure(i).Stlo = Stlo;
    responseStructure(i).H = H;
end
