function [S,modelNames] = transferWaveforms(S,varargin)
%
% transferWaveforms apply/remove response. so far, only remove works.
%
% S = transferWaveforms(S,lfc,hfc,npoles,newFs,units,direction,waFlag,CANUSEGPU,zeroPhaseFlag)
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Sunday, Jul 28, 2019
%

%%
nVarargin = length(varargin);
functionDefaults = {...
    -inf,...    % lfc
    -inf,...    % hfc
    4,...       % npoles
    false,...   % newFs
    "disp",...  % units
    true,...    % DECONFLAG (true = deconvolve)
    false,...   % WAFLAG (true => deconvolve instrument first, then convolve with WA response)
    false,...   % CANUSEGPU
    false};     % zeroPhaseFlag

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[lfc,hfc,npoles,newFs,units,DECONFLAG,WAFLAG,CANUSEGPU,zeroPhaseFlag] = deal(optsToUse{:});

%%
if isempty(lfc)
    lfc = -inf;
end

if isempty(hfc)
    hfc = -inf;
end

if isempty(npoles)
    npoles = 4;
end

if isempty(newFs)
    newFs = false;
end

if isempty(units)
    units = "disp";
end

%%
if WAFLAG
    [wapoles,wazeros,wasensitivity] = readSensorPolesZeros(10); %wood-anderson response
    %% force units to be displacement
    if ~strcmp(units,"disp")
        fprintf("forcing units to be displacement before applying WA\n")
        units = "disp";
    end
    wasensitivity = wasensitivity*1000; %convert to mm
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
[~,~,~,~,~,~,~,ic] = getPolesZeros(S);
uniqueic = unique(ic);
luic = length(uniqueic);
modelNames = repmat("",lS,1);

%%
for i = 1:luic
    thisI = i == ic;
    thisS = S(thisI);

    %%
    validTracesI = find(~isnat(pull(thisS,"ref")));
    sumGood = length(validTracesI);
    if ~sumGood
        fprintf("no valid traces, skipping...\n");
        continue;
    end

    Fs = 1./thisS(1).delta;
    fI = Fs >= 1;
    Fs(fI) = round(Fs(fI));

    t1 = thisS(1).ref;
    t2 = dateshift(t1,"end","day");

    thisSNCL = [thisS(1).kstnm,thisS(1).kcmpnm,thisS(1).knetwk,thisS(1).khole];
    if strcmp(thisS(1).kstnm,"BNAS") && t1 >= datetime(2023,01,06) && t1 <= datetime(2024,01,16)
        disp("switching BNAS to CO1V");
        thisSNCL = ["CO1V","HHZ",thisS(1).knetwk,thisS(1).khole];
        R = singleSNCLFreqResponse(thisSNCL,t1,t2,0,Fs,units);
    elseif strcmp(thisS(1).kstnm,"BTER") && t1 >= datetime(2022,01,05)
        disp("switching BTER to CO1V");
        thisSNCL = ["CO1V","HHZ",thisS(1).knetwk,thisS(1).khole];
        R = singleSNCLFreqResponse(thisSNCL,t1,t2,0,Fs,units);
    elseif strcmp(thisS(1).kstnm,"PULU") %&& t1 >= datetime(2020,12,03)
        fprintf("reducing instrument sensitivity at PULU\n");
        R = singleSNCLFreqResponse(thisSNCL,t1,t2,0,Fs,units);
        R.sensitivity = R.gain*419430; %new digitizer is already 4 times smaller than Andreas data
        R.constant = R.sensitivity * R.A0;
    elseif strcmp(thisS(1).kstnm,"COTA")
        fprintf("reducing instrument sensitivity at COTA\n");
        thisSNCL = ["COTA","HHZ",thisS(1).knetwk,thisS(1).khole];
        R = singleSNCLFreqResponse(thisSNCL,t1,t2,0,Fs,units);
        R.sensitivity = R.sensitivity/4;
        R.constant = R.sensitivity * R.A0;
    elseif strcmp(thisS(1).kstnm,"SUCR") || strcmp(thisS(1).kstnm,"SRAM") || strcmp(thisS(1).kstnm,"TOMA")
        disp("switching SUCR/SRAM/TOMA to BRRN...");
        thisSNCL = ["BRRN",thisS(1).kcmpnm,thisS(1).knetwk,thisS(1).khole];
        R = singleSNCLFreqResponse(thisSNCL,t1,t2,0,Fs,units);
    elseif strcmp(thisS(1).kcmpnm,"SHZ") || strcmp(thisS(1).kcmpnm,"SHN") || strcmp(thisS(1).kcmpnm,"SHE")
        disp("switching generic short period to VC1...");
        thisSNCL = ["VC1",thisS(1).kcmpnm,thisS(1).knetwk,thisS(1).khole];
        R = singleSNCLFreqResponse(thisSNCL,datetime(2023,01,01),datetime(2023,01,02),0,Fs,units);
        R.gain = 1;
        R.A0 = 1;
        R.sensitivity = (5/7)*1e8;
        R.constant = R.sensitivity * R.A0;
    else
        R = singleSNCLFreqResponse(thisSNCL,t1,t2,0,Fs,units);
    end

    zeroes = R.zeros;
    poles = R.poles;

    constant = R.constant;
    if strcmp(thisS(1).kstnm,"NAS2") || strcmp(thisS(1).kstnm,"PINO") || ...
            strcmp(thisS(1).kstnm,"JUA2") || strcmp(thisS(1).kstnm,"PITA")
        constant = -constant;
    end
    Gain = R.gain;
    Sensitivity = R.sensitivity;
    Normalization = R.A0;
    stla_ = R.Stla;
    stlo_ = R.Stlo;
    stel_ = R.Stel;
    modelName_ = R.instType;

    %%
    if WAFLAG
        if ~DECONFLAG
            fprintf("forcing decon option\n");
            DECONFLAG = true;
        end
        poles = [poles wazeros];    %#ok<AGROW>
        zeroes = [zeroes wapoles];  %#ok<AGROW>
        constant = constant/wasensitivity;
    end

    %%
    for tt = 1:sumGood
        thisTrace = thisS(validTracesI(tt));
        gapFlag = thisTrace.gapFlag;
        if gapFlag
            gapInfo = thisTrace.gapInfo;
            thisTrace = nanGapWaveforms(thisTrace,0); % clobber trace, replaces all gaps with 0's
            thisTrace = transfer(thisTrace,zeroes,poles,constant,lfc,hfc,...
                npoles,DECONFLAG,false,CANUSEGPU,zeroPhaseFlag); % 1 = deconvolve
            thisTrace.gapFlag = true;
            thisTrace.gapInfo = gapInfo;
        else
            thisTrace = transfer(thisTrace,zeroes,poles,constant,lfc,hfc,...
                npoles,DECONFLAG,false,CANUSEGPU,zeroPhaseFlag); % 1 = deconvolve
        end

        if newFs
            thisTrace = resampleWaveforms(thisTrace,newFs);
        end

        thisS(validTracesI(tt)) = thisTrace;
        thisS(validTracesI(tt)).stla = stla_;
        thisS(validTracesI(tt)).stlo = stlo_;
        thisS(validTracesI(tt)).stel = stel_;
        thisS(validTracesI(tt)).A0 = Normalization;
        thisS(validTracesI(tt)).poles = poles;
        thisS(validTracesI(tt)).zeroes = zeroes;
        thisS(validTracesI(tt)).constant = constant;
        thisS(validTracesI(tt)).gain = Gain;
        thisS(validTracesI(tt)).sensitivity = Sensitivity;
    end

    %%
    modelNames(thisI) = modelName_;
    S(thisI) = thisS;
end
S = reshape(S,sizeS);