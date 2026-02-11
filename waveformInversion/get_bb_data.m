function [S,Hdecon] = get_bb_data(stnm,component,yyyy,mm,days,units,finalFs,lfc,hfc,npoles,t1,noise,signal,knetwk)
if nargin < 14; knetwk = 'EC'; end

if nargin < 11 || isempty(t1) || isempty(noise) || isempty(signal)
    cutFlag = false;
else
    cutFlag = true;
end

S = [];
if finalFs < 1
    disp('final frequency must be greater than 1')
    return
end

bbdata = load('~/igdata/ecuador_bb_meta_data');
bbkstnm = bbdata.bbkstnm;
bbstla = bbdata.bbstla;
bbstlo = bbdata.bbstlo;
bbchantype = bbdata.bbchantype;
bbdigitizer = bbdata.bbdigitizer;
bbsensor = bbdata.bbsensor;

[lia,locb] = ismember(stnm,bbkstnm);
if lia
    stnm = bbkstnm{locb};
    if strcmp(stnm,'ISPT_')
        locID = '00';
    elseif strcmp(stnm,'PDNS') || strcmp(stnm,'LGCB') || strcmp(stnm,'CABP') %|| strcmp(stnm,'PTGL')
        locID = '01';
    elseif strcmp(stnm,'PAYG')
        locID = '10';
    else
        locID = '';
    end
    
    chantypeCode = bbchantype(locb);
    if chantypeCode == 1
        chantype = 'HH';
    else
        chantype = 'BH';
    end
    chan = strcat(chantype,component);
    digitizerCode = bbdigitizer(locb);
    sensorCode = bbsensor(locb);
    digitizerGain = readDigitizerGain(digitizerCode);
    [pole,zero,sensitivity,normalizationFactor] = readSensorPolesZeros(sensorCode);
    constant = sensitivity*digitizerGain*normalizationFactor;
    %S = loadSacData(stnm,chan,yyyy,mm,days,knetwk,locID);
    S = loadSacData(datetime(yyyy,mm,days),1,stnm,chan,knetwk,locID);
    S.stla = bbstla(locb);
    S.stlo = bbstlo(locb);
    d = S.d;
    d = detrend(d);
    d = demean(d);
    d = detrend(d);
    Fs = round(1/S.delta);
    
    if strcmpi(units,'acc') || strcmpi(units,'acceleration')
        zero = zero(3:end);
    elseif strcmpi(units,'vel') || strcmpi(units,'velocity')
        zero = zero(2:end);
    end
    
    d = detrend(d);
    d = demean(d);
    d = detrend(d);
    
    
    d = resample(d,finalFs,Fs);     %upsample
    Fs = finalFs;                   %Fs*finalFs;
    
    disp('removing instrument response...');
    newConstant = 1/constant;
    nfft = length(d);
    [nfft,d] = check_nfft(nfft,d);
    Hdecon = cmplxResp(nfft,pole,zero,newConstant,Fs);
    D = fft(d);
    norig = length(d);
    
    disp(' ')
    disp(['Total numberof points to process: ',num2str(norig)])
    disp(' ')
    
    if any(isfinite([lfc hfc]))
        Hbu = freqOperator(norig,lfc,hfc,Fs,npoles);
        Hdecon = Hdecon.*Hbu;
    else
        disp('no filtering requested');
    end
    D = D.*Hdecon;
    d = ifft(D,'symmetric');
    S.d = d;
    S.npts = norig;
    S.delta = 1/Fs;
    S.e = seconds((S.npts-1)*S.delta);
    
    if cutFlag
        t2 = S.ref + S.e;
        if t1 > t2
            disp('Cannot cut data. tstart surpasses file size')
            return;
        else
            S = cutSacData(S,t1,noise,signal,true,0,true);
        end
    end
else
    disp('something is weird');
end