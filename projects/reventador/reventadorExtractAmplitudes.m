clear; close all; clc;
load('~/research/now/reventador/ReventadorSubspaceDetectorResults_v10','tabs');
%tabs = tabs(tabs >= datetime(2018,07,18) & tabs < datetime(2018,07,19));
tabs = tabs(tabs >= datetime(2010,01,01) & tabs < datetime(2023,01,01));
%tabs = tabs(tabs >= datetime(2016,08,01) & tabs <= datetime(2019,04,01));

%
[N,edges] = histcounts(tabs,dateshift(min(tabs),'start','day'):days(1):dateshift(max(tabs),'end','day'));
edges = edges(1:end-1)';
[N,sortI] = sort(N,'descend');
N = N';
edges = edges(sortI);
load('~/research/now/reventador/CASC_BONI_basisFunctions_29DEC2021','kstnm','lfc','hfc');
uniqKstnms = unique(kstnm);

%
lIndiv = length(tabs);
ei = cumsum(N);
si = 1+cumsum([0; N]);
si = si(1:length(ei));
lN = length(N);
lAmps = 24;

%
[stla,stlo] = metaDataFromStationList(uniqKstnms);
refEllipse = referenceEllipsoid('wgs84');
[d_,azs] = distance(stla,stlo,-0.080850,-77.657995,refEllipse);
[d_,sortI] = sort(d_);
azs = azs(sortI);
uniqKstnms = uniqKstnms(sortI);

% import amplification curves
amplificationCurves = [];
load('~/research/now/hvsr/casc/casc2018.mat','hvMean','freqVec');
amplificationCurves = [amplificationCurves hvMean];
load('~/research/now/hvsr/boni/boni2018.mat','hvMean');
amplificationCurves = [amplificationCurves hvMean];
load('~/research/now/hvsr/ants/ants2018.mat','hvMean');
amplificationCurves = [amplificationCurves hvMean];
load('~/research/now/hvsr/antg/antg2018.mat','hvMean');
amplificationCurves = [amplificationCurves hvMean];

%
secDur = 120;
npoles = 4;
newFs = 12;
maxN = max(N);

%%
units = 'disp';
waFlag = true;
lUniqKstnms = length(uniqKstnms);
dayLoopN = 2; %lN;
tStarts = NaT(1,lUniqKstnms);
tEnds = NaT(dayLoopN,lUniqKstnms);
for i = 1%:lUniqKstnms
    tStarts(i) = dn2dt(now);
    stnm_ = uniqKstnms(i);
    amplificationCurve_ = [freqVec amplificationCurves(:,i)];
    az_ = azs(i);

    tmpAmps = NaN(maxN*lAmps,dayLoopN);
    lT = NaN(dayLoopN,1);

    for j = 1:dayLoopN
        day_ = edges(j);
        tI = tabs >= day_ & tabs < day_+1;
        t_ = tabs(tI);
        lT_ = sum(tI);
        lT(j) = lT_;

        %
        SZ_ = loadWaveforms(day_,1,stnm_,"HHZ");
        SN_ = loadWaveforms(day_,1,stnm_,"HHN");
        SE_ = loadWaveforms(day_,1,stnm_,"HHE");
        S_ = [SZ_; SN_; SE_];
        if any(isnat(pull(S_,'ref')))
            fprintf('no or incomplete data from %s on day %s, (%d)\n',stnm_,datestr(day_),j);
            tEnds(j,i) = dn2dt(now);
            continue;
        end

        %
        S_ = detrendWaveforms(S_);
        S_ = interpolateWaveforms(S_);
        S_ = syncWaveforms(S_);

        %
        Sorig = filterWaveforms(S_,lfc,hfc);
        S_ = transferWaveforms(S_,lfc,hfc,npoles,newFs,units,waFlag);
        %S_ = detrendWaveforms(S_);

        %
        SZ_ = S_(1);
        SN_ = S_(2);
        SE_ = S_(3);
        Srot = rotateWaveforms([SN_; SE_],az_,2,1,true);

        NEcorr = siteAmplificationCorrection([SN_;SE_],amplificationCurve_);
        Rotcorr = siteAmplificationCorrection(Srot,amplificationCurve_);

        tw = 0.08;
        try
            SZ = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(SZ_,t_,0,seconds(secDur))),tw)));
            SN = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(SN_,t_,0,seconds(secDur))),tw)));
            SE = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(SE_,t_,0,seconds(secDur))),tw)));

            Ncorr = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(NEcorr(1),t_,0,seconds(secDur))),tw)));
            Ecorr = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(NEcorr(2),t_,0,seconds(secDur))),tw)));
            R = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Srot(1),t_,0,seconds(secDur))),tw)));
            T = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Srot(2),t_,0,seconds(secDur))),tw)));
            Rcorr = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Rotcorr(1),t_,0,seconds(secDur))),tw)));
            Tcorr = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Rotcorr(2),t_,0,seconds(secDur))),tw)));

            SZorig = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Sorig(1),t_,0,seconds(secDur))),tw)));
            SNorig = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Sorig(2),t_,0,seconds(secDur))),tw)));
            SEorig = double(pull(taperWaveforms(detrendWaveforms(cutWaveforms(Sorig(3),t_,0,seconds(secDur))),tw)));

            amps_ = NaN(maxN,lAmps);
            amps_(1:lT_,1) = 0.5.*peak2peak(SZ)';
            amps_(1:lT_,2) = 0.5.*peak2peak(SN)';
            amps_(1:lT_,3) = 0.5.*peak2peak(SE)';
            amps_(1:lT_,4) = 0.5.*peak2peak(Ncorr)';
            amps_(1:lT_,5) = 0.5.*peak2peak(Ecorr)';
            amps_(1:lT_,6) = 0.5.*peak2peak(R)';
            amps_(1:lT_,7) = 0.5.*peak2peak(T)';
            amps_(1:lT_,8) = 0.5.*peak2peak(Rcorr)';
            amps_(1:lT_,9) = 0.5.*peak2peak(Tcorr)';
            amps_(1:lT_,10) = 0.5.*peak2peak(SZorig)';
            amps_(1:lT_,11) = 0.5.*peak2peak(SNorig)';
            amps_(1:lT_,12) = 0.5.*peak2peak(SEorig)';

            amps_(1:lT_,13) = rms(SZ)';
            amps_(1:lT_,14) = rms(SN)';
            amps_(1:lT_,15) = rms(SE)';
            amps_(1:lT_,16) = rms(Ncorr)';
            amps_(1:lT_,17) = rms(Ecorr)';
            amps_(1:lT_,18) = rms(R)';
            amps_(1:lT_,19) = rms(T)';
            amps_(1:lT_,20) = rms(Rcorr)';
            amps_(1:lT_,21) = rms(Tcorr)';
            amps_(1:lT_,22) = rms(SZorig)';
            amps_(1:lT_,23) = rms(SNorig)';
            amps_(1:lT_,24) = rms(SEorig)';


            tmpAmps(:,j) = amps_(:);
            fprintf('%s: <strong>%s</strong> (%d)\n',stnm_,datestr(day_),j);
            tEnds(j,i) = dn2dt(now);
        catch
            fprintf(2,'AMBIGUOUS ERROR: %s, <strong>%s</strong> (%d)\n',stnm_,datestr(day_),j);
            tEnds(j,i) = dn2dt(now);
            continue;
        end
    end

    tmpAmps2 = NaN(lIndiv,lAmps);
    for j = 1:dayLoopN
        lT_ = lT(j);
        si_ = si(j);
        ei_ = ei(j);
        amps_ = tmpAmps(:,j);
        amps_ = reshape(amps_,[maxN, lAmps]);
        amps_ = amps_(1:lT_,:);
        tmpAmps2(si_:ei_,:) = amps_;
    end

    if i == 1
        allCASCAmps = tmpAmps2;
        cd ~/research/now/reventador/
        save('correctedAmplitudesCASC_WA','allCASCAmps','lN','lT','tabs');
    elseif i == 2
        allBONIAmps = tmpAmps2;
        cd ~/research/now/reventador/
        save('correctedAmplitudesBONI_WA','allBONIAmps','lN','lT','tabs');
    elseif i == 3
        allANTSAmps = tmpAmps2;
        cd ~/research/now/reventador/
        save('correctedAmplitudesANTS_WA','allANTSAmps','lN','lT','tabs');
    else
        allANTGAmps = tmpAmps2;
        cd ~/research/now/reventador/
        save('correctedAmplitudesANTG_WA','allANTGAmps','lN','lT','tabs');
    end
end
