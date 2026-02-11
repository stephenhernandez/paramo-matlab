clear; close all;

%%
refEllipse = referenceEllipsoid('wgs84');
[stla,stlo] = metaDataFromStationList("VCH1");
load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag');
globalt = t;
globaleqlat = eqlat;
globaleqlon = eqlon;
globaleqdepth = eqdepth;
globaleqmag = eqmag;
clear t eq*

%%
dists = distance(stla,stlo,globaleqlat,globaleqlon,refEllipse)*1e-3;

%%
%Amplitude Constants
VEL = 3.5;                                  % surface wave velocity in km/sec
Elr0 = 1e-9;                                % minimum strain threshold
X0 = log10(Elr0*20*(VEL*1e9)/(2*pi)) + 2;   % strain threshold converted to log_10(amplitude)

%Spatial Constants
km2deg = 1/111.19;                        % coarse conversion factor
Dmin = 200;                               % distance threshold in km

Idist = (dists >= Dmin);
Ithresh = (globaleqmag > X0+1.66*log10(dists*km2deg)); %the +2 is hidden in X0
Itime = globalt >= datetime(2013,07,01) & globalt <= datetime(2019,07,01);
goodI = (Idist & Ithresh & Itime);

amps = 10.^(globaleqmag - 1.66*log10(dists.*km2deg) - 2);
strain = amps/(20*(VEL*1e9)/(2*pi));
SWAT = globalt + seconds(dists./5.5);       % large amplitude arrival time
TphaseTT = seconds(dists./1.5);             % duration object

%%
cd ~/research/now/sierra_negra/
load snSubspaceDetectorVCH1_v1.mat %<-- exactly the same as the ftext file ive sent you

ampI = tabs <= datetime(2015,01,057);
z2p(ampI,:) = z2p(ampI,:)/3.141950e+08; %BHZ data
z2p(~ampI,:) = z2p(~ampI,:)/1.256780e+09; %BHZ data
%nI = true(size(tabs));
nI = ~(nanmedian(z2p./p2rms,2) < 1e-7 | tabs <= datetime(2013,05,01)); % | nanmedian(kurt,2) > 50); %<-- z2p/p2rms/kurt are Nx3 columns corresponding to Z/N/E channels

figure();
ax2 = subplot(211); semilogy(tabs(nI),nanmedian(z2p(nI,:)./p2rms(nI,:),2),'.'); zoom on; grid on; title('amplitude');
ax2(2) = subplot(212); plot(tabs(nI),1:sum(nI),'.'); zoom on; grid on; title('cumulative number');
linkaxes(ax2,'x');

%%
localt = tabs(nI);
localamp = nanmedian(z2p(nI,:),2);

%%
goodStrain = strain(goodI);
goodt = globalt(goodI);
goodSWAT = SWAT(goodI);
goodTpTT = TphaseTT(goodI);

[goodStrain,sI] = sort(goodStrain,'descend');
goodt = goodt(sI);
goodSWAT = goodSWAT(sI);
goodTpTT = goodTpTT(sI);

% figure();
% semilogy(goodSWAT,goodStrain,'o'); zoom on; grid on;


%TB = days(30);
Nwindows = 1440;
maxMinutesAfter = 30;
minSTDs = 2;

maxTriggers = min([sum(goodI) 5e2]);
beta = NaN(maxTriggers,1);
ampRatio = beta;

%%
I = [];
t = [];
St = [];
NB = [];
NA = [];
BETA = [];
TA_ = [];
TB_ = [];
SUME = [];

%%
nn = 0;
for i = 1%:maxTriggers
    TA = minutes(maxMinutesAfter);
    tStart = goodSWAT(i);
    tEnd = min([tStart + TA tStart + goodTpTT(i)]);
    TA = tEnd - tStart; %duration object
    TB = Nwindows*TA;
    
    if i > 1
        nearness = tStart - goodSWAT(1:i-1);
        if  any(nearness > 0 & nearness <= days(1/2))
            %skipping because this potential trigger is too close in time
            %to a previous trigger
            %disp(['yikes, skipping: ',datestr(goodSWAT(i))])
            continue;
        end
    end
    t2I = find(localt > tStart,1);
    t2_ = seconds(localt(t2I)-tStart);
    
    Iafter = localt >= tStart & localt <= tEnd;
    Nafter = sum(Iafter);
    if ~Nafter
        %fprintf('no events for strain: %f\n',goodStrain(i));
        continue;
    end
    Ibefore = localt >= tStart-TB & localt <= tStart;
    Nbefore = sum(Ibefore);
    if ~Nbefore
        continue;
    end
    N = Nbefore + Nafter;
    T = days(TB) + days(TA);
    TaT = days(TA)/T;
    numerator_ = Nafter - N*TaT;
    denom_ = sqrt(N*(TaT)*(1-TaT));
    beta_ = numerator_/denom_;
    
    
    if beta_ > minSTDs
        nn = nn + 1;
        beta(i) = beta_;
        ampRatio(i) = nanmean(1e6*localamp(Iafter))./(nanmean(1e6*localamp(Ibefore)));
        I = [I; i];
        t = [t; datenum(tStart)];
        St = [St; log10(goodStrain(i))];
        NB = [NB; Nbefore];
        NA = [NA; Nafter];
        BETA = [BETA; beta_];
        TA_ = [TA_; minutes(minutes(TA))];
        TB_ = [TB_; days(days(TB))];
        SUME = [SUME;ampRatio(i)];
    end
end

%%
if nn
    t = dn2dt(t);
    TTtmp = timetable(t,NA,NB,BETA,TA_,TB_,St,SUME,I,'VariableNames',...
        {'Nafter','Nbefore','Beta','Tafter','Tbefore','logStrain','AmpRatio','TriggerRank'});
    disp(TTtmp)
end

%%
figure('units','normalized','outerposition',[0 0.25 1 0.7]);
axx = subplot(131);
plot(goodt(1:maxTriggers),beta,'s'); zoom on; grid on;
axx(2) = subplot(132);
semilogx(goodStrain(1:maxTriggers),beta,'o'); zoom on; grid on;
axx(3) = subplot(133);
semilogx(ampRatio,beta,'p'); zoom on; grid on;
linkaxes(axx,'y');
