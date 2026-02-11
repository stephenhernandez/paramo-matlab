clear; close all;
dayStart = datetime(2024,10,03);
dayEnd = datetime(2024,10,03);
dayInc = 1;
dayVec = (dayStart:dayInc:dayEnd)';
lDays = length(dayVec);

%
interp_method = "spline";
npoles = 4;
zeroPhaseFlag = true;
secDur = 2*10.24;
overlapPrcnt = 3/4;
lfc = 1/secDur/2;
hfc = -inf;
newFs = 100;
vecN = 2^13;
fVec = logspace(log10(lfc),log10(newFs/2),vecN)';

%
nw = 2;
totN = secDur*newFs;
[dpssSeq,lambda] = dpss(totN,nw,2*nw-1);
w = lambda./sum(lambda);
verboseFlag = false;
nfft = 2^nextpow2(2*totN-1);
detrendFlag = true;
deconOrig = NaN(nfft/2,24*lDays);
phaseOrig = deconOrig;
mscOrig = deconOrig;
decon_t = NaN(nfft,24*lDays);
tFilt = NaT(24*lDays,1);
[Hbu,Fbu] = freqOperator(nfft,lfc,hfc,newFs,npoles);
wl = 0;
CANUSEGPU = false; %canUseGPU();

%
tic;
n = 0;
for i = 1:lDays
    day_ = dayVec(i);
    %S = loadWaveforms(day_,dayInc,"SN14",["HAE";"HHE"],"EC","",0,0,"~/masa/backups"); S(2) = differentiateWaveforms(S(2));
    %S = loadWaveforms(day_,dayInc,"VC1",["HHN";"HHZ"],"EC",""); %S = [S; S];
    %S = loadWaveforms(day_,dayInc,"PINO",["SHE";"SHN"],"EC",""); %S = [S; S];
    %S = loadWaveforms(day_,dayInc,["FARN";"FARS"],"HHZ","EC","",0,0,"~/masa/backups");
    %S = loadWaveforms(day_,dayInc,["SN14";"VCH1"],"HHZ","EC","",0,0,"~/masa/backups");
    %S = loadWaveforms(day_,dayInc,"URBI",["HHN";"HHZ"],"EC","",0,0,"~/masa/backups");
    %S = loadWaveforms(day_,dayInc,"VCH1",["HHN";"HHZ"],"EC","",0,0,"~/masa/backups");
    %S = loadWaveforms(day_,dayInc,"SN12",["HHN";"HHZ"],"EC","",0,0,"~/masa/backups");

    PINO_SH = loadWaveforms(day_,dayInc,"PINO","SHN","EC",""); %,0,0,"~/masa/backups");
    refs = pull(PINO_SH,"ref");
    if sum(~isnat(refs)) < 1
        continue;
    end

    PINO_HH = loadWaveforms(day_,dayInc,"PINO",["HHE";"HHN";"HHZ"],"EC","",false,false);
    refs = pull(PINO_HH,"ref");
    if sum(~isnat(refs)) < 3
        continue;
    end

    PINO_HH = uvw2zne(PINO_HH,"reftek"); %convert UVW to _real_ ZNE
    PINO_HH = PINO_HH(2);
    S = [scaleWaveforms(PINO_SH,-1);PINO_HH];

    lS = length(S);
    if lS < 2
        fOrig_ = [];
        Sxy_dB_ = [];
        fprintf("day not good: %s\n",day_);
        continue;
    end

    %
    refs = pull(S,"ref");
    if sum(~isnat(refs)) < 2
        fOrig_ = [];
        Sxy_dB_ = [];
        continue;
    end

    S = S(1:2);
    fprintf("processing day: %s\n",day_); toc;
    tStart = dateshift(max(cat(1,S(:).ref)),"end","minute");
    tEnd = dateshift(min(cat(1,S(:).ref)+cat(1,S(:).e)),"start","minute");

    e = pull(S,"e");
    eI = e < seconds(secDur);
    if any(eI) || tStart >= tEnd
        fOrig_ = [];
        Sxy_dB_ = [];
        continue;
    end

    Sorig = detrendWaveforms(S);
    Sf = Sorig;
    Sf = syncWaveforms(Sf,false,true,true);
    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);

    refTime = dateshift(max(pull(Sf,"ref")),"start","hour");
    endTime = dateshift(min(pull(Sf,"ref")+pull(Sf,"e")),"end","hour");
    Sf = cutWaveforms(Sf,refTime,0,endTime-refTime,false,false);
    Sf(2) = transferWaveforms(Sf(2),newFs/nfft,-inf,4,newFs,"disp",true,false,...
        CANUSEGPU,zeroPhaseFlag);
    hourVec = (refTime:hours(1):endTime)';
    hourVec = hourVec(1:end-1);

    lVec = length(hourVec);
    for j = 1:lVec
        thisHour = hourVec(j);
        Sf_ = cutWaveforms(Sf,thisHour,0,hours(1),false,false);
        refs = pull(Sf_,"ref");
        goodI = ~isnat(refs);
        if sum(goodI) < 2
            continue;
        end

        %
        Sf_ = detrendWaveforms(Sf_);
        Sf_ = nanGapWaveforms(Sf_,0);
        [decon,mscohe] = fd_decon(Sf_,Hbu,totN,dpssSeq,w,nfft,...
            overlapPrcnt,detrendFlag,wl,CANUSEGPU); toc;
        if isempty(decon)
            continue;
        end

        n = n + 1;
        abs_decon = abs(decon(1:nfft/2,:));
        mscohe = mscohe(1:nfft/2,:);
        ncol = size(mscohe,2);
        mscohe2 = mscohe./sum(mscohe,2,"omitnan");
        mscohe2(1,:) = ones(1,ncol)/ncol;
        deconOrig(:,n) = sum(mscohe2.*abs_decon,2,"omitnan");
        phase_unwrap = unwrap(angle(decon(1:nfft/2,:)));
        meanx = sum(mscohe2.*cos(phase_unwrap),2,"omitnan");
        meany = sum(mscohe2.*sin(phase_unwrap),2,"omitnan");
        meanPhaseOrig = unwrap(atan2(meany,meanx));
        integer_pi = 2*fix(round((phase_unwrap - meanPhaseOrig)/pi)/2);
        new_phase = phase_unwrap - pi*integer_pi;
        phaseOrig(:,n) = mean(new_phase,2,"omitnan");
        mscOrig(:,n) = median(mscohe,2,"omitnan");
        tFilt(n) = thisHour;
        fprintf("\n");
    end
end
toc;

%%
fshift = (-nfft/2:nfft/2-1)'*(newFs/nfft);
fOrig = fftshift(fshift,1);
fCut = fOrig(1:nfft/2);
deconOrig = deconOrig(:,1:n);
phaseOrig = phaseOrig(:,1:n);
mscOrig = mscOrig(:,1:n);
tFilt = tFilt(1:n);

%%
ncol = size(mscOrig,2);
mscohe2 = mscOrig./sum(mscOrig,2,"omitnan");
mscohe2(1,:) = ones(1,ncol)/ncol;
phase_unwrap = unwrap(phaseOrig);
meanx = sum(mscohe2.*cos(phase_unwrap),2,"omitnan");
meany = sum(mscohe2.*sin(phase_unwrap),2,"omitnan");
meanPhaseOrig = unwrap(atan2(meany,meanx));
integer_pi = 2*fix(round((phase_unwrap - meanPhaseOrig)/pi)/2);
phaseOrig = phase_unwrap - pi*integer_pi;
decon = interp1(fCut,deconOrig,fVec,interp_method);
phase = interp1(fCut,phaseOrig,fVec,interp_method);

%%
LW = 3;
close all;

msc = interp1(fCut,mscOrig,fVec,interp_method);
fI = fVec >= 0.5 & fVec <= 1;
abs_decon = abs(decon);
med_abs_decon = median(abs_decon(fI,:));
mad_score = (med_abs_decon - median(med_abs_decon))./mad(med_abs_decon,1);
mI = abs(mad_score) <= 1;

msc = msc(:,mI);
ncol = size(msc,2);
mscohe2 = msc./sum(msc,2,"omitnan");
mscohe2(1,:) = ones(1,ncol)/ncol;

figure();
semilogy(tFilt,med_abs_decon,'o'); zoom on; grid on; hold on;
semilogy(tFilt(mI),med_abs_decon(mI),'p'); zoom on; grid on; hold on;

figure();
nSubplots = 2;
tile_spacing = "compact";
tiledlayout(nSubplots,1,"Padding","compact","TileSpacing",tile_spacing);
ax = gobjects(nSubplots,1);
ax(1) = nexttile();
loglog(fVec,abs_decon(:,mI),'-',"Color",[0.5,0.5,0.5]); zoom on; grid on; hold on;
loglog(fVec,median(abs_decon(:,mI),2,"omitnan"),'-',"linewidth",LW); zoom on;
loglog(fVec,sum(mscohe2.*abs_decon(:,mI),2,"omitnan"),'-',"linewidth",LW); zoom on;

ax(2) = nexttile();
semilogx(fVec,phase(:,mI),'-',"Color",[0.5,0.5,0.5]); zoom on; grid on; hold on;
semilogx(fVec,median(phase(:,mI),2,"omitnan"),'-',"linewidth",LW);
semilogx(fVec,sum(mscohe2.*phase(:,mI),2,"omitnan"),'-',"linewidth",LW); zoom on;
linkaxes(ax,"x");

AMP = sum(mscohe2.*abs_decon(:,mI),2,"omitnan");
AMP = zpkFilter(medfiltSH(AMP,9,true),-inf,1/101,1,1,1);
PHA = sum(mscohe2.*phase(:,mI),2,"omitnan");
PHA = zpkFilter(medfiltSH(PHA,9,true),-inf,1/101,1,1,1);
response = AMP.*exp(1j*PHA);
%sys0 = idtf([ones(1,14)],ones(1,19),0);
%sys0.Structure.Numerator.Free(1:2) = false;
G = idfrd(response,fVec,0,"FrequencyUnit","Hz");
%sys = tfest(G,19,14);

% parameters that work for HHN->SHN
sys = tfest(G,32); %19,14);
denom = sys.Denominator;
numer = sys.Numerator;
[z_,p_,k_] = tf2zpk(numer,denom);
p = p_;
si = find(real(z_),1);
z = sort(z_(si:end),"descend");
zOrig = z;
z(end-1:end) = 0;
z = sort(z);
%A0 = 148307462718; %for hhn->shn
A0 = 300934;
SENSITIVITY = round(abs(k_)/A0);
[b,a] = zp2tf(z,p,1);
sys = idtf(SENSITIVITY*A0*b(find(b,1):end),a);
CONSTANT = SENSITIVITY*A0;

[mag1,phase1,wout1] = bode(sys,2*pi*fVec);
figure('units','normalized','outerposition',[0 0.1 0.5 0.9]);
nSubplots = 2;
tile_spacing = "compact";
tiledlayout(nSubplots,1,"Padding","compact","TileSpacing",tile_spacing);
axP = gobjects(nSubplots,1);
axP(1) = nexttile();
loglog(fVec,AMP,'-',"linewidth",LW); zoom on; hold on;
loglog(wout1/2/pi,squeeze(mag1),'-',"linewidth",LW); zoom on;

axP(2) = nexttile();
semilogx(fVec,PHA,'-',"linewidth",LW); zoom on; hold on;
semilogx(wout1/2/pi,-4*pi + pi*squeeze(phase1)/180,'-',"linewidth",LW); zoom on;
linkaxes(axP,"x");

figure(2);
loglog(ax(1),wout1/2/pi,squeeze(mag1),'-',"linewidth",LW); zoom on;
semilogx(ax(2),wout1/2/pi,pi*squeeze(phase1)/180,'-',"linewidth",LW); zoom on;

%sys
[p,z] = pzmap(sys)
[H,f] = cmplxResp(2^15,z,p,CONSTANT,100); figure(); loglog(f,abs(H)); zoom on;
[Hinv,f] = cmplxResp(2^15,p,z,1e9/CONSTANT,100); figure(); loglog(f,abs(Hinv)); zoom on; %VEL
[HinvD,f] = cmplxResp(2^15,p,[0;z],1e9/CONSTANT,100); figure(); loglog(f,abs(HinvD)); zoom on; %DISP