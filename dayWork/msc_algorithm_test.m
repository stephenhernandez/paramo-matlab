%msc_algorithm_test
clear; close all;
dayStart = datetime(2026,01,28);
%dayStart = datetime(2023,05,12);
%dayStart = datetime(2018,06,26);
dayInc = 1;
interp_method = "spline";
npoles = 4;
zeroPhaseFlag = true;
secDur = 40;
overlapPrcnt = 0.5;
lfc = 1/secDur/2;
hfc = -inf;
newFs = 50;
nw = 4;
nTapers = 2*nw - 1;
totN = secDur*newFs;
[dpssSeq,lambda] = dpss(totN,nw,nTapers);
dpssSeq_(:,1,:) = dpssSeq;
dpssSeq = dpssSeq_;
nfft = 2^(nextpow2(2*totN-1)+1);
detrendFlag = true;

Hbu1 = freqOperator(nfft,lfc,hfc,newFs,npoles);
if zeroPhaseFlag
    Hbu1 = Hbu1.*conj(Hbu1);
end
GAINCORRECTION = false;
HILBERTFLAG = false;
CANUSEGPU = false;
%S = loadWaveforms(dayStart,dayInc,"BTAM",["BHZ";"BHN";"BHE"],"EC","",0,0); %,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,"BREF",["BHZ";"BHN";"BHE"],"EC","",0,0); %,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,["BREF";"BTAM";"BVC2"],"BHZ","EC","",0,0); %,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,"TAM4",["DPZ";"DPN";"DPE"],"EC","",0,0,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,"SN12",["HHZ";"HHN";"HHE"],"EC","",0,0,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7"],["DPZ"],"EC","",0,0,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,["TAM1";"TAM2";"TAM3";"TAM4";"TAM5"],"DPZ","EC","",0,0,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,["CTX1";"CTX6"],"DPZ","EC","",0,0,"~/masa/backups");
S = loadWaveforms(dayStart,dayInc,["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"COSE";"SUCR";"SLOR"],["HHZ";"BHZ"],"EC","",0,0);
%S = loadWaveforms(dayStart,dayInc,["CTX1";"CTX2";"CTX3";"CTX4";"CTX5";"CTX6";"CTX7"],"DPZ","EC","",0,0,"~/masa/backups");
%S = loadWaveforms(dayStart,dayInc,"CTX1",["DPZ";"DPN";"DPE"],"EC","",0,0,"~/masa/backups");

Sf = detrendWaveforms(S);
Sf = syncWaveforms(Sf,false,true,true);
Sf = resampleWaveforms(Sf,newFs);
Sf = filterWaveforms(Sf,lfc,-inf,npoles,[],zeroPhaseFlag);
Sf = nanGapWaveforms(Sf,0);
Sf = padWaveforms(Sf);
fshift = (-nfft/2:nfft/2-1)'*(newFs/nfft);
fOrig = fftshift(fshift,1);
fCut = fOrig(1:nfft/2);

%%
tic; [Amp,SpecOrig,Phase] = spectral_estimation(Sf,Hbu1,...
    totN,dpssSeq,nfft,lambda,overlapPrcnt,detrendFlag,...
    GAINCORRECTION,HILBERTFLAG,CANUSEGPU); toc;

%%
[Hbu,Fbu] = freqOperator(nfft,lfc,hfc,newFs,npoles);
Hbu = Hbu.*conj(Hbu);
lS = length(Sf);
NWIN = 2*3600/secDur;
NHOUR = 3600/secDur;
Nfilt = 5;
close all;
n = 0;
mean_msc = zeros(size(Amp(:,:,1)));
[stla,stlo] = metaDataFromStationList(pull(Sf,"kstnm"));
refEllipse = referenceEllipsoid('wgs84');
[d_,azs] = distance(-0.684099,-78.436745,stla,stlo,refEllipse);
for j = 1:lS-1
    kstnm1 = Sf(j).kstnm;
    dist1 = d_(j);
    Amp1 = Amp(:,:,j);
    for k = j+1:lS
        n = n+1;
        fCut_ = fCut;
        dist2 = d_(k);
        kstnm2 = Sf(k).kstnm;
        Amp2 = Amp(:,:,k);
        fprintf("%d %d %d\n",n,j,k);
        deconvolution = SpecOrig(:,:,n)./Amp1;
        xcoherence = deconvolution.*sqrt(Amp1);
        xcoherence = xcoherence./sqrt(Amp2);
        msc = xcoherence.*conj(xcoherence);
        msc(~isfinite(msc)) = 0;
        % % %msc = atanh(msc);
        % % %msc_f = detrend(diff(msc)); fCut_ = fCut_(1:end-1);
        % % %msc_f = detrend(diff(real(xcoherence))); fCut_ = fCut_(1:end-1);
        % % msc_f = Hbu.*[msc; msc(1,:); flipud(conj(msc(2:end,:)))];
        % % msc_f = fftshift(ifft(msc_f,[],1,"symmetric"),1);
        % % %msc_f = msc_f(1:nfft/2,:);
        % % msc_f = normalizeWaveforms(msc_f);
        % % msc_f = (((pws(msc_f,~true,~true,1,NWIN))));
        % % d2 = zpkFilter(msc_f(:,NWIN:NHOUR:end),-inf,1/Nfilt,1,1,1);
        % % %d2 = tanh(d2);
        % % %axL = linkedPlot([],(d2),"-","none"); zoom on;
        % % %axL = linkedPlot(fCut_,(d2),"-","none"); zoom on;
        % % % for i = 1:length(axL)
        % % %     axL(i).XScale = "log"; xlim([0.1 10]);
        % % % end
        % % %sgtitle(sprintf("pair: %d - %d",j,k));
        % % %xlim([0.2 2])
        % % %close all;
        figure();
        %msc = tanh(msc);
        imagesc((1:size(msc,2)-1)',fCut,10*log10(msc)); zoom on; colorbar; axis xy;
        clim([-10 0]);
        title(sprintf("pair: %s (%f) - %s (%f), interstation distance: %f",kstnm1,dist1,kstnm2,dist2,dist2-dist1));
        mean_msc = mean_msc + atanh(msc);
    end
end
mean_msc = mean_msc/n;

%%
figure();
mean_msc = tanh(mean_msc);
imagesc((1:size(mean_msc,2)-1)',fCut,10*log10(mean_msc)); zoom on; colorbar; axis xy;
clim([-10 0]);
title(sprintf("mean msc - all %d pairs",n));




