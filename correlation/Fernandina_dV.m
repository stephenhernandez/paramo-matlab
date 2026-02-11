% function [CCpairs,kstnms,elementLats,elementLons,elementElevs,newFs,shortWin,...
%     longWin,nOverlap1,nOverlap2,tMain,fshift,fOrig,SW,Ddiag] = ...
%     ArrayAndNetworkSVD(dayVec)
clear; close all; clc;
prewhiten = true;

%Long Window seismic
newFs = 10;
shortWin = 60;  nOverlap1 = 1/2;
longWin = 600;  nOverlap2 = 9/10;

refEllipse = referenceEllipsoid('wgs84');
lfc = 1/(shortWin);
dW = 0; % window length (in Hz) for smoothing spectrum
detrendFlag = true;
gapFillValue = 0;
hfc = -inf; %floor(0.8*0.5*newFs);
npoles = 4;
lfc2 = 0.2;
hfc2 = 0.8;

kstnms = ["CBDG"; "FER3"];
transferFlag = false;
dayStart = datetime(2024,01,01);
dayEnd = datetime(2024,07,31);
dayInc = 1;
dayVec = (dayStart:dayInc:dayEnd)';
lDays = length(dayVec);

%%
Hd = zpkOperator(lfc2,hfc2,newFs,npoles);
winlen = shortWin*newFs;
nfft = 2^(nextpow2(winlen)+1);

tw = shortWin*newFs*10;
df = NaN(nfft,lDays);

%%
for l = 1:lDays
    day_ = dayVec(l);
    C = loadWaveforms(day_,dayInc,kstnms,"HHZ","EC","",true,true,'~/data/Galapagos_August2024_ServiceRun');
    kstnms = pull(C,'kstnm');

    tDayStart = dateshift(C(1).ref,'start','day');
    tBegs = pull(C,'ref');
    tStops = tBegs + pull(C,'e');
    tStart = dateshift(max(tBegs),'end','minute');
    tEnd = dateshift(min(tStops),'start','minute');
    Ccut = cutWaveforms(C,tStart,0,tEnd-tStart);

    %

    Cf = detrendWaveforms(...
        (...
        detrendWaveforms(...
        taperWaveforms(...
        interpolateWaveforms(...
        resampleWaveforms(...
        detrendWaveforms(...)
        (Ccut)),...
        newFs*1),gapFillValue),tw))));
    Cf = padWaveforms(Cf);

    if ~prewhiten
        Cf = normalizeWaveforms(Cf,detrendFlag,false,true);
    end

    %%
    if transferFlag
        Cf = resampleWaveforms(Cf,newFs);
        Cf = padWaveforms(Cf);
        Cf = nanGapWaveforms(Cf,gapFillValue);
        Cf = transferWaveforms(Cf,lfc,hfc,npoles,newFs,"vel",1,false);
        Cf = scaleWaveforms(Cf,1e9);
        Cf = padWaveforms(Cf);
    end

    lC = length(Cf);
    nXpairs = 0.5*lC*(lC-1);
    numCombos = lC + nXpairs;

    %% get data cube
    dCube = [];
    for i = 1:lC
        d_ = Cf(i).d;
        dcut_ = cutWindows(d_,longWin*newFs,nOverlap2,detrendFlag);
        dCube = cat(3,dCube,dcut_);
    end

    [~,~,endIndex] = cutWindows(d_,longWin*newFs,nOverlap2,detrendFlag);
    tMain = getTimeVec(Cf(1));
    tMain = tMain(endIndex);
    [~,lT,lC] = size(dCube);

    %w = rectwin(winlen);
    %w = kaiser(winlen,10);
    w = blackmanharris(winlen); %
    %w = tukeywin(winlen,2*(1-nOverlap1)); %

    fshift = (-nfft/2:nfft/2-1)'*(newFs/nfft);
    fOrig = fftshift(fshift,1);
    Ddiag = NaN(nfft,lT);
    SW = Ddiag;

    CCpairs = NaN(nfft*1,lT,numCombos);
    D3 = NaN(nfft,numCombos);

    %%
    for i = 1:lT
        disp(i);
        d_ = squeeze(dCube(:,i,:));
        dCube2 = [];
        for kk = 1:lC
            d__ = d_(:,kk);
            if prewhiten
                d__ = tdNorm(d__,shortWin,1,newFs);
                d__ = fdWhiten(d__,lfc,hfc,dW,newFs,true,true);
            end
            dcut_ = cutWindows(d__,shortWin*newFs,nOverlap1,~detrendFlag);
            dcut = normalizeWaveforms(dcut_,true);
            dcut_ = dcut_.*w;
            D_ = fft(dcut_,nfft);
            dCube2 = cat(3,dCube2,D_);
        end

        %% now loop through each frequency
        svdFlag = true;
        lT2 = size(dCube2,2);
        for j = 1:nfft
            W = squeeze(dCube2(j,:,:))'; %nRows = nChans, nColumns = nTimeSlices
            [Uorig,Ddiag_,~] = svd(W,'econ');
            Ddiag_ = diag(Ddiag_);
            SW_ = 1-((2*Ddiag_(end)^2)/(Ddiag_(1)^2 + Ddiag_(2)^2));  %planarity
            %SW_ = 1-(Ddiag_(3)^2/Ddiag_(1)^2);                     %linearity (v1)
            %SW_ = 1-(Ddiag_(end)^2/Ddiag_(1)^2);                   %linearity (v2)
            %SW_ = sum((0:lC-1)'.*Ddiag_)./sum(Ddiag_);
            SW(j,i) = SW_;
            Ddiag(j,i) = Ddiag_(1);
            U = Uorig(:,1);
            CC = U*U';
            triuCC = triu(CC)';
            diagtriuCC = diag(triuCC);
            triuCC = triuCC - diag(diagtriuCC);

            D3(j,1:nXpairs) = triuCC(triuCC ~= 0)';
            D3(j,nXpairs+1:end) = diagtriuCC';
        end
        CCpairs(:,i,:) = D3;
    end

    %%
    df_ = ifft(CCpairs(:,1),[],1,'symmetric');
    df_ = fftshift(df_,1);
    df_ = filter(Hd,df_);
    df_ = flipud(df_);
    df_ = filter(Hd,df_);
    df_ = flipud(df_);
    df(:,l) = df_;
end


