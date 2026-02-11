%clear;
%close all; clc;
cd ~/masa/supermouv/profil010/

%
suFiles = dir('*.su');
lFiles = length(suFiles);
secDur = 3;
maxSeconds = 5;
FsOrig = 2000;
maxWinlen = 1 + FsOrig*secDur;
nTraces = 96;
bigMatrix = NaN(maxWinlen,nTraces,lFiles);

Nmed = 1;
tw = 40;
dx = 6.25; %spacing between hydrophones
dx2 = dx/2;
Dshot2hydrophone1 = dx; %this is an assumption, i have no clue if its right!

npoles = 4;
dumbTraces = 4;
cmap = [bone(256); sky(256)];
zeroPhaseFlag = false;
clipLevel = 2e0;
upsampleFactor = 1;

lfc1 = 20;
hfc1 = 320;

close all;
dtShot = 6; % 6 seconds between shots
traceToShift = 1; %variable, depends on ship velocity
nMidpoints = nTraces;
midpointIndex = (1:nMidpoints)';
lastMidpointIndex = max(midpointIndex);
close all
for i = 620%:5:lFiles %loop through each shot
    disp(i);
    thisFile = char(suFiles(i).name);
    [Data,SegyTraceHeaders,SegyHeader]=ReadSu(thisFile,'endian','b');
    Data = Data(:,dumbTraces+1:dumbTraces+nTraces);
    [nT,nR] = size(Data);
    nw = 4;
    totN = nT;
    [dpssSeq,lambda] = dpss(totN,nw,2*nw-1);
    w = lambda./sum(lambda);
    lS = length(w);
    Data = taper(detrend(Data),4*FsOrig/lfc1);
    Data = normalizeWaveforms(Data);
    nfft = 2*nT-1;

    [Hbu,Fbu] = freqOperator(nfft,lfc1,hfc1,FsOrig,npoles);
    Data = repmat(Data,1,1,lS);
    for j = 1:lS
        dpssSeq_ = dpssSeq(:,j);
        Data(:,:,j) = dpssSeq_.*Data(:,:,j);
    end

    D = fft(detrend((Data)),nfft,1);
    Syy = sum(D.*conj(D),3);
    Syy = ones(size(Syy));
    Sxy = sum(D.*conj(D),3); %cross spectrum (cross correlation)

    D = Hbu.*(Sxy./Syy).*conj(Hbu);
    d3 = ifft(D,[],1,"symmetric");
    d3 = fftshift(taper(detrend(d3),4*FsOrig/lfc1),1);
    d3 = abs(hilbert(normalizeWaveforms((tdNorm(detrend(diff(d3)),0.1,1,FsOrig)))));
    d3 = flipud(d3(1:10000,:));
    tmp_stack = plot_family(d3,1:size(d3,2),40,FsOrig);
    D = fft(detrend(d3),nfft,1);

    autoconv = ifft(D.*conj(D),[],1,"symmetric");
    %autoconv = ifft(D.*D,[],1,"symmetric");
    %autoconv = fftshift(ifft(D,[],1,"symmetric"),1);
    reflection_profile = tdNorm(d3,0.1,1,FsOrig); %d3;
    figure('units','normalized','outerposition',[0.08 0 0.2 0.95]);
    %imagesc((fliplr(wiener2(fliplr(normalizeWaveforms(wiener2(reflection_profile(1:19000,:),[41,7]))),[41,7]))));
    %imagesc(fliplr(normalizeWaveforms(wiener2(fliplr(reflection_profile(1:5000,:)),[31,7]))));
    imagesc(wiener2(reflection_profile,[15,5]));
    zoom on;
    colormap(cmap);
    title(sprintf("shot: %d",i));
    colorbar;
    %clim(0.04*[-1 1]);
    pause(0.9);
    %close all;
end

%%
d2 = detrend(diff(Data));
d2 = taper(d2,tw);
D1 = d2(:,1);
D2 = d2(:,5);
nw = 4;
[npts,ntimewin] = size(D1);
[dpsSeq,lambda] = dpss(npts,nw,2*nw-1);
w = lambda./sum(lambda);
nfft = 2^nextpow2(npts);
Hbu = freqOperator(nfft,lfc1,hfc1,FsOrig,npoles);

%%
Sxx = NaN(nfft,ntimewin);
Syy = Sxx;
Sxy = Sxx;
for i = 1:ntimewin
    dcut_1 = detrend(D1(:,i));
    dcut_1 = dcut_1.*dpsSeq;
    dcut_2 = detrend(D2(:,i));
    dcut_2 = dcut_2.*dpsSeq;
    D_1 = fft(dcut_1,nfft);
    D_2 = fft(dcut_2,nfft);
    D_1 = w'.*D_1;
    D_2 = w'.*D_2;
    Sxx_1 = sum(D_1.*conj(D_1),2);
    Sxx(:,i) = Sxx_1;
    Syy_ = sum(D_2.*conj(D_2),2);
    Syy(:,i) = Syy_;
    Sxy_ = sum(conj(D_1).*D_2,2);
    Sxy(:,i) = Sxy_;
end

%%
uOrig = ifft(conj(conj(D_2.*Hbu).*Hbu),[],1,'symmetric');
uOrig = sum(uOrig,2);

X = ifft(conj(conj(D_1.*Hbu).*Hbu),[],1,'symmetric');
X = sum(X,2);

decon = Sxy./Sxx;
g = ifft(conj(conj(decon.*Hbu).*Hbu),[],1,'symmetric');

f = g(1:npts);
f(f<0) = 0;
figure(); plot(f); zoom on;
F = fft(f.*dpsSeq,nfft);
F = w'.*F;

uEst = ifft(conj(conj((D_1.*F).*Hbu).*Hbu),[],1,'symmetric');
uEst = sum(uEst,2);
figure(); plot([X uOrig uEst]); zoom on; grid on;