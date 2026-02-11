clear; close all; %clc;

%%
tic;
S(1) = readMiniSeed('II.MBAR.00.BHZ.2021.05.20.mseed');
S(2) = readMiniSeed('II.MBAR.00.BHZ.2021.05.21.mseed');
S(3) = readMiniSeed('II.MBAR.00.BHZ.2021.05.22.mseed');
S(4) = readMiniSeed('II.MBAR.00.BHZ.2021.05.23.mseed');
S(5) = readMiniSeed('II.MBAR.00.BHZ.2021.05.24.mseed');
S(6) = readMiniSeed('II.MBAR.00.BHZ.2021.05.25.mseed');
S(7) = readMiniSeed('II.MBAR.00.BHZ.2021.05.26.mseed');
toc;

%%
S = mergeWaveforms(S);
S = resampleWaveforms(S,20);
toc;

%%
close all;
Sf = filterWaveforms(S,0.75,3);
df = Sf.d;
tf = getTimeVec(Sf);
toc;

Scut = detrendWaveforms(cutWaveforms(Sf,datetime(2021,05,23,11,37,00),0,minutes(1)));
tmplate = Scut.d;
tmplate = flipud(tmplate/norm(tmplate));
winlen = length(tmplate);
box = ones(winlen,1);
toc;

cc = fftfilt(tmplate,df);
norms = fftfilt(box,df.^2);
norms = sqrt(abs(norms));
ccnorm = cc./norms;
toc;


thresh = 0.3;
[pks,locs] = findpeaks(ccnorm,'MINPEAKHEIGHT',thresh,'MINPEAKDISTANCE',winlen*2);
locsI = locs > winlen;
toc;

%%
if ~sum(locsI)
    fprintf('no matches found\n');
    return;
end

locs = locs(locsI) - winlen + 1;
pks = pks(locsI);

%%
tmatch = tf(locs);

figure();
plot(tmatch,1:length(tmatch),'.'); zoom on; grid on;

%%
Scut = cutWaveforms(Sf,tmatch,0,seconds(60),true);
dcut = pull(Scut);

