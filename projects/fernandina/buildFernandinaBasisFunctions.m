clear
close all;

%
cd ~/research/now/fernandina/
load fernandinaSubspaceDetector_v1.mat

%
goodI = NCC>=0.1 & nanmedian(z2p,2)>=200 & Neff > 1;
clearvars -except tabs goodI

%
lfc = 4;
hfc = 12;

noiseWin = 0;
secDur = 40;

tic;
FER1Z = extractWaveforms(tabs(goodI)-seconds(noiseWin),seconds(secDur),"FER1","BHZ","EC","",true,true);
FER1N = extractWaveforms(tabs(goodI)-seconds(noiseWin),seconds(secDur),"FER1","BHN","EC","",true,true);
FER1E = extractWaveforms(tabs(goodI)-seconds(noiseWin),seconds(secDur),"FER1","BHE","EC","",true,true);

%
FER2Z = extractWaveforms(tabs(goodI)-seconds(noiseWin),seconds(secDur),"FER2","HHZ","EC","",true,true);
FER2N = extractWaveforms(tabs(goodI)-seconds(noiseWin),seconds(secDur),"FER2","HHN","EC","",true,true);
FER2E = extractWaveforms(tabs(goodI)-seconds(noiseWin),seconds(secDur),"FER2","HHE","EC","",true,true);
toc;

%
tic;
newFs = 40;
FER1Z = interpolateWaveforms(resampleWaveforms(FER1Z,newFs));
FER1N = interpolateWaveforms(resampleWaveforms(FER1N,newFs));
FER1E = interpolateWaveforms(resampleWaveforms(FER1E,newFs));

FER2Z = interpolateWaveforms(resampleWaveforms(FER2Z,newFs));
FER2N = interpolateWaveforms(resampleWaveforms(FER2N,newFs));
FER2E = interpolateWaveforms(resampleWaveforms(FER2E,newFs));
toc;

%
S = [FER1Z; FER1N; FER1E; FER2Z; FER2N; FER2E];

%
S = differentiateWaveforms(S);
S = detrendWaveforms(S);

%
S = intWaveforms(taperWaveforms(detrendWaveforms(filterWaveforms(S,lfc,hfc)),0.1));
S = reshape(S,[4870 6]);
toc;

% here i need to figure out which streams are bad because they are either empty or they're full of NaNs
tic;
goodStreams = true(1,size(S,1));
for i = 1:size(S,2)
    S_ = S(:,i);
    d = double(pull(S_));
    t_isgood = ~isnat(pull(S_,'ref'));
    ii = logical(sum(~isfinite(d)) == 0 & t_isgood');
    goodStreams = goodStreams & ii;
end

%
FER1Z = S(goodStreams,1);
FER1N = S(goodStreams,2);
FER1E = S(goodStreams,3);

FER2Z = S(goodStreams,4);
FER2N = S(goodStreams,5);
FER2E = S(goodStreams,6);
toc;

%%
tic;
S = [FER1Z; FER1N; FER1E; FER2Z; FER2N; FER2E];
S = reshape(S,[length(FER1Z) 6]);
dMaster = [];
for i = 1:size(S,2)
    S_ = S(:,i);
    d = double(pull(S_));
    d = d(1:1600,:);
    d = d./rssq(d);
    d = detrend(d);
    d = taper(d,0.1);
    d = d./rssq(d);
    dMaster = [dMaster; d];
end
toc;

%%
[maxccp_,plags_,maxccn_,nlags_] = doccFreqCircShift(dMaster,true);
toc;


%%
close all; 
tic;
[family,l_uniq_indices,Nsingletons,tree] = generateFamilies(maxccp_,0.5,'weighted',5,true);
designSet = [];
winlen = 100;
maxFams = 200;
pickTimes = NaN(maxFams,6);
for i = 1:maxFams
    data = dMaster(:,family{i});
    data = apply_vdcc(data);
    data = detrend(data./rssq(data));
    tmpStack = pws(data);
    stackComponents = detrend(reshape(tmpStack,[1600,6]));
    
    for j = 1:6
        st_ = stackComponents(:,j);
        envst = zpkFilter(abs(hilbert(st_)),-inf,1/10,1,1,false);
        envst = fftfilt(ones(100,1)/100,envst);
        snr = zpkFilter(taper(envst(winlen:end)./envst(1:end-winlen+1),0.1),-inf,1/10,1,1,1);
        [~,maxI] = max(snr);
        pickTimes(i,j) = min(maxI);
    end
    designSet = [designSet normalizeWaveforms(detrend(tmpStack))];
end
toc;

%
pickTimesF = floor(nanmedian(pickTimes,2)) - 2*newFs;
dSet2 = [];
for i = 1:maxFams
    ds_ = designSet(:,i);
    ds_ = detrend(reshape(ds_,[1600,6]));
    ds_ = circshift(ds_,-pickTimesF(i),1);
    ds_ = demean(detrend(normalizeWaveforms(ds_)));
    ds_ = demean(detrend(normalizeWaveforms(ds_(:))));
    dSet2 = [dSet2 ds_];
end
toc;

[U,SVs,V] = svd(dSet2,0); % i should try on a per-component basis, but too late now perhaps for version 3??





