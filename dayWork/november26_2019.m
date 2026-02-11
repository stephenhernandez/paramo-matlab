clear; close all;
tic;
cd ~/research/now/sangay/
load maxccp
toc;
load PUYOWaveformsFromSangay.mat
N = length(S);
lfc = 0.6;
hfc = 1.2;
npoles = 6;
newFs = 10;
tw = 0.004;

S = differentiateWaveforms(S(1:N));
S = taperWaveforms(detrendWaveforms(S),tw);
S = filterWaveforms(S,lfc,hfc,npoles);
S = taperWaveforms(detrendWaveforms(S),tw);
puyo = pull(S);

%%
SAGAFlag = true;
method = 'complete';
thresh = 0.7;

if SAGAFlag
    load waveformDetector_S.mat
    S = differentiateWaveforms(S(1:N));
    S = taperWaveforms(detrendWaveforms(S),tw);
    S = filterWaveforms(S,lfc,hfc,npoles);
    S = taperWaveforms(detrendWaveforms(S),tw);
    S = resampleWaveforms(S,newFs);
    indiv_events = pull(S);
else
    indiv_events = puyo;
    clear puyo;
end

indiv_events = pull(S);
maxAmpRMS = rms(indiv_events)';
clear S;

load waveformDetector
tI = snr >= 100;
tabs = tabs(tI);
toc;
ttmp = datenum(tabs);
winlen = size(indiv_events,1);
indiv_events = normalizeWaveforms(indiv_events);

%% convert cc vectors to squareform (2D) matrix
maxccp1 = squareform(maxccp1);
maxccn1 = squareform(maxccn1);
plags1 = squareform(plags1);
nlags1 = squareform(nlags1);

maxccp1 = 1-maxccp1; %turn similarity matrix into distance matrix
dI = maxccp1 > 0.99999;
maxccp1(dI) = 0;

clear dI
if ~SAGAFlag
    plags1 = 10*plags1;
end

%%

markersize = 12;
%[family,l_uniq_indices] = returnClusts(maxccp1,thresh,method);
load('progressiveClusteringOnSAGAWaveforms_8_75_7');
family = newFamilies;
%l_uniq_indices =
ClustTot = length(newFamilies); %(l_uniq_indices);
%mostMembers = l_uniq_indices(1);
minNumElements = 10;
nFam = 1e3; %find(l_uniq_indices < minNumElements,1) - 1;
%Nsingletons = ClustTot - nFam;

l_families = 0;
maxNFam = 1999;
nFam = min([nFam maxNFam]);
linear_stack = puyo(:,1:nFam);


%%
close all;
% figh(8) = figure('units','normalized','outerposition',[0 0 1 1]);
% plot(1:nFam,l_uniq_indices(1:nFam),'.','markersize',12)
% xlabel('Family Number')
% title('Family Members per Family (min. of 2 members)');
% text(0.8*nFam,0.95*mostMembers,['Number of Singletons: ',num2str(Nsingletons)]);
% zoom on;
MAXN = 200;
winlen = size(puyo,1);
for i = 1:15%nFam
    disp(i)
    family1 = sort(family{i});
    lfam = length(family1);
    maxNumberToActuallyStack = min([lfam,MAXN]);
    
    p = randperm(lfam,maxNumberToActuallyStack);
    family1 = family1(p);
    lfam = length(family1);
    
    indiv_tmp = indiv_events(:,family1);
    indiv_tmp = resample(double(indiv_tmp),10,1);
    puyoTmp = puyo(:,family1);
    [maxccptmp,plagstmp] = doccFreqCircShift(double(indiv_tmp));
    ccData = [maxccptmp plagstmp];
    disp(' ');
%     
%     maxccptmp = maxccp1(:,family1);
%     maxccptmp = maxccptmp(family1,:);
%     maxccptmp = squareform(maxccptmp)'; %column vector
%     
%     plagstmp = plags1(:,family1);
%     plagstmp = plagstmp(family1,:);
%     plagstmp = squareform(plagstmp)';   %column vector
%     
%     ccData = [maxccptmp plagstmp]; % maxccntmp nlagstmp];
%     disp(['ccData: ',num2str(size(ccData))])
    
    
    tic;
    disp(' ')
    disp(['stacking family: ',num2str(i)])
    shifted_data = apply_vdcc(normalizeWaveforms(puyoTmp),ccData);
    toc;
    linear_stack(:,i) = normalizeWaveforms(nanmedian(shifted_data,2));
    if lfam >= minNumElements
        l_families = l_families + 1;
        fig = figure(1000+i);
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1 1];
        if SAGAFlag
            Fs = newFs;
        else
            Fs = 100;
        end
        linear_stack(:,i) = plot_family_and_amps(shifted_data,tabs(family1),maxAmpRMS(family1),1:lfam,5,Fs); %,-inf,-inf,false,[refTime refTime+length(days)],[plotMinAmp plotMaxAmp]);
        zoom on;
        
        figure(1000+i);
        figure(10000+i);
        close(10000+i);
        
        fig = figure(10000+i);
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1 1];
        
        ha(1) = subplot(211);
        plot(tabs(family1),maxAmpRMS(family1),'o','color',[0.5 0.5 0.5],'markersize',markersize,'linewidth',2); hold on;
        [P_,S_,MU_] = polyfit(ttmp(family1),log10(maxAmpRMS(family1)),1);
        plot(tabs(family1),10.^(polyval(P_,ttmp(family1),S_,MU_)),'k-','linewidth',2);
        title(['Family ',num2str(i),', N = ',num2str(length(family1))])
        
        ha(2) = subplot(212);
        tdum = (0:winlen-1)/(Fs*10);
        plot(tdum,normalizeWaveforms(shifted_data,1,0),'-','color',[0.5 0.5 0.5],'linewidth',1);
        hold(ha(2),'on');
        plot(tdum,normalizeWaveforms(linear_stack(:,i),1,0),'k','linewidth',2);
        xlabel('time [sec.]')
        ylabel('Normalized Amplitude');
        axis(ha(2),'tight');
        zoom on;
        
        figure(10000+i);
        pause(1);
        
        close(1000+i);
        close(10000+i);
    end
end