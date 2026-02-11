clear; close all; clc;
tStart = datetime(2020,07,01);
tEnd = datetime(2021,03,31);
%tStart = datetime(2020,10,01);
%tEnd = datetime(2021,03,15);
days = (tStart:tEnd)';
ldays = length(days);
maxN = 300;

%%
tabs = NaT(maxN,1);
NCC = NaN(maxN,1);
z2p = NaN(maxN,3);
Neff = NCC;
p2rms = z2p;
kurt = z2p;

%%
plotFlag = 0;
verboseFlag = true;
threshold = 0.2;

%%
n = 1;
for i = 1:ldays
    day_ = days(i);
    
    %%
    %S_ = loadWaveforms(day_,1,"OTAV",["BHZ";"BH1";"BH2"],"IU","10");
    S_ = loadWaveforms(day_,1,"RVRD",["HHZ";"HHN";"HHE"],"EC","");
    
    %%
    if ~isnat(S_(1).ref)
        %[~,tabs_,NCC_,z2p_,Neff_,p2rms_,kurt_] = subspaceDetector(S_,...
        %    threshold,'~/research/now/esmeraldas/OTAV_3cmp_basisFunctions',5,55,100,20,...
        %    false,true,plotFlag,verboseFlag);
        
    %     [indiv_events,tabs,energyRatio,z2p,Neff,p2rms,kurt,eratio,t] = ...
    % subspaceDetector(S,threshold,basisFunctionFileName,maxBasisFunctions,...
    % recordLength,maxN,mpd,diffFlag,linearccnorm,plotFlag,verboseFlag,smoothFlag);

        % [~,tabs_,NCC_,z2p_,Neff_,p2rms_,kurt_] = subspaceDetector(S_,...
        %     threshold,'~/research/now/esmeraldas/RVRD_3cmp_basisFunctions',3,35,100,10,...
        %     false,true,plotFlag,verboseFlag);
        [~,tabs_,NCC_,z2p_,Neff_,p2rms_,kurt_] = subspaceDetector(S_,...
            threshold,'~/research/now/esmeraldas/RVRD_3cmp_basisFunctions',3,35,100,10,...
            false,true,plotFlag,verboseFlag);
        
        if ~isempty(tabs_)
            lnew = length(tabs_);
            tabs(n:n+lnew-1) = tabs_;
            NCC(n:n+lnew-1) = NCC_;
            z2p(n:n+lnew-1,:) = z2p_;
            Neff(n:n+lnew-1) = Neff_;
            p2rms(n:n+lnew-1,:) = p2rms_;
            kurt(n:n+lnew-1,:) = kurt_;
            n = n + lnew;
        end
    end
end
gI = ~isnat(tabs);
tabs = tabs(gI);
NCC = NCC(gI);
z2p = z2p(gI);
Neff = Neff(gI);
p2rms = p2rms(gI);
kurt = kurt(gI);

%%
%R = extractWaveforms(tabs,seconds(55),"OTAV","BHZ","IU","10",true,true);
%Rf = detrendWaveforms(R);
%Rf = intWaveforms(detrendWaveforms(filterWaveforms(differentiateWaveforms(Rf),2,4,2,0.01,false)));
lfc = 1;
hfc = 4;
newFs = 50;

if ~exist('R','var')
    tic;
    R = extractWaveforms(tabs,seconds(35),"RVRD","HHZ","EC","",true,true);
    Rf = detrendWaveforms(R);
    Rf = intWaveforms(detrendWaveforms(filterWaveforms(differentiateWaveforms(Rf),lfc,hfc,4,0.02,true)));
    Rf = resampleWaveforms(Rf,newFs);
    pullRf = normalizeWaveforms(detrend(double(pull(Rf))));
    
    R = extractWaveforms(tabs,seconds(35),"RVRD","HHN","EC","",true,true);
    Rf = detrendWaveforms(R);
    Rf = intWaveforms(detrendWaveforms(filterWaveforms(differentiateWaveforms(Rf),lfc,hfc,4,0.02,true)));
    Rf = resampleWaveforms(Rf,newFs);
    pullRf = [pullRf; normalizeWaveforms(detrend(double(pull(Rf))))];
    
    R = extractWaveforms(tabs,seconds(35),"RVRD","HHE","EC","",true,true);
    Rf = detrendWaveforms(R);
    Rf = intWaveforms(detrendWaveforms(filterWaveforms(differentiateWaveforms(Rf),lfc,hfc,4,0.02,true)));
    Rf = resampleWaveforms(Rf,newFs);
    pullRf = [pullRf; normalizeWaveforms(detrend(double(pull(Rf))))];
    toc;
% else
%     Rf = detrendWaveforms(R);
%     Rf = intWaveforms(detrendWaveforms(filterWaveforms(differentiateWaveforms(Rf),lfc,hfc,4,0.02,true)));
%     Rf = resampleWaveforms(Rf,newFs);
%     pullRf = detrend(double(pull(Rf)));
end

%%
close all;
[maxccp,plags] = doccFreqCircShift(pullRf);
[family,l_uniq_indices,Nsingletons,tree] = generateFamilies(maxccp,0.6,'weighted',2,false);
logamp = log10(max(z2p,[],2)) - 1.;
for i = 1:length(family)
    if length(family{i}) > 1
        figure(1);
        ax(1) = subplot(211); hold on;
        pp = plot(tabs(family{i}),1:length(family{i}),'o-','linewidth',6);
        pp.Color(4) = 0.6; zoom on;
        
        ax(2) = subplot(212); hold on;
        pp = plot(tabs(family{i}),logamp(family{i}),'o','linewidth',6);
        pp.Color(4) = 0.6; zoom on;
        
        data = detrend(apply_vdcc(detrend(pullRf(:,family{i}))));
        tmp_stack = plot_family(data,(1:size(data,2))',8,newFs);
    end
end
linkaxes(ax,'x');
%[CN,I,families,lf] = plotClusterTimeSeries(family,tabs,10.^logamp);

%
%figure(); plot(tabs,(1:length(tabs))','.'); zoom on;

legend(ax(1),{'Family 1','Family 2','Family 3','Family 4','Family 5','Family 6', 'Family 7'},'Location','NorthWest');
ylabel(ax(1),'Numero Acumulativo');
ylabel(ax(2),'Magnitud Estimada');

