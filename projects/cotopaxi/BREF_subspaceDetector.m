%function BREF_subspaceDetector()
clear; close all; clc;

%%
tStart = datetime(2006,01,01);
tEnd = datetime(2020,08,26);
days = (tStart:tEnd)';

%%
maxN = 1e7;
lieOrig = maxN*1e2;
tabs = NaT(maxN,1);
scaledCC = NaN(maxN,1);
z2p = scaledCC;
p2rms = scaledCC;
NCC = scaledCC;
Neff = scaledCC;
maxBasisFunctions = 20;
secDuration = 60; %value for BREF
thresh = 3.5;
maxN = 5e2;
mpd = 5;

%%
cd ~/research/now/cotopaxi/
basisFunctionFileName = '~/research/now/cotopaxi/bref_svd_basis_functions';

%%
kstnms = "BREF";
maxSensors = length(kstnms);
z2p = repmat(z2p,1,maxSensors);

%%F
pullWaveformsFlag = false;
if pullWaveformsFlag
    indiv_events = NaN(lieOrig,1);
end

%%
ntot = 1;
for i = 1:length(days)
    S_ = loadWaveforms(days(i),1,kstnms,"BHZ");
    refs = pull(S_,'ref');
    gI = ~isnat(refs);
    if sum(gI)
        S_ = S_(gI);
        if pullWaveformsFlag
            [indiv_events_,tabs_,scaledCC_,z2p_,NCC_,Neff_,p2rms_,kurt_,ccnorm,t] = ...
                subspaceDetector(S_,thresh,basisFunctionFileName,maxBasisFunctions,secDuration,maxN,mpd);
        else
            [~,tabs_,scaledCC_,z2p_,NCC_,Neff_,p2rms_,kurt_,ccnorm,t] = ...
                subspaceDetector(S_,thresh,basisFunctionFileName,maxBasisFunctions,secDuration,maxN,mpd);
        end
        lpks = length(scaledCC_);
        if lpks
            if pullWaveformsFlag
                if i == 1
                    winlen = size(indiv_events_,1);
                    indiv_events = reshape(indiv_events(1:(winlen*floor(lieOrig/winlen))),[winlen floor(lieOrig/winlen)]);
                end
            end
            
            %%
            if pullWaveformsFlag
                indiv_events(1:winlen,ntot:ntot+lpks-1) = indiv_events_;
            end
            tabs(ntot:ntot+lpks-1) = tabs_;
            scaledCC(ntot:ntot+lpks-1) = scaledCC_;
            
            nc = size(z2p_,2);
            if nc > maxSensors
                nc = maxSensors;
            end
            z2p(ntot:ntot+lpks-1,1:nc) = z2p_(:,1:nc);
            p2rms(ntot:ntot+lpks-1,1:nc) = p2rms_(:,1:nc);
            kurt(ntot:ntot+lpks-1,1:nc) = kurt_(:,1:nc);
            NCC(ntot:ntot+lpks-1) = NCC_;
            Neff(ntot:ntot+lpks-1) = Neff_;
            ntot = ntot+lpks;
        end
    else
        disp('data not found');
    end
end

%%
if pullWaveformsFlag
    clear indiv_events_
end
clear tabs_ pks_ pksOrig_ z2p_ Neff_ p2rms_ kurt_

%%
ntot = ntot - 1;
if pullWaveformsFlag
    indiv_events = indiv_events(:,1:ntot);
end
tabs = tabs(1:ntot);
scaledCC = scaledCC(1:ntot);
NCC = NCC(1:ntot);
kurt = kurt(1:ntot,:);
z2p = z2p(1:ntot,:);
p2rms = p2rms(1:ntot,:);
Neff = Neff(1:ntot);

%%
discr = z2p.*kurt./p2rms; % kurtosis * rms
rmsEstimate = z2p./p2rms;
discr2 = kurt./rmsEstimate;
goodI = z2p >= 200; % | discr2 < 0.1;

%%
clearvars -except indiv_events tabs scaledCC NCC z2p p2rms Neff ...
    kurt discr discr2 rmsEstimate goodI ccnorm t secDuration S_

%%
plotFlag = 0;
if plotFlag
    %%
    close all;
    pI = true(size(kurt)); 
    
    %%
    figure('units','normalized','outerposition',[0 0 0.8 1]);
    ax(1) = subplot(611);
    semilogy(tabs,nanmedian(z2p,2),'.'); zoom on; grid on;
    title('zero-to-peak');
    ax(2) = subplot(612);
    plot(tabs,nanmedian(p2rms,2),'.'); zoom on; grid on;
    title('peak-to-rms');
    ax(3) = subplot(613);
    plot(tabs,nanmedian(NCC,2),'.'); zoom on; grid on;
    title('gamma');
    ax(4) = subplot(614);
    plot(tabs,nanmedian(kurt,2),'.'); zoom on; grid on;
    title('kurtosis');
    ax(5) = subplot(615);
    plot(tabs,rmsEstimate./kurt,'.'); zoom on; grid on;
    title('rms / kurtosis');
    ax(6) = subplot(616);
    plot(tabs,error,'.'); zoom on; grid on;
    title('error');
    linkaxes(ax,'x');
    
    %%
    figure(); plot(tabs(pI),1:length(tabs(pI)),'.-'); 
    zoom on; grid on;
    
    %%
    figure(); 
    [N,edges] = histcounts(tabs(pI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
    N = N';
    edges = edges(1:end-1)';
    plot(edges,N,'.'); zoom on; grid on;
    
    %%
    figure();
    ss = scatter(1./discr2,z2p,3*exp(log10(z2p)),NCC,'filled');
    ax2 = gca;
    ax2.XScale = 'log';
    ax2.YScale = 'log';
    
    cc = colorbar; zoom on;
    ss.MarkerFaceAlpha = 0.4;
    ss.MarkerEdgeColor = 'k';
    ss.MarkerEdgeAlpha = 0.2;
    grid on;
    
    %%
    dI = ~goodI;
    figure(1);
    hold(ax(1),'on'); semilogy(ax(1),tabs(dI),nanmedian(z2p(dI),2),'p'); zoom on; grid on;
    
    figure(); 
    plot(t,ccnorm); zoom on;
    hold on;
    plot(tabs+seconds(secDuration),NCC,'o'); zoom on;rm
    plot(tabs(dI)+seconds(secDuration),NCC(dI),'p'); zoom on;
    grid on;
    
    figure(); 
    plot(tabs(goodI),1:sum(goodI),'.'); zoom on; grid on;
else
    clear S_ t ccnorm
end

%%
%save('cotopaxiSubspaceDetectorBREF','-v7.3');
