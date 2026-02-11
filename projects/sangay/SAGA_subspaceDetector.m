%function PUYO_subspaceDetector()
clear; close all; clc;

%%
tStart = datetime(2020,08,01);
tEnd = datetime(2020,08,02);
days = (tStart:tEnd)';

%%
maxN = 1e7;
lieOrig = maxN*1e2;
tabs = NaT(maxN,1);
pks = NaN(maxN,1);
z2p = pks;
p2rms = pks;
pksOrig = pks;
Neff = pks;
maxBasisFunctions = 10;
thresh = 0.1;

%%
cd ~/research/now/sangay/
%basisFunctionFileName = '~/research/now/sangay/puyo_svd_basis_functions';
%basisFunctionFileName = '~/research/now/sangay/sangay_svd_basis_functions';
%basisFunctionFileName = '~/research/now/sangay/saga_svd_basis_functions_2';
basisFunctionFileName = '~/research/now/sangay/saga_svd_basis_functions_4';

%%
kstnms = "SAGA";
maxSensors = length(kstnms);
z2p = repmat(z2p,1,maxSensors);


%%
pullWaveformsFlag = false;
if pullWaveformsFlag
    indiv_events = NaN(lieOrig,1);
end

%%
ntot = 1;
for i = 1:length(days)
    S_ = loadWaveforms(days(i),1,kstnms,"HHZ");
    refs = pull(S_,'ref');
    gI = ~isnat(refs);
    if sum(gI)
        S_ = S_(gI);
        if pullWaveformsFlag
            [indiv_events_,tabs_,pks_,z2p_,pksOrig_,Neff_,p2rms_,kurt_,ccnorm,t] = ...
                subspaceDetector(S_,thresh,basisFunctionFileName,maxBasisFunctions,180);
        else
            [~,tabs_,pks_,z2p_,pksOrig_,Neff_,p2rms_,kurt_,ccnorm,t] = ...
                subspaceDetector(S_,thresh,basisFunctionFileName,maxBasisFunctions,180);
        end
        lpks = length(pks_);
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
            pks(ntot:ntot+lpks-1) = pks_;
            nc = size(z2p_,2);
            if nc > maxSensors
                nc = maxSensors;
            end
            z2p(ntot:ntot+lpks-1,1:nc) = z2p_(:,1:nc);
            p2rms(ntot:ntot+lpks-1,1:nc) = p2rms_(:,1:nc);
            kurt(ntot:ntot+lpks-1,1:nc) = kurt_(:,1:nc);
            pksOrig(ntot:ntot+lpks-1) = pksOrig_;
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
clear tabs_ pks_ pksOrig_ z2p_ Neff_ p2rms_

%%
ntot = ntot - 1;
if pullWaveformsFlag
    indiv_events = indiv_events(:,1:ntot);
end
tabs = tabs(1:ntot);
pks = pks(1:ntot);
pksOrig = pksOrig(1:ntot);
z2p = z2p(1:ntot,:);
p2rms = p2rms(1:ntot,:);
Neff = Neff(1:ntot);

%%
plotFlag = false;
if plotFlag
    close all;
    pI = (kurt >= 9.5 & kurt < 40) | (p2rms >=7 & p2rms < 12); %pksOrig >= 0.2 | (p2rms >= 7 & p2rms < 12);
    figure('units','normalized','outerposition',[0 0 1 1]);
    ax(1) = subplot(411);
    semilogy(tabs,nanmedian(z2p,2),'.'); zoom on;
    title('zero-to-peak');
    ax(2) = subplot(412);
    plot(tabs,nanmedian(p2rms,2),'.'); zoom on;
    title('peak-to-rms');
    ax(3) = subplot(413);
    plot(tabs,nanmedian(pksOrig,2),'.'); zoom on;
    title('gamma');
    ax(4) = subplot(414);
    plot(tabs,nanmedian(kurt,2),'.'); zoom on;
    title('kurtosis');
    linkaxes(ax,'x');
    
    %%
    figure(); plot(tabs(pI),1:length(tabs(pI)),'.-'); zoom on;
    [N,edges] = histcounts(tabs(pI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
    N = N';
    edges = edges(1:end-1)';
    figure(); plot(edges,N,'.'); zoom on;
end
