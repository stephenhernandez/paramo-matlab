%function PUYO_subspaceDetector()
clear; close all; clc;

%%
tStart = datetime(2020,08,11);
tEnd = datetime(2020,08,11);
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
maxBasisFunctions = 20;
thresh = 0.25;

%%
cd ~/research/now/sangay/
%basisFunctionFileName = '/home/shernandez/research/now/sangay/puyo_svd_basis_functions';
basisFunctionFileName = '/home/shernandez/research/now/sangay/sangay_svd_basis_functions';

%%
kstnms = "BMAS"; %;"BULB";"TAIS";"TAMH";"BMAS";"BPAT";"BREF";"PKYU";"CHSH"];
maxSensors = length(kstnms);
z2p = repmat(z2p,1,maxSensors);

%%
pullWaveformsFlag = true;
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
plotFlag = true;
if plotFlag
    close all;
    % pI = z2p >= 8 & kurt < 25 & p2rms < 11.5 & ( kurt >= 5 | p2rms >= 6 ); % PUYO, gamma >= 0.45
    % pI = z2p >= 8 & kurt < 25 & p2rms < 11.5 & ( kurt >= 3.5 | p2rms >= 4.5 ); % TAIS, gamma >= 0.45
    % pI = z2p >= 8 & kurt < 25 & p2rms < 11.5 & ( kurt >= 3.5 | p2rms >= 4 ); % PKYU, gamma >= 0.52
    % pI = z2p >= 30 & kurt < 25 & p2rms < 10 & ( kurt >= 4 | p2rms >= 5 ); % TAMH, gamma >= 0.4
    % pI = z2p >= 30 & kurt < 12 & p2rms < 11 & ( kurt >= 4.5 | p2rms >= 5 ); % BULB, gamma >= 0.25
    pI = z2p >= 10 & kurt < 10 & p2rms < 10 & ( kurt >= 5.5 | p2rms >= 5.5 ); % BMAS, gamma >= 0.25
    % pI = z2p >= 30 & kurt < 14 & p2rms < 12 & ( kurt >= 4 | p2rms >= 5.5 ); % BPAT, gamma >= 0.45
    %pI = z2p >= 25 & kurt < 25 & p2rms < 12.2 & (kurt >= 4.5 | p2rms >= 5.5 ); %CHSH, gamma >= 0.3
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    ax(1) = subplot(411);
    semilogy(tabs,nanmedian(z2p,2),'.'); zoom on;
    title('zero-to-peak');
    grid on;
    ax(2) = subplot(412);
    plot(tabs,nanmedian(p2rms,2),'.'); zoom on;
    title('peak-to-rms');
    grid on;
    ax(3) = subplot(413);
    plot(tabs,nanmedian(pksOrig,2),'.'); zoom on;
    title('gamma');
    grid on;
    ax(4) = subplot(414);
    plot(tabs,nanmedian(kurt,2),'.'); zoom on;
    title('kurtosis');
    grid on;
    linkaxes(ax,'x');
    
    %%
    figure(); plot(tabs(pI),1:length(tabs(pI)),'.-'); zoom on;
    [N,edges] = histcounts(tabs(pI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
    N = N';
    edges = edges(1:end-1)';
    figure(); plot(edges,N,'.'); zoom on;
end
