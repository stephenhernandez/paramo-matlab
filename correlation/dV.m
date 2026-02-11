function [dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,...
    dVsymStack,ccCaus,ccAcaus,ccSymmetric,allFlag,caus,acaus,symStack,...
    maxPrcnt,rmsPrcntCaus,rmsPrcntAcaus,rmsPrcntSym] = dV(varargin)
%
% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2018,06,20),...        % tStart
    datetime(2018,06,30),...        % tEnd
    datetime(2018,06,25),...        % referenceTime
    ["SN06","HHZ","9D",""],...      % SNCL1
    ["SN02","HHZ","9D",""],...      % SNCL2
    0.5,...                         % dW
    0.2,...                         % lfc
    1,...                           % hfc
    4,...                           % cwiStart
    12,...                          % cwiDur
    false,...                       % allFlag
    32,...                          % newFs
    2^10,...                        % secDur
    128,...                         % maxLag
    false,...                       % dailyFlag
    '~/rawdata/'};                  % pathToWaveformServerOrDirectory

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;

[tStart,tEnd,referenceStartTime,SNCL1,SNCL2,dW,lfc,hfc,cwiStart,cwiDur,...
    allFlag,newFs,secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory] = deal(optsToUse{:});

%%
%cd ~/research/now/turrialba_irazu/VTCE/;
stackMethod = 'pws';
technique = 'stretch';
force32bit = false;

%%
stnm1 = SNCL1(1);
stnm2 = SNCL2(1);
chan1 = SNCL1(2);
chan2 = SNCL2(2);
cwiEnd = cwiStart + cwiDur;
maxPrcnt = 5;
SSAACC = sum(SNCL1 == SNCL2) == 4; % Single _Station _Auto _And _Cross _Correlation

%%
charchan1 = char(chan1);
charchan1 = charchan1(1:2);
charchan2 = char(chan2);
charchan2 = charchan2(1:2);

%%
if SSAACC
    if ~allFlag
        % just get the single chans over with
        chans1 = chan1;
        chans2 = chan2;
        
        lchans = length(chans1);
        for i = 1:lchans
            chan1 = chans1(i);
            chan2 = chans2(i);
            SNCL1(2) = chan1;
            SNCL2(2) = chan2;
            
            tic;
            [caus,acaus,symStack,t,lags] = getDV(tStart,tEnd,SNCL1,SNCL2,dW,lfc,hfc,newFs,secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory);
            [dVcaus,dVacaus,dVsymmetric,ccCaus,ccAcaus,ccSymmetric,...
                rmsPrcntCaus,rmsPrcntAcaus,rmsPrcntSym] ...
                = pii(t,caus,acaus,symStack,referenceStartTime,cwiStart,...
                cwiEnd,newFs,stnm1,stnm2,chan1,chan2,maxPrcnt,stackMethod,technique,force32bit,lfc,hfc);
            toc;
        end
    else
        % for SSAACC, for each day read all three components
        chans1 = [string([charchan1,'Z']); string([charchan1,'N']); string([charchan1,'E']); ...
            string([charchan1,'Z']); string([charchan1,'Z']); string([charchan1,'E'])];
        chans2 = [string([charchan2,'Z']); string([charchan2,'N']); string([charchan2,'E']); ...
            string([charchan1,'E']); string([charchan2,'N']); string([charchan2,'N'])];
        
        dVcaus = [];
        dVacaus = [];
        dVsymmetric = [];
        
        ccCaus = [];
        ccAcaus = [];
        ccSymmetric = [];
        
        rmsPrcntCaus = [];
        rmsPrcntAcaus = [];
        rmsPrcntSym = [];
        
        lchans = length(chans1);
        for i = 1:lchans
            chan1 = chans1(i);
            chan2 = chans2(i);
            SNCL1(2) = chan1;
            SNCL2(2) = chan2;
            
            % get correlation functions
            [caus,acaus,symStack,t,lags] = getDV(tStart,tEnd,SNCL1,SNCL2,dW,...
                lfc,hfc,newFs,secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory);
            
            % extract velocity changes
            [dVcaus_,dVacaus_,dVsymmetric_,ccCaus_,ccAcaus_,ccSymmetric_,...
                rmsPrcntCaus_,rmsPrcntAcaus_,rmsPrcntSym_] ...
                = pii(t,caus,acaus,symStack,referenceStartTime,cwiStart,...
                cwiEnd,newFs,stnm1,stnm2,chan1,chan2,maxPrcnt,stackMethod,technique,force32bit,lfc,hfc);
            
            % velocity change
            dVcaus = [dVcaus dVcaus_];
            dVacaus = [dVacaus dVacaus_];
            dVsymmetric = [dVsymmetric dVsymmetric_];
            
            % peak cc
            ccCaus = [ccCaus ccCaus_];
            ccAcaus = [ccAcaus ccAcaus_];
            ccSymmetric = [ccSymmetric ccSymmetric_];
            
            % associated error
            rmsPrcntCaus = [rmsPrcntCaus rmsPrcntCaus_];
            rmsPrcntAcaus = [rmsPrcntAcaus rmsPrcntAcaus_];
            rmsPrcntSym = [rmsPrcntSym rmsPrcntSym_];
        end
    end
else
    if allFlag
        % do not use allFlag=true for sensors from two different locations
        charchan1 = char(chan1);
        orient1 = charchan1(3);
        charchan1 = charchan1(1:2);
        
        charchan2 = char(chan2);
        orient2 = charchan2(3);
        charchan2 = charchan2(1:2);
        
        
        chans1 = [string([charchan1,'Z']); string([charchan1,'N']); string([charchan1,'E']); ...
            string([charchan1,'Z']); string([charchan1,'Z']); string([charchan1,'E'])];
        chans2 = [string([charchan2,'Z']); string([charchan2,'N']); string([charchan2,'E']); ...
            string([charchan1,'E']); string([charchan2,'N']); string([charchan2,'N'])];
        
        dVcaus = [];
        dVacaus = [];
        dVsymmetric = [];
        ccCaus = [];
        ccAcaus = [];
        ccSymmetric = [];
        rmsPrcntCaus = [];
        rmsPrcntAcaus = [];
        rmsPrcntSym = [];
        
        lchans = length(chans1);
        for i = 1:lchans
            chan1 = chans1(i);
            chan2 = chans2(i);
            SNCL1(2) = chan1;
            SNCL2(2) = chan2;
            
            tic;
            [dayStack,t,lags] = greenFromNoise(tStart,tEnd,SNCL1,SNCL2,dW,lfc,hfc,newFs,secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory);
            toc;
            [caus,acaus,symStack] = decomposeNoiseCorrMatrix(dayStack);
            toc;
            [dVcaus_,dVacaus_,dVsymmetric_,ccCaus_,ccAcaus_,ccSymmetric_] = ...
                getDV(t,caus,acaus,symStack,referenceStartTime,cwiStart,cwiEnd,newFs,stnm1,stnm2,chan1,chan2,maxPrcnt);
            toc;
            
            dVcaus = [dVcaus dVcaus_];
            dVacaus = [dVacaus dVacaus_];
            dVsymmetric = [dVsymmetric dVsymmetric_];
            ccCaus = [ccCaus ccCaus_];
            ccAcaus = [ccAcaus ccAcaus_];
            ccSymmetric = [ccSymmetric ccSymmetric_];
        end
    else
        % just get the single chans over with
        chans1 = chan1;
        chans2 = chan2;
        
        lchans = length(chans1);
        for i = 1:lchans
            chan1 = chans1(i);
            chan2 = chans2(i);
            SNCL1(2) = chan1;
            SNCL2(2) = chan2;
            
            tic;
            [caus,acaus,symStack,t,lags] = getDV(tStart,tEnd,SNCL1,SNCL2,dW,lfc,hfc,newFs,secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory);
            [dVcaus,dVacaus,dVsymmetric,ccCaus,ccAcaus,ccSymmetric,...
                rmsPrcntCaus,rmsPrcntAcaus,rmsPrcntSym] ...
                = pii(t,caus,acaus,symStack,referenceStartTime,cwiStart,...
                cwiEnd,newFs,stnm1,stnm2,chan1,chan2,maxPrcnt,stackMethod,technique,force32bit,lfc,hfc);
            toc;
        end
    end
end

%%
figure();
ax(1) = subplot(311);
plot(t,dVcaus); zoom on;
ax(2) = subplot(312);
plot(t,dVacaus); zoom on;
ax(3) = subplot(313);
plot(t,dVsymmetric); zoom on;
linkaxes(ax,'x');

figure();
imagesc(sign(normalizeWaveforms(caus)'));
colorbar; zoom on;

% dVstack = nanmedian([nanmedian(dVcaus,2) nanmedian(dVacaus,2) nanmedian(dVsymmetric,2)],2);
% dayStack = normalizeWaveforms(dayStack);
% dayStack = normalizeWaveforms(dVstack);
%
% refTrace = normalizeWaveforms(nanmedian(dayStack,2));
% refTrace = repmat(refTrace,1,size(dayStack,2));
% ccMaster = single(sum(refTrace .* dayStack)');
% refTrace = single(refTrace(:,1));

%%
%dVstack = nanmedian([nanmedian(dVcaus,2) nanmedian(dVacaus,2) nanmedian(dVsymmetric,2)],2);
lchans = size(dVcaus,2);
if lchans > 3
    %dVstack = nanmean([dVcaus dVacaus(:,4:end)],2); % 9 measurements % dVsymmetric(:,4:end)],2); %12 measurements here
    
    %     ccSum = sum([ccCaus ccAcaus(:,4:end)],2);
    %     dVstack = sum([dVcaus dVacaus(:,4:end)].*([ccCaus ccAcaus(:,4:end)]./ccSum),2); %weighted sum based on cc value
    
    ccMaster = [ccCaus ccAcaus(:,4:end) ccSymmetric(:,4:end)];
    ccSum = sum(ccMaster,2);
    ccWeights = (ccMaster./ccSum);
    dVstack = sum([dVcaus dVacaus(:,4:end) dVsymmetric(:,4:end)].*ccWeights,2); %weighted sum based on cc value
else
    dVstack = nanmean(dVcaus,2);
end

ccSum = sum(ccSymmetric,2);
disp(max(ccSum))
dVsymStack = nanmedian(dVsymmetric,2); %sum(dVsymmetric.*(ccSymmetric./ccSum),2); %weighted sum based on cc value

%dVsymStack2 = nanmean(dVsymmetric,2);
%maxPrcnt = 3;

maxWindows = floor(86400/secDur);
Nmed = max([4 2^(nextpow2(maxWindows)-6)]);
fprintf('Nmed: %d\n',Nmed);
plotFlag = true;
%ccGood = ccMaster >= 0.2 &
ccGood = abs(dVsymStack) <= maxPrcnt+1 & abs(dVstack) <= maxPrcnt+1;
if plotFlag
    close all;
    minCCLim = 0.25;
    maxCCLim = 0.75;
    
    % figure 1
    figure('units','normalized','outerposition',[0 0 1 1]);
    stack1 = subplot(211);
    plot(stack1,t(ccGood),dVstack(ccGood),'.'); zoom on; title('$dV_{stack}$');
    hold on;
    stackMed = medfiltSH(dVstack,Nmed,true);
    stackMed(ccGood) = medfiltSH(dVstack(ccGood),Nmed,true);
    stackMed(~ccGood) = NaN;
    plot(stack1,t,stackMed,'linewidth',2);
    ylabel('$dv/v[\%]$');
    plot(stack1,t,dVstack+rmsPrcntSym,'k','linewidth',0.5);
    plot(stack1,t,dVstack-rmsPrcntSym,'k','linewidth',0.5);
    plot(stack1,[min(t) max(t)],[0 0],'--','Color',[0.5 0.5 0.5]);
    
    stack2 = subplot(212);
    plot(stack2,t(ccGood),dVsymStack(ccGood),'.'); zoom on; title('$dV_{symStackOnly}$');
    hold on;
    symMed = medfiltSH(dVsymStack,Nmed,true);
    symMed(ccGood) = medfiltSH(dVsymStack(ccGood),Nmed,true);
    symMed(~ccGood) = NaN;
    plot(stack2,t,symMed,'linewidth',2);
    ylabel('$dv/v[\%]$');
    linkaxes([stack1 stack2],'x');
    plot(stack2,t,dVsymStack+nanmean(rmsPrcntSym,2),'k','linewidth',0.5);
    plot(stack2,t,dVsymStack-nanmean(rmsPrcntSym,2),'k','linewidth',0.5);
    plot(stack2,[min(t) max(t)],[0 0],'--','Color',[0.5 0.5 0.5]);
    
    if allFlag
        lchans = size(dVcaus,2);
        axCaus = gobjects(lchans,1);
        axAcaus = axCaus;
        axSym = axCaus;
        
        titles = ["ZZ" "NN" "EE" "ZE" "ZN" "EN"];
        
        % figure 2
        figure('units','normalized','outerposition',[0 0 1 1]);
        for f = 1:lchans
            axCaus(f) = subplot(lchans/3,3,f);
            scatter(t(ccGood),dVcaus(ccGood,f),[],ccCaus(ccGood,f),'.');
            zoom on; title(titles(f));
            ylabel('$dv/v[\%]$'); colorbar; caxis([minCCLim maxCCLim]); ylim([-1 1]*maxPrcnt);
        end
        linkaxes(axCaus,'xy');
        figure(2); suptitle('$dV_{causal}$');
        colormap(flipud(winter));
        
        % figure 3
        figure('units','normalized','outerposition',[0 0 1 1]);
        for f = 1:lchans
            axAcaus(f) = subplot(lchans/3,3,f);
            scatter(t(ccGood),dVacaus(ccGood,f),[],ccAcaus(ccGood,f),'.');
            zoom on; title(titles(f));
            ylabel('$dv/v[\%]$'); colorbar; caxis([minCCLim maxCCLim]); ylim([-1 1]*maxPrcnt);
        end
        linkaxes(axAcaus,'xy');
        figure(3); suptitle('$dV_{acausal}$');
        colormap(flipud(winter));
        
        % figure 4
        figure('units','normalized','outerposition',[0 0 1 1]);
        for f = 1:lchans
            axSym(f) = subplot(lchans/3,3,f);
            scatter(t(ccGood),dVsymmetric(ccGood,f),[],ccSymmetric(ccGood,f),'.');
            zoom on; title(titles(f));
            ylabel('$dv/v[\%]$'); colorbar; caxis([minCCLim maxCCLim]); ylim([-1 1]*maxPrcnt);
        end
        linkaxes(axSym,'xy');
        figure(4);
        suptitle('$dV_{symmetric}$');
        colormap(flipud(winter));
    end
    
    % figure 5
    figure('units','normalized','outerposition',[0 0 1 1]);
    signax(1) = subplot(2,1,1);
    imagesc(signax(1),lags(lags >= 0),(1:length(t))',sign(caus'));
    colorbar;
    zoom on;
    title('causal');
    
    signax(2) = subplot(2,1,2);
    imagesc(signax(2),lags(lags >= 0),(1:length(t))',sign(acaus'));
    colorbar;
    zoom on;
    title('acausal');
    linkaxes(signax,'x');
    
    % figure 6
    figure('units','normalized','outerposition',[0 0 1 1]);
    stack1 = subplot(211);
    stackMed = medfiltSH(dVstack,Nmed,true);
    stackMed(ccGood) = medfiltSH(dVstack(ccGood),Nmed,true);
    stackMed(~ccGood) = NaN;
    plot(stack1,t,dVstack-stackMed,'.');
    zoom on; title('noise $dV_{stack}$');
    ylabel('$dv/v[\%]$');
    
    stack2 = subplot(212);
    symMed = medfiltSH(dVsymStack,Nmed,true);
    symMed(ccGood) = medfiltSH(dVsymStack(ccGood),Nmed,true);
    symMed(~ccGood) = NaN;
    plot(stack2,t,dVsymStack-symMed,'.');
    ylabel('$dv/v[\%]$');
    linkaxes([stack1 stack2],'x');
    
    % figure 7
    figure();
    if lchans > 3
        plot(t,ccWeights,'.'); zoom on; %title('mean weight based on cc values');
    else
        plot(t,[ccCaus ccAcaus(:,4:end) ccSymmetric(:,4:end)],'.'); zoom on;
    end
end
